
#define _USE_MATH_DEFINES
#include <math.h>
#define DTOR (M_PI/180.0)
#define RTOD (180.0/M_PI)

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <list>
#include <stdlib.h>
#ifdef USE_BOOST_REGEX
#include <boost/regex.hpp>
#else
#include <regex>
#endif
#ifdef USE_BOOST_STRINGS
#include <boost/algorithm/string.hpp>
#define lower_case boost::algorithm::to_lower
#else
#define lower_case(x) std::transform((x).begin(), (x).end(), (x).begin(), ::tolower);
#endif
#include "okada.h"
#include "get_image_path.h"
// extern "C" {
#include "tmproj.h"
// }

using namespace std;
#ifdef USE_BOOST_REGEX
using namespace boost;
#endif

// Abstract base class for projections.  
// Three are supported, Ian Reilly's "seismologists projection", as used
// by John Beavan, and TM projections, as used by Ian Hamling's software,
// and a null projection for which the spatial units are the same as the 
// model units.

class Projection
{
    public:
    Projection();
    virtual ~Projection();
    virtual Projection *clone() const=0;
    virtual void XY( double lon, double lat, double &x, double &y )=0;
    virtual void LatLon( double x, double y, double &lon, double &lat )=0;
    virtual void SfConv( double x, double y, double &sf, double &cnv );
    virtual bool IsValid()=0;
    virtual bool HasSfConv(){ return false; }
};

Projection::Projection(){}
Projection::~Projection(){}
void Projection::SfConv( double x, double y, double &sf, double &cnv )
{
    sf=1.0;
    cnv=0.0;
}

// The NullProjection is just that - lon=x, lat=y..

class NullProjection : public Projection
{
    public:
    NullProjection() : Projection() {}
    ~NullProjection() {};
    virtual Projection *clone() const{ return new NullProjection; }
    virtual void XY( double lon, double lat, double &x, double &y ){ x=lon; y=lat; }
    virtual void LatLon( double x, double y, double &lon, double &lat ){ lon=x; lat=y; }
    virtual bool IsValid() { return true; }
};

// Comments from John Beavan...
// "Seismologist's projection" as Ian Reilly disparagingly called it.
//
// Choose (lat0, lon0) near the centre of the area of interest.
//
// y = (lat - lat0)*111.2 km
// x = (lon - lon0)/cos(lat0)*111.2 km
//
// I have never carefully investigated how quickly the errors grow with distance from (lat0, lon0) 
// compared to a standard projection.  Usually, the geodetic data are sparse enough and the fault 
// poorly enough known that such differences are unlikely to be important.
//
// For very big earthquakes (such as Sumatra) the half-space approximation breaks down 
// and it is necessary to use  models of dislocations in a sphere.  Not the case for Darfield, fortunately.

class GNSProjection : public Projection
{
public:
    GNSProjection()
        : lat0(0), lon0(0), isvalid(false) {}
    GNSProjection( double lon0, double lat0 ){ SetOrigin(lon0,lat0); }
    virtual ~GNSProjection(){}
    virtual Projection *clone() const;
    virtual void XY( double lon, double lat, double &x, double &y );
    virtual void LatLon( double x, double y, double &lon, double &lat );
    virtual bool IsValid() { return isvalid; }
    void SetOrigin( double lon0, double lat0 )
    {
        this->lon0 = lon0; this->lat0 = lat0; this->clat0=cos(DTOR*lat0); isvalid = true;
    }
private:
    GNSProjection( const GNSProjection &proj )
        : lon0(proj.lon0), lat0(proj.lat0), clat0(proj.clat0), isvalid(proj.isvalid) {}
    void setup(double lon0, double lat0);
    static const double scale;
    double lat0;
    double lon0;
    double clat0;
    bool isvalid;
};

const double GNSProjection::scale = 111200.0;

Projection *GNSProjection::clone() const
{
    return new GNSProjection( *this );
}

void GNSProjection::XY( double lon, double lat, double &x, double &y )
{
    y = (lat - lat0)*scale;
    x = (lon - lon0)*scale*clat0;
};

void GNSProjection::LatLon( double x, double y, double &lon, double &lat )
{
    lat = y/scale + lat0;
    lon = x/(scale * cos(DTOR*lat)) + lon0;
};

class TMProjection : public Projection
{
public:
    TMProjection() : isvalid(false) {}
    TMProjection( double a, double rf, double cm, double sf, double lto, double fe, double fn ) 
    {
        SetParameters( a, rf, cm, sf, lto, fe, fn );
    }
    virtual Projection *clone() const;
    virtual void XY( double lon, double lat, double &x, double &y );
    virtual void LatLon( double x, double y, double &lon, double &lat );
    virtual void SfConv( double x, double y, double &sf, double &cnv );
    virtual bool IsValid() { return isvalid; }
    virtual bool HasSfConv(){ return true; }
    void SetParameters( double a, double rf, double cm, double sf, double lto, double fe, double fn );
private:
    tmprojection tm;
    bool isvalid;
};

void TMProjection::SetParameters( double a, double rf, double cm, double sf, double lto, double fe, double fn )
{
    define_tmprojection( &tm, a, rf, cm*DTOR, sf, lto*DTOR, fe, fn, 1.0 );
    isvalid=true;
}

Projection *TMProjection::clone() const
{
    if( isvalid ) return new TMProjection( tm.a, tm.rf, tm.meridian*RTOD, tm.scalef, tm.orglat*RTOD, tm.falsee, tm.falsen);
    return new TMProjection();
}

void TMProjection::XY( double lon, double lat, double &x, double &y )
{
    geod_tm( &tm, lon*DTOR, lat*DTOR, &x, &y );
}

void TMProjection::LatLon( double x, double y, double &lon, double &lat )
{
    tm_geod( &tm, x, y, &lon, &lat );
    lon *= RTOD;
    lat *= RTOD;
}

void TMProjection::SfConv( double x, double y, double &sf, double &cnv )
{
    tm_sf_conv( &tm, x, y, &sf, &cnv );
    cnv *= RTOD;
}

// 

class FaultSet
{
public:
    FaultSet() : proj(new NullProjection()), factor(1.0), applyConvergence(true) {};
    FaultSet(Projection *proj) : proj(proj ? proj->clone() : new NullProjection()), factor(1.0), applyConvergence(true) {};
    ~FaultSet();
    bool ReadGNSDefinition( istream &str, int nskip = 0 );
    void setFactor( double fctr ){ factor=fctr; }
    void setApplyConvergence( bool apply ){ applyConvergence = apply; }
    bool CalcOkada( double lon, double lat, double *dislocation, double *strain, double *tilt );
    bool AddOkada( double lon, double lat, double *dislocation, double *strain, double *tilt );
    ostream & write( ostream &os, int style=0, bool header=true );
private:
    Projection *proj;
    double factor;
    bool applyConvergence;
    list<SegmentedFault *>faults;
    list<string> names;
};

FaultSet::~FaultSet()
{
    for( list<SegmentedFault *>::const_iterator i = faults.begin(); i != faults.end(); i++ )
    {
        delete( *i );
    }
    if( proj ) delete proj;
}


bool FaultSet::ReadGNSDefinition( istream &str, int nskip )
{
    string buffer;
    bool ok = true;
    bool started = false;
    bool prjcrds = false;
    bool dsclast = false; // Description is last field
    list<string> fields;

    regex re = regex("^(\\w[\\w\\s]+):\\s+(.*?)\\s*$");
    regex fieldre = regex("^(\\w+)\\=(.*?)$");

    int type, fnum;
    double strike=0.0,dip=0.0,rake=0.0,length=0.0,width=0.0;
    double slip=0.0,Uts=0.0,depth=0.0,lat=0.0,lon=0.0,bottom=0.0;
    double dummy=0.0;
    string name("");
    string fieldlist("");
    type = 3;
    fnum = 0;

    while( ok )
    {
        getline(str,buffer);
        if( ! str.good()) break;
        if( nskip > 0 )
        {
            nskip--;
            continue;
        }

        smatch match;
        if( regex_match(buffer,match,re))
        {
            string command(match[1].str());
            lower_case(command);
            stringstream s(match[2].str());
            if(command == "origin" )
            {
                double lon, lat;
                s >> lon >> lat;
                if( s )
                {
                    if( proj ) delete proj;
                    proj=new GNSProjection( lon, lat );
                }
                prjcrds=false;
            }
            if(command == "projection" )
            {
                string projcode;
                if( proj ) delete proj;
                proj=0;
                s >> projcode;
                if( ! s ) continue;
                lower_case(projcode);
                if( projcode == "nztm")
                {
                    proj=new TMProjection(6378137.0,298.257222101,173.0,0.9996,0.0,1600000.0,10000000.0);
                }
                else if( projcode == "utm59" )
                {
                    proj=new TMProjection(6378137.0,298.257222101,171.0,0.9996,0.0,500000.0,10000000.0);
                }
                else if( projcode == "utm60" )
                {
                    proj=new TMProjection(6378137.0,298.257222101,177.0,0.9996,0.0,500000.0,10000000.0);
                }
                else if( projcode == "tm" )
                {
                    double a,rf,cm,sf,lto,fe,fn;
                    s >> a >> rf >> cm >> sf >> lto >> fe >> fn;
                    if( s )
                    {
                        proj=new TMProjection(a, rf, cm,sf,lto,fe,fn);
                    }
                }
                else if( projcode == "none" )
                {
                    proj=new NullProjection();
                }
                if( proj ) prjcrds=true;
            }
            if(command == "vector orientation")
            {
                string vorient;
                s >> vorient;
                lower_case(vorient);
                if( vorient == "geographic" )
                {
                   applyConvergence=false; 
                }
                else if( vorient == "projection" )
                {
                   applyConvergence=true; 
                }
                else
                {
                    cerr << "Vector orientation must be \"geographic\" or \"projection\"\n" << endl;
                    ok=false;
                    continue;
                }
            }
            if(command == "columns" )
            {
                fieldlist=match[2].str();
            }
            if( command == "fault parameter type" )
            {
                s >> type;
            }
            if(command == "additional data" )
            {
                string fielddata;
                while( s >> fielddata )
                {
                    smatch fieldmatch;
                    if( ! regex_match(fielddata,fieldmatch,fieldre))
                    {
                        cerr << "Cannot interpret additional data " << fielddata << endl;
                        ok=false;
                        continue;
                    }
                    string fieldname(fieldmatch[1].str());
                    lower_case(fieldname);
                    stringstream fieldvalstr(fieldmatch[2].str());

                    if( fieldname == "fault_num" ) fieldvalstr >> fnum;
                    else if( fieldname == "fault_type" ) fieldvalstr >> type;
                    else if( fieldname == "strike_deg" ) fieldvalstr >> strike;
                    else if( fieldname == "dip_deg" ) fieldvalstr >> dip;
                    else if( fieldname == "rake_deg" ) fieldvalstr >> rake;
                    else if( fieldname == "length_km" ) fieldvalstr >> length;
                    else if( fieldname == "width_km" ) fieldvalstr >> width;
                    else if( fieldname == "depth_km" ) fieldvalstr >> depth;
                    else if( fieldname == "top_depth_km" ) fieldvalstr >> depth;
                    else if( fieldname == "bottom_depth_km" ) fieldvalstr >> bottom;
                    else if( fieldname == "slip_m" ) fieldvalstr >> slip;
                    else if( fieldname == "opening_m" ) fieldvalstr >> Uts;
                    else if( ! prjcrds && fieldname == "lat_deg" ) fieldvalstr >> lat;
                    else if( ! prjcrds && fieldname == "lon_deg" ) fieldvalstr >> lon;
                    else if( prjcrds && fieldname == "north_m" ) fieldvalstr >> lat;
                    else if( prjcrds && fieldname == "east_m" ) fieldvalstr >> lon;
                    else
                    {
                        cerr << "Invalid field name " << fieldname << " in additional data " << fielddata << endl;
                    }

                }
            }
            continue;
        }

        // Read field header line
        if( buffer.find("slip_m") != string::npos )
        {
            fieldlist=buffer;
            continue;
        }

        if( fieldlist != "" )
        {
            string field;
            stringstream s(fieldlist);
            fields.clear();
            dsclast=false;
            while( s >> field )
            {
                dsclast = field=="fault_description";
                fields.push_back(string(field));
            }
        }

        // Trim comment and leading spaces
        int p = buffer.find('#');
        if( p != string::npos) { buffer.erase(p); }

        // Skip empty lines
        if( buffer.find_first_not_of(' ') == string::npos ) continue;

        // If fields not defined, then set up default values
        if( fields.size() == 0 )
        {
            fields.push_back("strike_deg");
            fields.push_back("dip_deg");
            fields.push_back("rake_deg");
            fields.push_back("length_km");
            fields.push_back("width_km");
            fields.push_back("slip_m");
            fields.push_back("opening_m");
            fields.push_back("depth_km");
            if( prjcrds )
            {
                fields.push_back("east_m");
                fields.push_back("north_m");
            }
            else
            {
                fields.push_back("lat_deg");
                fields.push_back("lon_deg");
            }
        }
        if( prjcrds &&
                find( fields.begin(), fields.end(), "lat_deg" ) != fields.end() &&
                find( fields.begin(), fields.end(), "east_m" ) == fields.end() )
        {
            prjcrds=false;
        }

        stringstream s(buffer);
        fnum += 1;

        for( list<string>::iterator it=fields.begin(); it != fields.end(); it++ )
        {
            if( *it == "fault_num" ) s >> fnum;
            else if( *it == "fault_type" ) s >> type;
            else if( *it == "strike_deg" ) s >> strike;
            else if( *it == "dip_deg" ) s >> dip;
            else if( *it == "rake_deg" ) s >> rake;
            else if( *it == "length_km" ) s >> length;
            else if( *it == "width_km" ) s >> width;
            else if( *it == "depth_km" ) s >> depth;
            else if( *it == "top_depth_km" ) s >> depth;
            else if( *it == "bottom_depth_km" ) s >> bottom;
            else if( *it == "slip_m" ) s >> slip;
            else if( *it == "opening_m" ) s >> Uts;
            else if( ! prjcrds && *it == "lat_deg" ) s >> lat;
            else if( ! prjcrds && *it == "lon_deg" ) s >> lon;
            else if( prjcrds && *it == "north_m" ) s >> lat;
            else if( prjcrds && *it == "east_m" ) s >> lon;
            else if( *it == "fault_x_km") s >> dummy;
            else if( *it == "fault_y_km") s >> dummy;
            else if( *it == "fault_description" )
            {
                if( s && dsclast )
                {
                    getline(s,name);
                    string::size_type pos=name.find_first_not_of(" \t");
                    if( pos != string::npos ) name=name.substr(pos);
                }
                else
                {
                    s >> name;
                }
            }
            else
            {
                cerr << "Error reading fault definition - unrecognized field " << *it << endl;
                ok = false;
                continue;
            }
        }
        if( ! ok ) continue;
        if( !s )
        {
            cerr << "Error reading fault definition\n"
                     << buffer << endl;
            ok = false;
        }

        string::size_type pos = name.find_first_not_of(" \t");
        if( pos == string::npos )
        {
            name = "";
        }
        else
        {
            name = name.substr(pos);
        }

        // Convert to standard units
        length *= 1000.0;
        width *= 1000.0;
        depth *= 1000.0;
        bottom *= 1000.0;
        strike = 90 - strike;
        dip = -dip;

        if( ! proj )
        {
            cerr << setiosflags( ios::fixed ) << setprecision(5);
            cerr << "Setting projection origin to " << lon << " " << lat << endl;
            proj = new GNSProjection( lon, lat );
            prjcrds=false;
        }

        double x, y;
        // If input coordinate are on a projection, then assume order is east,north
        // Otherwise convert lat/lon to x,y
        if( prjcrds ) 
        {
            x=lon;
            y=lat;
        }
        else
        {
            proj->XY(lon,lat,x,y);
        }

        double fs0=0, fs1=length;
        double fd0=0, fd1=width;

        // Calculating the start and end locations of the fault segment in
        // the fault plane defined by the x,y,depth (reference point) and
        // the strike and dip.  The start and end along the strike is defined
        // by fs[0],fs[1], and the start and end down dip are fd[0],fd[1]
        //
        // x, y, depth defines the reference point from which the fault
        // is defined.

        switch( type )
        {
        case 0:
            fs0 = -fs1; fs1 = 0.0;
            fd0 = -fd1; fd1 = 0.0;
            break;
        case 1:
            fs0 = -fs1; fs1 = 0.0;
            fd0 = depth/sin(dip*DTOR);
            fd1 += fd0;
            depth = 0.0;
            break;
        case 2:
            fs0 = -fs1; fs1 = 0.0;
            break;
        case 3:
            fs1 /= 2; fs0 = -fs1;
            fd1 /= 2; fd0 = -fd1;
            break;
        case 4:
            fs1 /= 2; fs0 = -fs1;
            fd0 = depth/sin(dip*DTOR);
            fd1 += fd0;
            depth = 0.0;
            break;
        case 5:
            fs1 /= 2; fs0 = -fs1;
            break;
        case 6:
            fs1 /= 2; fs0 = -fs1;
            dip=dip+180.0;
            fd0 = depth/sin(dip*DTOR);
            fd1 = bottom/sin(dip*DTOR);
            depth = 0.0;
            rake=180.0-rake;
            break;
        case 7:
            fs1 /= 2; fs0 = -fs1;
            dip=dip+180.0;
            depth = (depth+bottom)/2.0;
            fd1 = (bottom-depth)/sin(dip*DTOR);
            fd0 = -fd1;
            rake=180.0-rake;
            break;
            // Next time this is edited look at rewriting as per 
            // github issue 2.
        default:
            cerr << "Invalid fault type " << type << endl;
            ok=false;
            continue;
        }

        double Uss = slip * cos(rake*DTOR);
        double Uds = slip * sin(rake*DTOR);

        double fs[2] = { fs0, fs1 };
        double fd[2] = { fd0, fd1 };

        SegmentedFault *f = new SegmentedFault( x, y, depth, strike, dip );
        f->SetId( fnum );
        f->SetStrikeSegBreaks( 1, fs );
        f->SetDipSegBreaks( 1, fd );
        f->SetFaultSlip( 0, 0, Uss, Uds, Uts );

        faults.push_back( f );
        names.push_back(name);
    }
    if( faults.size() == 0 )
    {
        cerr << "Fault definition file contains no fault definitions" << endl;
        ok = false;
    }
    return ok;
}

bool FaultSet::CalcOkada( double lon, double lat, double *dislocation, double *strain, double *tilt )
{
    dislocation[0] = dislocation[1] = dislocation[2] = 0.0;
    if( strain ) { strain[0] = strain[1] = strain[2] = strain[3] = 0.0; }
    if( tilt ) { tilt[0] = tilt[1] = 0.0; }
    AddOkada( lon, lat, dislocation, strain, tilt );
}

bool FaultSet::AddOkada( double lon, double lat, double *dislocation, double *strain, double *tilt )
{
    double x, y;
    double denu[3];
    bool ok = true;
    proj->XY(lon,lat,x,y);
    for( list<SegmentedFault *>::const_iterator i = faults.begin(); i != faults.end(); i++ )
    {
        if( ! (*i)->AddOkada(x,y,denu,strain,tilt,factor)) ok = false;
    }
    if( dislocation )
    {
        if( applyConvergence && proj->HasSfConv() )
        {
            double sf,conv;
            proj->SfConv(x,y,sf,conv);
            sf=1.0/sf;
            double ccnv=cos(DTOR*conv);
            double scnv=-sin(DTOR*conv);
            double de = (denu[0]*ccnv+denu[1]*scnv)/sf;
            double dn = (denu[1]*ccnv-denu[0]*scnv)/sf;
            denu[0]=de;
            denu[1]=dn;
        }
        dislocation[0] += denu[0];
        dislocation[1] += denu[1];
        dislocation[2] += denu[2];
    }

    return ok;
}

ostream &FaultSet::write( ostream &os, int style, bool header )
{
    bool havez = (style & 1) ? true : false;
    bool linestring = (style & 2) ? true : false;

    string prefix, zstring, suffix;
    zstring = havez ? " Z" : "";
    if( linestring ) { prefix="MULTILINESTRING"+zstring+"(("; suffix="))"; }
    else { prefix="POLYGON"+zstring+"(("; suffix="))"; }

    if( header )
        os << "id\tstrike\tdip\tUss\tUds\tUts\tslip\tmin_depth\tmax_depth\tname\tshape" << endl;

    int nflt = 0;

    // Note: this code is written for the case that every fault has just
    // one segment - GNS fault model loaded as multiple faults rather than
    // single segmented fault.

    list<string>::const_iterator in = names.begin();
    list<SegmentedFault *>::const_iterator i = faults.begin();
    for( in = names.begin(), i = faults.begin(); i != faults.end(); in++, i++ )
    {
        SegmentedFault *f = *i;
        string name = *in;
        replace(name.begin(),name.end(),'\t',' ');
        nflt++;

        os << setiosflags( ios::fixed )
           << setprecision(1)
           << nflt << "\t" << f->strike() << "\t" << f->dip() << "\t";

        double *sv = f->SlipVector(0,0);
        double slip=sqrt(sv[0]*sv[0]+sv[1]*sv[1]+sv[2]*sv[2]);
        os << setprecision(4);
        os << sv[0] << "\t" << sv[1] << "\t" << sv[2] << "\t" << slip << "\t";
        os << setprecision(7);
        double x, y, z, minz, maxz;
        f->FaultLocation(0, 0, x, y, minz );
        f->FaultLocation(0, f->NDipSegments(), x, y, maxz );
        if( minz > maxz ) { z = minz; minz = maxz; maxz = z; }
        os << minz << "\t" << maxz << "\t";
        os << name << "\t";
        os << prefix;

        bool firstpt = true;
        for( int ipt = 0; ipt < 5; ipt++ )
        {
            int is = 0, id = 0;
            if( ipt == 2 || ipt == 3 ) id = f->NDipSegments();
            if( ipt == 1 || ipt == 2 ) is = f->NStrikeSegments();
            f->FaultLocation(is, id, x, y, z );
            z = -z;
            double lat, lon;
            proj->LatLon(x,  y, lon, lat );
            if( firstpt ) firstpt = false;
            else os << ", ";
            os << lon << " " << lat;
            if( havez ) os << " " << z;

        }
        // For line string representation show the surface projection
        // of the fault.
        if( linestring )
        {
            double x, y, z, lon, lat;
            os << "),(";
            f->FaultLocation(0, -1, x, y, z );
            z = -z;
            proj->LatLon(x,  y, lon, lat );
            os << lon << " " << lat;
            if( havez ) os << " " << z;
            f->FaultLocation(f->NStrikeSegments(), -1, x, y, z );
            z = -z;
            proj->LatLon(x,  y, lon, lat );
            os << ", " << lon << " " << lat;
            if( havez ) os << " " << z;
        }
        os << suffix << endl;
    }
    return os;
}

static void calcStrainComponents( double *strain, double &dil, double &rot, double &shear, double &err )
{
    double uxx = strain[0]*1.0e6;
    double uxy = strain[1]*1.0e6;
    double uyx = strain[2]*1.0e6;
    double uyy = strain[3]*1.0e6;

    dil=(uxx+uyy)/2.0;
    rot=(uyx-uxy)/2.0;
    shear=sqrt((uxx-uyy)*(uxx-uyy) + (uxy+uyx)*(uxy+uyx))/2;

    double ea = (uxx*uxx+uyx*uyx)/2;
    double eb = (uxy*uxy+uyy*uyy)/2;
    double ec = uxx*uxy+uyx*uyy;

    err = sqrt(fabs(ea+eb) + sqrt((ea-eb)*(ea-eb)+ec*ec));
}

void split_string( const string &source, const string &delim, list<string> &parts )
{
    string::size_type pos = 0;
    while( 1 )
    {
        string::size_type end= source.find(delim,pos);
        if( end == string::npos ) break;
        if( end > pos ) parts.push_back( source.substr(pos,end-pos) );
        pos = end+1;
    }
    parts.push_back( source.substr(pos) );
}

void help()
{
    string helpfile(get_image_path());
    helpfile += ".help";
    ifstream hf(helpfile.c_str());
    if( hf)
    {
        string line;
        while(hf)
        {
            getline(hf,line);
            cout << line << endl;
        }
        hf.close();
        exit(0);
    }
}

int main( int argc, char *argv[] )
{
    char *wktfile = 0;
    int wkttype = 3;
    bool isvalid = true;
    bool compare = false;
    bool havenames = false;
    bool showlength = false;
    bool calcdiff = false;
    bool calcstrain = false;
    bool calctilt = false;
    bool applysfconv = false;
    char delim='\t';
    char idelim;
    int nskip = 0;
    int llprecision = 6;
    int dxyprecision = 4;
    int strnprecision = 4;
    Projection *proj=0;

    while( argc > 1 && argv[1][0] == '-' && argv[1][1] )
    {
        double lon, lat;
        switch( argv[1][1] )
        {
        case 'w':
        case 'W':
            if( argc > 2 )
            {
                wktfile = argv[2];
                if( argv[1][2] == 'p' || argv[1][2]=='P') wkttype=1;
                argv++;
                argc--;
            }
            else
            {
                cerr << "-w requires wkt_file argument" << endl;
                isvalid = false;
            }
            break;
        case 'p':
        case 'P':
            if( argc > 2 && _stricmp(argv[2],"none") == 0 )
            {
                if( proj ) delete proj;
                proj=new NullProjection();
                argv++;
                argc--;
            }
            else if( argc > 3
                     && sscanf(argv[2],"%lf",&lon) == 1
                     && sscanf(argv[3],"%lf",&lat) == 1  )
            {
                if( proj ) delete proj;
                proj=new GNSProjection(lon,lat);
                argv+=2;
                argc-=2;
            }
            else
            {
                cerr << "-p requires lon and lat arguments or \"none\"" << endl;
                isvalid = false;
                break;
            }
            break;
        case 'f':
        case 'F':
            applysfconv=true;
            break;
        case 'n':
        case 'N':
            havenames = true;
            break;
        case 'l':
        case 'L':
            showlength = true;
            break;
        case 'c':
        case 'C':
            compare = true;
            break;
        case 'd':
        case 'D':
            calcdiff = true;
            break;
        case 's':
        case 'S':
            calcstrain = true;
            break;
        case 't':
        case 'T':
            calctilt = true;
            break;
        case 'v':
        case 'V':
            delim=',';
            break;
        case 'h':
        case 'H':
            if( argc < 2 || sscanf(argv[2],"%d",&nskip) != 1 )
            {
                cerr << "-h requires number of lines to skip" << endl;
                isvalid = false;
            }
            else
            {
                argv++;
                argc--;
            }
            break;
        case 'x':
        case'X':
            llprecision += 4;
            dxyprecision += 4;
            strnprecision += 4;
            break;
        default:
            cerr << "Invalid argument " << argv[1] << endl;
            isvalid = false;
        }
        argv++;
        argc--;
    }

    if( ! isvalid )
    {
        cerr << "Failed to run - invalid arguments" << endl;
        return 0;
    }

    if( argc < 2 )
    {
        help();
    }

    if( ! (argc == 4 || (argc ==2 && wktfile)))
    {
        cerr << "Require parameters: [-w[l|p] wktfile] fault_model_file test_point_file output_file\n";
        return 0;
    }

    list<FaultSet *> faultlist;
    list<string> filenames;
    split_string(argv[1],"+",filenames);
    
    FaultSet faults(proj);

    // Read the fault model

    for( list<string>::iterator i = filenames.begin(); i != filenames.end(); i++ )
    {
        string filename(*i);
        if( filename.length() == 0 ) continue;

        double factor=1.0;
        string::size_type end= filename.find('*');
        if( end != string::npos ) 
        {
            string factorstr=filename.substr(0,end);
            filename=filename.substr(end+1);
            if( sscanf(factorstr.c_str(),"%lf",&factor) != 1 )
            {
                cerr << "Invalid fault model file scale factor " << factor << " defined for " << filename << "\n";
                return 0;
            }
            // cout << "Scale factor = " << factor << "\n";
        }

        ifstream f(filename.c_str());

        if( ! f.good() )
        {
            filename += ".model";
            f.open(filename.c_str());
        }

        if( ! f.good() )
        {
            const char *buf = getenv("FAULT_MODEL_DIR");
            if( buf ) 
            {
                string filename(buf);
                filename += "/";
                filename += argv[1];
                f.open(filename.c_str());
                if( ! f.good())
                {
                    filename += ".model";
                    f.open(filename.c_str());
                }
            }
        }

        if( ! f.good() )
        {
            cerr << "Cannot read fault model file " << argv[1] << endl;
            return 0;
        }

        FaultSet *faults = new FaultSet( proj );
        faults->setApplyConvergence(applysfconv);
        if( ! faults->ReadGNSDefinition(f, nskip) ) { return 0; }
        faults->setFactor(factor);
        faultlist.push_back(faults);
    }

    if( wktfile )
    {
        ofstream wkt(wktfile);
        if( ! wkt.good())
        {
            cerr << "Cannot open wkt output file " << argv[4] << endl;
            return 0;
        }
        bool header = true;
        for( list<FaultSet *>::iterator f = faultlist.begin(); f != faultlist.end(); f++ )
        {
            (*f)->write(wkt,wkttype,header);
            header=false;
        }
        wkt.close();
    }
    if( argc == 2 ) exit(0);

    // Open the input and output files..
    
    istream *in;
    string fname(argv[2]);
    idelim=delim;
    if( fname.substr(0,5) == "grid:" || fname.substr(0,6) == "point:" )
    {
      replace(fname.begin(),fname.end(),':',' ');
      if( fname.substr(0,5) == "point" ) fname=fname.substr(6);
      in = new istringstream(fname);
      idelim=' ';
      havenames=false;
      compare=false;
    }
    else if( fname == "-" )
    {
      in = &cin;
    }
    else
    {
      in = new ifstream(argv[2]);
      if( ! in->good() )
      {
	  cerr << "Cannot open input test point file " << argv[2] << endl;
	  return 0;
      }
    }

    double demax = 0.0;
    double dnmax = 0.0;
    double dumax = 0.0;

    ofstream *outs = 0;
    string fnameout(argv[3]);
    if( fnameout != "-" )
    {
        outs = new ofstream(argv[3]);
        if( ! outs->good() )
        {
            cerr << "Cannot open output file " << argv[3] << endl;
            return 0;
        }
    }
    ostream &out =  outs ? *outs : cout;

    if( havenames ) out << "name" << delim;
    out << "lon" << delim << "lat" << delim << "de" << delim << "dn" << delim << "du";
    if( showlength ) out << delim << "ds";
    if( calcdiff ) out << delim << "ddede" << delim << "ddedn" << delim << "ddnde" << delim << "ddndn" << delim << "ddude" << delim << "ddudn" ;
    if( calcstrain ) out << delim << "dil" << delim << "rot" << delim << "shear" << delim << "err";
    if( calctilt ) out << delim << "tilte" << delim << "tiltn" << delim << "tiltmax";
    if( compare )
    {
        out << "" << delim << "obs_de" << delim << "obs_dn" << delim << "obs_du";
        if( showlength) out << "" << delim << "obs_ds";
        out << "" << delim << "dif_de" << delim << "dif_dn" << delim << "dif_du";
        if( showlength) out << delim << "dif_ds";
    }
    out << endl;
    out << setiosflags( ios::fixed );
    string buffer;
    int lineno=0;
    double ppm=1.0e6;
    while( getline(*in, buffer) )
    {
        lineno++;
        // Skip comments
        int p = buffer.find('#');
        if( p != string::npos) { buffer.erase(p); }
        if( buffer.find_first_not_of(' ') == string::npos ) continue;

        double lon,lat, uxyz[3], duxy[4], dz[2];
        double *strain = (calcstrain || calcdiff) ? duxy : 0;
        double *tilt = (calctilt || calcdiff) ? dz : 0;
        double &ux = uxyz[0];
        double &uy = uxyz[1];
        double &uz = uxyz[2];

        double lon0, lat0, lon1, lat1, dlon, dlat;
        int nln, nlt;
        replace(buffer.begin(),buffer.end(),idelim,' ');
        stringstream s(buffer);
        string gridstr;

        if( ! compare && ! havenames)
        {
            s >> gridstr >> lon0 >> lat0 >> lon1 >> lat1 >> nln >> nlt;
            if( s ) std::transform(gridstr.begin(),gridstr.end(),gridstr.begin(),::tolower);
        }
        if( ! s || gridstr != "grid" )
        {
            s.clear();
            s.seekg(0);
            string name;
            if( havenames ) s >> name;
            s >> lon0 >> lat0;
            double oux, ouy, ouz;
            if( compare )
            {
                s >> oux >> ouy >> ouz;
            }
            if( s.fail() )
            {
                if( lineno > 1 ) cerr << "Invalid test point at line " 
                                      << lineno << ": " << buffer << endl;
                continue;
            }

            uxyz[0] = uxyz[1] = uxyz[2] = 0.0;
            duxy[0] = duxy[1] = duxy[2] = duxy[3] = 0.0;
            dz[0] = dz[1] = 0.0;

            for( list<FaultSet *>::iterator f = faultlist.begin(); f != faultlist.end(); f++ )
            {
                (*f)->AddOkada( lon0, lat0, uxyz, strain, tilt );
            }
            if( havenames ) out << name << delim;
            out << setprecision(llprecision)
                << lon0 <<  "" << delim << lat0
                << setprecision(dxyprecision)
                << delim << ux << delim << uy << delim << uz;
            if( showlength )
            {
                double us = sqrt(ux*ux+uy*uy);
                out << delim << us;
            }
            if( calcdiff )
            {
                out << setprecision(strnprecision)
                    << delim << strain[0]*1000000.0
                    << delim << strain[1]*1000000.0
                    << delim << strain[2]*1000000.0
                    << delim << strain[3]*1000000.0
                    << delim << tilt[0]*1000000.0
                    << delim << tilt[1]*1000000.0;
            }
            if( calcstrain )
            {
                double dil, rot, shear, err;
                calcStrainComponents(strain,dil,rot,shear,err);
                out << setprecision(strnprecision)
                    << delim << dil << delim << rot
                    << delim << shear << delim << err;
            }
            if( calctilt )
            {
                double tiltmax=sqrt(tilt[0]*tilt[0]+tilt[1]*tilt[1])*ppm;
                out << setprecision(strnprecision)
                    << delim << tilt[0]*ppm << delim << tilt[1]*ppm << delim << tiltmax;
            }
            if( compare )
            {
                out << setprecision(dxyprecision)
                    << delim << oux << delim << ouy << delim << ouz;
                if( showlength )
                {
                    double us = sqrt(oux*oux+ouy*ouy);
                    out << delim << us;
                }
                oux -= ux; ouy -= uy; ouz -= uz;
                if( fabs(oux) > demax ) demax = fabs(oux); 
                if( fabs(ouy) > dnmax ) dnmax = fabs(ouy); 
                if( fabs(ouz) > dumax ) dumax = fabs(ouz); 
                out << delim << oux << delim << ouy << delim << ouz;
                if( showlength )
                {
                    double us = sqrt(oux*oux+ouy*ouy);
                    out << delim << us;
                }
            }
            out << endl;
        }
        else
        {
            dlon = (lon1 - lon0)/nln;
            dlat = (lat1 - lat0)/nlt;
            nln++; nlt++;
            for( lat = lat0; nlt > 0; nlt--, lat += dlat )
            {
                int nlnrow = nln;
                for( lon = lon0; nlnrow > 0; nlnrow--, lon += dlon )
                {
                    uxyz[0] = uxyz[1] = uxyz[2] = 0.0;
                    duxy[0] = duxy[1] = duxy[2] = duxy[3] = 0.0;
                    dz[0] = dz[1] = 0.0;

                    for( list<FaultSet *>::iterator f = faultlist.begin(); f != faultlist.end(); f++ )
                    {
                        (*f)->AddOkada( lon, lat, uxyz, strain, tilt );
                    }
                    out << setprecision(llprecision)
                        << lon <<  delim << lat << delim
                        << setprecision(dxyprecision)
                        << ux << delim << uy << delim << uz;
                    if( showlength )
                    {
                        double us = sqrt(ux*ux+uy*uy);
                        out << delim << us;
                    }
                    if( calcdiff )
                    {
                        out << setprecision(strnprecision)
                            << delim << strain[0]*1000000.0
                            << delim << strain[1]*1000000.0
                            << delim << strain[2]*1000000.0
                            << delim << strain[3]*1000000.0
                            << delim << tilt[0]*1000000.0
                            << delim << tilt[1]*1000000.0;
                    }
                    if( calcstrain )
                    {
                        double dil, rot, shear, err;
                        calcStrainComponents(strain,dil,rot,shear,err);
                        out << setprecision(strnprecision)
                            << delim << dil << delim << rot
                            << delim << shear << delim << err;
                    }
                    if( calctilt )
                    {
                        double tiltmax=sqrt(tilt[0]*tilt[0]+tilt[1]*tilt[1])*ppm;
                        out << setprecision(strnprecision)
                            << delim << tilt[0]*ppm << delim << tilt[1]*ppm << delim << tiltmax;
                    }
                    out << endl;
                }
            }
        }
    }
    if( in != &cin ) delete(in);
    if( outs ) outs->close();
    if( compare )
    {
        cout << "Max displacement differences (" << demax << ", " << dnmax << ", " << dumax << ")" << endl;
    }
    return 0;
}
