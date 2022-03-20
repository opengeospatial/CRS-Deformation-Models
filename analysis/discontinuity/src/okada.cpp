
// Make M_PI available to MS compiler
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>
#include <iostream>

#include "okada.h"

// Added for using cout << ...
#include <iomanip>
using namespace std;

const double SMALLCOS = 0.000001;
const double SMALLSLIP = 0.000001;
const double SMALLOFFSET = 0.000001;

// Implementation of FaultRefSys - defines the fault based coordinate system,

FaultRefSys::FaultRefSys() :
    x0(0), y0(0), d0(0), strike(0), dip(90.0)
{
    Setup();
}

FaultRefSys::FaultRefSys( double x0, double y0, double d0, double strike, double dip )
    : x0(x0), y0(y0), d0(d0), strike(strike), dip(dip)
{
    Setup();
}

FaultRefSys::FaultRefSys( const FaultRefSys &src )
    : x0(src.x0), y0(src.y0), d0(src.d0), strike(src.strike), dip(src.dip)
{
    Setup();
}

FaultRefSys & FaultRefSys::operator=( const FaultRefSys &src )
{
    x0 = src.x0;
    y0 = src.y0;
    d0 = src.d0;
    strike = src.strike;
    dip = src.dip;
    Setup();
    return *this;
}

void FaultRefSys::Setup()
{
    // Calculate required trigonometric values
    double conv = M_PI/180;
    double rstrike = strike * conv;
    coss = cos(rstrike);
    sins = sin(rstrike);

    // Note: Convert dip from measuring down to the left, 
    // to measuring up to the left, to be consistent with Okada model
    double rdip = -dip * conv;
    cosd = cos(rdip);
    sind = sin(rdip);
    vertical = fabs(cosd) < SMALLCOS;

    poisson = 0.25;
    medium  = 1.0 - 2.0*poisson;
}

void FaultRefSys::XYD( double fs, double fd, double &x, double &y, double &d )
{
    double yr = fd*cosd;
    x = (fs*coss - yr*sins)+x0;
    y = (fs*sins + yr*coss)+y0;
    d = d0-fd*sind;
}


class OkadaPoint
{
    public:
        // Initiallize the point with the fault plane reference system
        OkadaPoint( FaultRefSys &fsys ) : fsys(fsys) {}
        // Set the surface point location
        void SetSurfacePoint( double xs, double ys ){ this->xs = xs; this->ys = ys; }
        // Set the fault point location (strike offset, dip offset) from
        // reference point x0, y0, d0
        void SetFaultPoint( double fs, double fd ){ this->fs = fs; this->fd = fd; }

        void AddOkada( double U1, double U2, double U3, double *dislocation, double *strain, double *tilt, double factor );
    private:
        // Input parameters
        FaultRefSys &fsys;
        double xs;
        double ys;
        double fs;
        double fd;

        // Fault coordinate system parameters
        double xi;
        double eta;
        double q;
        double y;
        double d;
        double R;
        double X;

        void SetPosition();
};


// Calculate the fault based coordinate system terms xi, eta, p, q, y, d

void OkadaPoint::SetPosition()
{
    // Convert xs, ys to points xr, yr relative to the surface point x0, y0 in

    double tmp = xs - fsys.x0;
    double yr = ys - fsys.y0;
    double xr = tmp*fsys.coss+yr*fsys.sins;
    yr = yr*fsys.coss - tmp*fsys.sins;

    // Convert to fault plane aligned coordinate system (xi, eta, q), relative 
    // to the evaluation point on the fault plane.
    
    xi = xr - fs;
    y = yr - fd*fsys.cosd;
    d = fsys.d0 - fd*fsys.sind;
    
    eta = y*fsys.cosd + d*fsys.sind;
    q = y*fsys.sind - d*fsys.cosd;

    static int id = 0;
    id++;

    // Calculate the X and R terms
    
    X = xi*xi+q*q;
    R = X + eta*eta;
    X = sqrt(X);
    R = sqrt(R);
}

void OkadaPoint::AddOkada( double U1, double U2, double U3, 
        double *dislocation, double *strain, double *tilt, double factor )
{

    SetPosition();
    if( R < SMALLOFFSET ) return;

    bool calcStrains = (strain != 0);
    bool calcTilts = (tilt != 0);
    bool strikeslip = fabs(U1) > SMALLSLIP;
    bool dipslip = fabs(U2) > SMALLSLIP;
    bool tensile = fabs(U3) > SMALLSLIP;
    if( !(strikeslip || dipslip || tensile )) return;

    double ux = 0, uy = 0, uz = 0;
    double uxx = 0, uxy=0, uyx=0, uyy=0;
    double uzx = 0, uzy=0;

    double I1=0.0, I2=0.0, I3=0.0, I4=0.0, I5=0.0;
    double J1=0.0, J2=0.0, J3=0.0, J4=0.0;
    double K1=0.0, K2=0.0, K3=0.0;
    double Axi = 0.0, Aeta=0.0;
    double cosd = fsys.cosd;
    double sind = fsys.sind;
    double scd = sind*cosd;
    double medium = fsys.medium;

    double Rd = R+d;
    double RRd = R*Rd;
    double qRd2 = q/(Rd*Rd);
    double Reta = R+eta;
    double RReta = R*Reta;
    double RRxi = R*(R+xi);
    // Handle potential numerical instability in calculated R+xi
    if( xi < 0.0 && RRxi < R*0.000001 )
    {
        RRxi=R*(eta*eta+q*q)/(-2.0*xi);
    }
    double R3 = R*R*R;
    double lnReta = log(Reta);

    // If the fault is vertical, then some terms become ...
    if( fsys.vertical )
    {
        I5 = - medium * (xi*sind/Rd);
        I4 = - medium * (q/Rd);
        I3 =   (medium/2.0)*(eta/Rd + y*qRd2 - log(R + eta));
        I1 = - (medium/2.0)*xi*qRd2;

        if( calcStrains || calcTilts )
        {
            K1 = medium*eta*qRd2/R;
            K3 = medium*(sind/Rd)*(eta*eta/RRd-1);
            K2 = medium*(q*cosd/RReta-sind/R)-K3;
            J1 = medium*qRd2*(xi*xi/(RRd)-0.5);
            J2 = medium*xi*(sind/(Rd*Rd))*(q*q/(RRd) - 0.5);
        }
    }
    else
    {

        I5 = 2.0 * (medium/cosd) * atan((eta*(X+q*cosd)+X*(R+X)*sind)/(xi*(R+X)*cosd) );
        I4 = (medium/cosd) * (log(Rd) - sind*lnReta);
        I3 = medium * (y/(Rd*cosd) - lnReta) + sind*I4/cosd;
        I1 = -medium*(xi/(cosd*Rd)) - sind*I5/cosd;

        if( calcStrains || calcTilts )
        {
            K1 = medium*xi*(1.0/RRd-sind/RReta)/cosd;
            K3 = medium*(q/RReta-y/RRd)/cosd;
            K2 = medium*(q*cosd/RReta-sind/R)-K3;
            J1 = (medium*(xi*xi/RRd-1.0)/Rd - sind*K3)/cosd;
            J2 = (medium*xi*y/(RRd*Rd) - sind*K1)/cosd;
        }
    }

    // Only need I2 for strike-slip component
    if( strikeslip ) I2 = -medium*lnReta - I3;

    if( calcStrains || calcTilts )
    {

        J3 = - medium*xi/RReta - J2;
        J4 = -medium*(cosd/R + q*sind/RReta) - J1;

        Axi = (2*R+xi)/(R*RRxi*RRxi);
        Aeta = (2*R+eta)/(R*RReta*RReta);
    }

    double atpeqR = 0.0;
    if( q != 0.0 ) atpeqR = atan( (xi*eta)/(q*R) );

    if( strikeslip )
    {
        double U2pi = U1/(2.0*M_PI);
        ux -= U2pi*(xi*q/RReta + atpeqR + I1*sind);
        uy -= U2pi*(y*q/RReta+q*cosd/Reta+I2*sind);
        uz -= U2pi*(d*q/RReta+q*sind/Reta+I4*sind);

        if( calcStrains )
        {
            double xi3 = xi*xi*xi;
            uxx += U2pi*(xi*xi*q*Aeta - J1*sind);
            uxy += U2pi*(xi3*d/(R3*(eta*eta+q*q))-sind*(xi3*Aeta+J2));
            uyx += U2pi*(xi*q*cosd/R3+sind*(xi*q*q*Aeta-J2));
            uyy += U2pi*(y*q*cosd/R3 + (q*q*q*Aeta*sind-2*q*sind/RReta-(xi*xi+eta*eta)*cosd/R3-J4)*sind);
        }

        if( calcTilts )
        {
            uzx += U2pi*(sind*(xi*q/R3-K1)-cosd*xi*q*q*Aeta);
            uzy += U2pi*(cosd*d*q/R3+sind*(xi*xi*q*Aeta*cosd-sind/R+y*q/R3-K2));
        }
    }
    if( dipslip )
    {
        double U2pi = U2/(2.0*M_PI);
        ux -= U2pi*(q/R - I3*scd);
        uy -= U2pi*(y*q/RRxi + cosd*atpeqR - I1*scd);
        uz -= U2pi*(d*q/RRxi + sind*atpeqR - I5*scd);

        if( calcStrains )
        {
            uxx += U2pi*(xi*q/R3+J3*scd);
            uxy += U2pi*(y*q/R3-sind/R+J1*scd);
            uyx += U2pi*(y*q/R3+q*cosd/RReta+J1*scd);
            uyy += U2pi*(y*y*q*Axi-(2*y/RRxi+xi*cosd/RReta)*sind+J2*scd);
        }

        if( calcTilts )
        {
            uzx += U2pi*(d*q/R3+sind*q/RReta+K3*scd);
            uzy += U2pi*(y*d*q*Axi-sind*(2*d/RRxi+xi*sind/RReta)+K1*scd);
        }
    }
    if( tensile )
    {
        double U2pi = U3/(2.0*M_PI);
        double term = xi*q/RReta-atpeqR;
        ux += U2pi*(q*q/RReta - I3*sind*sind);
        uy -= U2pi*(d*q/RRxi+sind*term+I1*sind*sind);
        uz += U2pi*(y*q/RRxi+cosd*term-I5*sind*sind);

        if( calcStrains )
        {
            double sind2 = sind*sind;
            uxx -= U2pi*(xi*q*q*Aeta+J3*sind2);
            uxy += U2pi*(d*q/R3+xi*xi*q*Aeta*sind-J1*sind2);
            uyx -= U2pi*(q*q*cosd/R3+q*q*q*Aeta*sind+J1*sind2);
            uyy -= U2pi*((y*cosd-d*sind)*q*q*Axi - 2*q*scd/RRxi-(xi*q*q*Aeta-J2)*sind2);
        }

        if( calcTilts )
        {
            uzx += -U2pi*(q*q*sind/R3-q*q*q*Aeta*cosd+K3*sind*sind);
            uzy += -U2pi*((y*sind+d*cosd)*q*q*Axi+xi*q*q*Aeta*scd-(2*q/RRxi-K1)*sind*sind);
        }
    }

    // Convert ux, uy from strike oriented to world oriented coordsys
    
    double coss = fsys.coss;
    double sins = fsys.sins;
    dislocation[0] += (ux*coss-uy*sins)*factor;
    dislocation[1] += (uy*coss+ux*sins)*factor;
    dislocation[2] += uz*factor;

    if( calcStrains )
    {

        double cs=coss*sins;
        double c2=coss*coss;
        double s2=sins*sins;
        strain[0] += (c2*uxx+cs*(uxy+uyx)+s2*uyy)*factor;
        strain[1] += (cs*(uyy-uxx)+c2*uxy-s2*uyx)*factor;
        strain[2] += (cs*(uyy-uxx)+c2*uyx-s2*uxy)*factor;
        strain[3] += (s2*uxx-cs*(uxy+uyx)+c2*uyy)*factor;
    }

    if( calcTilts )
    {
        tilt[0] += (uzx*coss-uzy*sins)*factor;
        tilt[1] += (uzy*coss+uzx*sins)*factor;
    }
}

SegmentedFault::SegmentedFault()
{
    Setup();
}

SegmentedFault::SegmentedFault( const FaultRefSys &fsys ) :
    fsys(fsys)
{
    Setup();
}

SegmentedFault::SegmentedFault( 
                double x0, double y0, double d0, 
                double strike, double dip ) :
    fsys( x0, y0, d0, strike, dip )
{
    Setup();
}

SegmentedFault::SegmentedFault( 
                double x0, double y0, double d0, 
                double strike, double dip,
                double lens, double lend,
                double Uss, double Uds, double Uts ) :
    fsys( x0, y0, d0, strike, dip )
{
    Setup();
    double seg[2];
    seg[0] = 0;
    seg[1] = lens;
    SetStrikeSegBreaks( 1, seg );
    seg[1] = lend;
    SetDipSegBreaks( 1, seg );
    SetFaultSlip( 0, 0, Uss, Uds, Uts );
}

void SegmentedFault::Setup()
{
    nsegs = nsegd = 0;
    offsets = 0;
    offsetd = 0;
    slipvector = 0;
}

void SegmentedFault::FreeSlipComponents()
{
    if( slipvector ) delete [] slipvector;
    slipvector = 0;
}

void SegmentedFault::FreeAll()
{
    FreeSlipComponents();
    if( offsets ) delete [] offsets;
    offsets = 0;
    if( offsetd ) delete [] offsetd;
    offsetd = 0;
    nsegs = nsegd = 0;
}

void SegmentedFault::SetStrikeSegBreaks( int pnsegs, double *poffsets )
{
    FreeSlipComponents();
    if( offsets ) delete [] offsets;
    nsegs = pnsegs;
    offsets = new double[nsegs+1];
    memcpy(offsets, poffsets, sizeof(double)*(nsegs+1));
}

void SegmentedFault::SetDipSegBreaks( int pnsegd, double *poffsetd )
{
    FreeSlipComponents();
    if( offsetd ) delete [] offsetd;
    nsegd = pnsegd;
    offsetd = new double[nsegd+1];
    memcpy(offsetd, poffsetd, sizeof(double)*(nsegd+1));
}

void SegmentedFault::SetFaultSlip( int isegs, int isegd, double Uss, double Uds, double Uts )
{
    if( ! slipvector ) slipvector = new double [nsegs*nsegd*3];
    double *vector = SlipVector(isegs,isegd);
    vector[0] = Uss;
    vector[1] = Uds;
    vector[2] = Uts;
}


bool SegmentedFault::AddOkada( double x, double y, double *dislocation, double *strain, double *tilt, double factor )
{
    if( ! slipvector ) return false;

    OkadaPoint pt(fsys);
    pt.SetSurfacePoint( x, y );
    for( int i = 0; i <= nsegs; i++ ) for( int j = 0; j <= nsegd; j++ )
    {
        pt.SetFaultPoint( offsets[i], offsetd[j] );
        double U1 = 0.0, U2 = 0.0, U3 = 0.0;
        for( int i1 = -1; i1 < 1; i1++ ) 
        {
            int isegs = i + i1;
            if( isegs < 0 || isegs == nsegs ) continue;
            for( int j1 = -1; j1 < 1; j1++ )
            {
                int isegd = j + j1;
                if( isegd < 0 || isegd == nsegd ) continue;
                double *vector = SlipVector(isegs, isegd );
                int ufactor = i1 == j1 ? 1 : -1;
                U1 += ufactor * vector[0];
                U2 += ufactor * vector[1];
                U3 += ufactor * vector[2];
                double pux, puy, puz;
                pt.AddOkada( U1, U2, U3, dislocation, strain, tilt, factor );
            }
        }
    }
    return true;
}

bool SegmentedFault::OkadaDislocation( double x, double y, double &ux, double &uy, double &uz )
{
    double uxyz[3] = {0.0, 0.0, 0.0};
    bool result = AddOkada(x, y, uxyz );
    ux = uxyz[0];
    uy = uxyz[1];
    uz = uxyz[2];
    return result;
}

bool SegmentedFault::OkadaStrain( double x, double y, double &uxx, double &uxy, double &uyx, double &uyy )
{
    double uxyz[3] = {0.0, 0.0, 0.0};
    double strain[4] = {0.0, 0.0, 0.0, 0.0};
    bool result = AddOkada(x, y, uxyz, strain );
    uxx = strain[0];
    uxy = strain[1];
    uyx = strain[2];
    uyy = strain[3];
    return result;
}

void SegmentedFault::FaultLocation( int isegs, int isegd, double &x, double &y, double &d )
{
    x = y = d = 0.0;
    if( isegs < 0 || isegs > nsegs || isegd < -1 || isegd > nsegd ) return;

    if( isegd >= 0 )
    {
        fsys.XYD(offsets[isegs],offsetd[isegd],x,y,d);
    }
    else if( fabs(fsys.sind) >= 0.02 )
    {
        double x0, y0, d0, x1, y1, d1;
        fsys.XYD(offsets[isegs],offsetd[nsegd],x0,y0,d0);
        fsys.XYD(offsets[isegs],offsetd[0],x1,y1,d1);
        if( fabs(fsys.sind) > 0.02 && fabs(d0-d1) != 0.0 )
        {
            x = x0 - d0*((x1-x0)/(d1-d0));
            y = y0 - d0*((y1-y0)/(d1-d0));
        }
        else
        {
            x = x0; y = y0; d = d0;
        }
    }
    return;
}
