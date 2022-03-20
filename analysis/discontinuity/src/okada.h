#ifndef OKADA_H
#define OKADA_H

// Code for calculating the surface displacement due to a dislocation on
// a rectangular fault.  Based on:
//
//	Okada Y., Surface deformation due to shear and tensile faults in a
//	   half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.

// Reference system used for defining a fault plane (or segment of one).
// Defines a reference point, bearing, and dip

class FaultRefSys
{
    public:
        // Define the fault based reference system in terms of a point
        // on the fault plance x0, y0, d0 (depth) in metres, and a 
        // strike direction and dip in degrees.  
        // Strike is measured from the X axis towards the Y axis.
        // Dip is measured down to the left of the fault, viewing along strike 
        // (assuming a left hand coordinate system)
        
        FaultRefSys();
        FaultRefSys( double x0, double y0, double d0, double strike, double dip );
        FaultRefSys( const FaultRefSys &src );
        FaultRefSys& operator= (const FaultRefSys &src );
        void XYD( double fs, double fd, double &x, double &y, double &d );
        

    public:
        void Setup();
        // Coordinates of reference point on top of fault segment 
        double x0, y0, d0; 
        // Strike of fault (degrees), cos and sin */
        double strike, coss, sins; 
        // Dip (degrees) measured as per Okada, 
        // up to left of fault or down to right, 
        // viewing along strike */
        double dip, cosd, sind;  
        // True if the fault will be treated as vertical
        bool vertical;
        // Poisson ratio
        double poisson;
        // Medium constant ( mu/(lambda+mu) )
        double medium;
};

// Defines a segmented fault definition
// 

class SegmentedFault
{
    public:

        SegmentedFault();
        // Define a gridded rectangular plane fault
        // Initiallize with coordinate system
        SegmentedFault( 
                double x0, double y0, double d0, 
                double strike, double dip );
        SegmentedFault( const FaultRefSys &fsys );
        // Set an id to identify the particular fault
        void SetId( int id ){ this->id=id; }
        // Set the fault reference system
        void SetFaultRefSys( const FaultRefSys &fsys ){ this->fsys = fsys; }
        // Set number of segments along the fault and offsets of segment
        // boundaries (offsets indexed from 0 to nsegs = nsegs+1 values)
        void SetStrikeSegBreaks( int nsegs, double *offsets );
        // Set number of segments down dip and offsets of segment boundaries
        void SetDipSegBreaks( int nsegd, double *offsetd );
        // Set the slip vector on the isegs,isegd segment as Uss (strike slip),
        // Uds (dip slip), and Uts (tensile) components.  The isegs segment
        // runs from offsets[isegs] to offset[isegs+1] (ie numbered from
        // 0 to isegs-1).  Similarly for isegd.
        void SetFaultSlip( int isegs, int isegd, double Uss, double Uds, double Uts );

        // One hit definition of a simple rectangular fault
        // x0,y0,d0,strike,dip as for FaultRefSys
        // lens and lend are distance along strike and dip from reference point
        // Uss, Uds, Uts are strike slip, dip slip, and tensile movement
        SegmentedFault( 
                double x0, double y0, double d0, 
                double strike, double dip,
                double lens, double lend,
                double Uss, double Uds, double Uts
                );

        // Destructor
        ~SegmentedFault(){ FreeAll(); }

        // Add the dislocations and strains from the formulae to cumulative
        // values.  strain can be 0 to skip strain calcs
        //
        bool AddOkada( double x, double y, double *dislocation, double *strain=0, double *tilt=0, double factor=1.0 );
        // Calculate dislocation using Okada formulae at a specified point
        // Return false if this cannot be done (eg fault not completely defined)

        bool OkadaDislocation( double x, double y, double &ux, double &uy, double &uz );

        // Calculate strain using Okada formulae at a specified point
        // Return false if this cannot be done (eg fault not completely defined)

        bool OkadaStrain( double x, double y, double &uxx, double &uxy, double &uyx, double &uyy );

        int NStrikeSegments(){ return nsegs; }
        int NDipSegments() { return nsegd; }
        int Id() { return id; }

        // Get the coordinates of a point defined  by distance along and across fault
        // Note: use isegd=-1 for surface projection of fault

        void FaultLocation( int isegs, int isegd, double &x, double &y, double &d );
        // Calculate the slip vector for a fault segment
        double *SlipVector(int isegs,int isegd){ return slipvector + (isegd*nsegs + isegs)*3; }
        double strike(){ return fsys.strike; }
        double dip(){ return fsys.dip; }

    private:
        void Setup();
        void FreeSlipComponents();
        void FreeAll();
        FaultRefSys fsys;
        int id;
        int nsegs;
        int nsegd;
        double *offsets;
        double *offsetd;
        double *slipvector;
};


#endif
