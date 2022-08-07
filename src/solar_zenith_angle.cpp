#include <cmath>
#include "param.h"
#include "funcdef.h"

using namespace std;

/*************************************************************************
  Function solar_zenith
      calculate solar zenith angle in rad

  Input: parms - structure for system parameters
 *       t     - time in sec relative to the start of the simulation
 * Output: none. Result is solar zenith angles on all the grids in the global
 *               pointer *zenith to an one dimensional array

  By Jiannan Tu
  March 1, 2011, 1/13/2014, 3/26/2020
*************************************************************************/

void RECALC_08(int,int,int,int,int,double,double,double,double[],double[],
        double[],double[]);
void GEOMAG_08(double&, double&, double&, double&, double&, double&, int, 
        double[]);
void SPHCAR_08(double&, double&, double&, double&, double&, double&, int);

void solar_zenith(AppCtx *params, int ys, int ym, int zs, int zm)
{
    int    k, j;
    double eot, bt, lt, lst, shr;
    double gst, slong, srasn, sdec, hr, lat, pih=0.5*pi;
    int    imin, isec, iday;
    double x, y, z, xgeo, ygeo, zgeo;
    double AAP[106],G[106],H[106],REC[106];
    double rg, thetag, phig;
    
/* day number of the year*/
    iday=dayno(params->iyr, params->mon, params->idate);

    hr=params->sec/3600.0;
    imin=int((params->sec-int(hr)*3600.0)/60.0);
    isec=int(params->sec-int(hr)*3600.0-imin*60.0);

    RECALC_08(params->iyr,iday,(int)hr,imin,isec,params->vgsex,params->vgsey,params->vgsez,AAP,G,H,REC);

/* solar declination angle */
    SUN_08(params->iyr, iday, int(hr), imin, isec, gst, slong, srasn, sdec);

    //equation of time in minutes (bt in radians & constant rai = 360/365/deg)
    bt=rad*(iday-81.0);
    eot=9.87*sin(bt)-7.53*cos(bt)-1.5*sin(bt);

//#pragma omp parallel for private (i,phin,lt,lst,sha)
    for (k = zs; k < zs+zm; k++) {
        for (j = ys; j < ys+ym; j++) {
            SPHCAR_08(rr[0], theta[j], phi[k], x, y, z, 1);
            GEOMAG_08 (xgeo, ygeo, zgeo, x, y, z, -1, AAP);
            SPHCAR_08(rg, thetag, phig, xgeo, ygeo, zgeo, -1);

            //local time in hours = UT+longitude/15
            lt = hr + phig*deg/15.0;
            if (lt >= 24.0) lt=lt-24.0;

            /* local solar time in hours by adjusting local time with time correction
             * factor */
            lst=lt+eot/60.0;

            //solar hour angle in radians
            shr=0.2617993878*(lst-12.0);

            //latitude in radian
            lat=pih-thetag;

            //finally calculate solar zenith angle in radians
            zenith[k-zs][j-ys]=acos(sin(sdec)*cos(lat)+cos(sdec)*sin(lat)*cos(shr));

            //cos of zenith angle above which the sun is visible
            //for (i = xs; i < xs+xm; i++) 
            //    coschid[k-zs][j-ys][i-xs]=cos(pi-asin(1.0/(1.0+rr[i]*r0/Re)));
        }
    }
    
    return;
}
