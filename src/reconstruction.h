#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

#include "param.h"

inline void ambient_mfd(double rg, double thetag, double &Br0, double &Btheta0) 
{
    const double Mz=8.0e15;
    double r3;

    r3=rg*rg*rg;
    Br0=-2.0*Mz*cos(thetag)/r3;  //normalized
    Btheta0=-Mz*sin(thetag)/r3;  //normalized
}

inline double minmod(double aa, double bb){
    if (aa*bb < 0) return 0.0;
    else if(fabs(aa) < fabs(bb)) return aa;
    else return bb;
}

inline double limited_slope_r(Field ***xx, int i, int j, int k, int s)
{
    double aa, bb, rrb, rrt, y1, y2;
    int    Nrm1=Nr-1, Nrm2=Nr-2;

    if (i == 0) {
        aa = (xx[k][j][1].fx[s]-xx[k][j][i].fx[s])/(rC[1]-rC[i]);

        rrb=(rCm1-rC[1])/(rC[2]-rC[1]);  // rrb = (r_{-1} - r1)/(r2 - r1)
        if ((s >= 7 && s <= 9) || (s >= 19 && s <= 21)) {
            bb = (xx[k][j][i].fx[s]-(xx[k][j][1].fx[s]+rrb*(xx[k][j][2].fx[s]-xx[k][j][1].fx[s])))/(rC[i]-rCm1);
        }
        else { //y = y1 + (r -r1)/(r2 - r1) (y2 - y1) ---> y_{-1} = y1 + (r_{-1} - r1)/(r2 - r1) (y2 - y1)
            y1 = log(xx[k][j][1].fx[s]);
            y2 = log(xx[k][j][2].fx[s]);

            //bb = (xx[k][j][0].fx[s] - xx[k][j][-1].fx[s]) / (r0 - r_{-1})
            bb = (xx[k][j][i].fx[s]-exp(y1+rrb*(y2-y1)))/(rC[i]-rCm1);
        }
    }
    else if (i > 0 && i < Nr) {
        aa = (xx[k][j][i+1].fx[s]-xx[k][j][i].fx[s])/(rC[i+1]-rC[i]);
        bb = (xx[k][j][i].fx[s]-xx[k][j][i-1].fx[s])/(rC[i]-rC[i-1]);
    }
    else {
        rrt=(rC[a1]-rC[Nrm1])/(rC[Nrm1]-rC[Nrm2]);
        if ((s >= 7 && s <= 9) || (s >= 19 && s <= 21)) {
            aa = ((xx[k][j][Nrm1].fx[s]+rrt*(xx[k][j][Nrm1].fx[s]-xx[k][j][Nrm2].fx[s]))-xx[k][j][i].fx[s])
                 /(rC[a1]-rC[i]);
        }
        else {
            y1= log(xx[k][j][Nrm1].fx[s]);
            y2= log(xx[k][j][Nrm2].fx[s]);
            aa = (exp(y1+rrt*(y1-y2)) - xx[k][j][i].fx[s])/(rC[a1]-rC[i]);
        }

        bb = (xx[k][j][i].fx[s]-xx[k][j][Nrm1].fx[s])/(rC[i]-rC[Nrm1]);
    }

    return minmod(aa, bb);
}

inline double limited_slope_theta(Field ***xx, int i, int j, int k, int s)
{
    int kc;
    double cc, dd, sgn;

    if (s == 8 || s ==20) sgn = -1.0;
    else sgn = 1.0;

    if (j == 0) {
        kc = (k+a3/2) % a3;

        cc = (xx[k][j+1][i].fx[s]-xx[k][j][i].fx[s])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[s]-sgn*xx[kc][0][i].fx[s])/(2.0*thetaC[j]);
    }
    else if (j > 0 && j < Nthm) {
        cc = (xx[k][j+1][i].fx[s]-xx[k][j][i].fx[s])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[s]-xx[k][j-1][i].fx[s])/(thetaC[j]-thetaC[j-1]);
    }
    else {
        kc = (k+a3/2) % a3;

        cc = (sgn*xx[kc][Nthm][i].fx[s]-xx[k][j][i].fx[s])/(2.0*thetaC[j]);
        dd = (xx[k][j][i].fx[s]-xx[k][j-1][i].fx[s])/(thetaC[j]-thetaC[j-1]);
    }

    return minmod(cc, dd);
}

inline double limited_slope_phi(Field ***xx, int i, int j, int k, int s)
{
    int km, kp;
    double ee, ff;

    if (k == 0) km = Np; else km = k-1;
    if (k < Np) kp = k+1; else kp = 0;

    ee = (xx[kp][j][i].fx[s]-xx[k][j][i].fx[s])/dph;
    ff = (xx[k][j][i].fx[s]-xx[km][j][i].fx[s])/dph;

    return minmod(ee, ff);
}

/* reconstructed conservative or primitive fluid quantities within the cell (i, j, k) */
inline double reconstructed(Field ***xx, int i, int j, int k, int s, double rg, double thetag, double phig)
{
    double Uijk, slope_r = 0.0, slope_theta = 0.0, slope_phi = 0.0;

    if (rg - rC[i] != 0.0) slope_r = limited_slope_r(xx, i, j, k, s)*(rg - rC[i]);
    if (thetag - thetaC[j] != 0.0) slope_theta = limited_slope_theta(xx, i, j, k, s)*(thetag - thetaC[j]);
    if (phig - phi[k] != 0.0) slope_phi = limited_slope_phi(xx, i, j, k, s)*(phig - phi[k]);

    Uijk = xx[k][j][i].fx[s] + slope_r + slope_theta + slope_phi;

    return Uijk;
}

inline double flux_limited_slope_r(Field ***vv, int i, int j, int k, int s)
{
    double aa, bb, rrb;
    int    Nrm2=Nr-2;

    if (i == 0) {
        aa = (vv[k][j][1].fx[s]-vv[k][j][i].fx[s])/(rC[1]-rC[i]);

        rrb=(rCm1-rC[1])/(rC[2]-rC[1]);  // rrb = (r_{-1} - r1)/(r2 - r1)
        bb = (vv[k][j][i].fx[s]-(vv[k][j][1].fx[s]+rrb*(vv[k][j][2].fx[s]-vv[k][j][1].fx[s])))/(rC[i]-rCm1);
    }
    else if (i > 0 && i < Nr) {
        aa = (vv[k][j][i+1].fx[s]-vv[k][j][i].fx[s])/(rC[i+1]-rC[i]);
        bb = (vv[k][j][i].fx[s]-vv[k][j][i-1].fx[s])/(rC[i]-rC[i-1]);
    }
    else {
        aa = (vv[k][j][Nrm].fx[s]-vv[k][j][Nrm2].fx[s])/(rC[Nrm]-rC[Nrm2]);

        bb = (vv[k][j][Nrm2].fx[s]-vv[k][j][Nrm-3].fx[s])/(rC[Nrm2]-rC[Nr-3]);
    }

    return minmod(aa, bb);
}

/* reconstructed fluxes within the cell (i, j, k) */
inline double reconstructed_flux(Field ***vv, int i, int j, int k, int s, double rg, double thetag, double phig)
{
    double flux, slope_r = 0.0, slope_theta = 0.0, slope_phi = 0.0;
    
    if (rg - rC[i] != 0.0) slope_r = flux_limited_slope_r(vv, i, j, k, s)*(rg - rC[i]);
    if (thetag - thetaC[j] != 0.0) slope_theta = limited_slope_theta(vv, i, j, k, s)*(thetag - thetaC[j]);
    if (phig - phi[k] != 0.0) slope_phi = limited_slope_phi(vv, i, j, k, s)*(phig - phi[k]);

    flux = vv[k][j][i].fx[s] + slope_r + slope_theta + slope_phi;

    return flux;
}

/* limited slopes for Br, Btheta, and Bphi */
inline double limited_slope_Br_theta(Field ***xx, int i, int j, int k)
{
    double cc, dd;

    if (j == 0) {
        cc = (xx[k][j+1][i].fx[23]-xx[k][j][i].fx[23])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[23]-xx[(k+a3/2) % a3][0][i].fx[23])/(2.0*thetaC[j]);
    }
    else if (j > 0 && j < Nthm) {
        cc = (xx[k][j+1][i].fx[23]-xx[k][j][i].fx[23])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[23]-xx[k][j-1][i].fx[23])/(thetaC[j]-thetaC[j-1]);
    }
    else {
        cc = (xx[(k+a3/2) % a3][Nthm][i].fx[23]-xx[k][j][i].fx[23])/(2.0*thetaC[j]);
        dd = (xx[k][j][i].fx[23]-xx[k][j-1][i].fx[23])/(thetaC[j]-thetaC[j-1]);
    }

    return minmod(cc, dd);
}

inline double efd_limited_slope_r(Field ***xx, int i, int j, int k, int s)
{
    double aa=0.0, bb=0.0, rrb, rrt;
    int    Nrm1=Nr-1, Nrm2=Nr-2;

    if (i == 0) {
        aa = (xx[k][j][1].fx[s]-xx[k][j][i].fx[s])/(rC[1]-rC[i]);

        rrb=(rCm1-rC[1])/(rC[2]-rC[1]);  // rrb = (r_{-1} - r1)/(r2 - r1)
        bb = (xx[k][j][i].fx[s]-(xx[k][j][1].fx[s]+rrb*(xx[k][j][2].fx[s]-xx[k][j][1].fx[s])))/(rC[i]-rCm1);
    }
    else if (i > 0 && i < Nr) {
        aa = (xx[k][j][i+1].fx[s]-xx[k][j][i].fx[s])/(rC[i+1]-rC[i]);
        bb = (xx[k][j][i].fx[s]-xx[k][j][i-1].fx[s])/(rC[i]-rC[i-1]);
    }
    else {
        rrt=(rC[a1]-rC[Nrm1])/(rC[Nrm1]-rC[Nrm2]);
        aa = ((xx[k][j][Nrm1].fx[s]+rrt*(xx[k][j][Nrm1].fx[s]-xx[k][j][Nrm2].fx[s]))-xx[k][j][i].fx[s])
             /(rC[a1]-rC[i]);

        bb = (xx[k][j][i].fx[s]-xx[k][j][Nrm1].fx[s])/(rC[i]-rC[Nrm1]);
    }

    return minmod(aa, bb);
}

inline double efd_limited_slope_theta(Field ***xx, int i, int j, int k, int s)
{
    int kc;
    double cc, dd, sgn;

    if (j == 0) {
        kc = (k+a3/2) % a3;

        if (s == 1) sgn = -1.0; else sgn = 1.0;

        cc = (xx[k][j+1][i].fx[s]-xx[k][j][i].fx[s])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[s]-sgn*xx[kc][0][i].fx[s])/(2.0*thetaC[j]);
    }
    else if (j > 0 && j < Nth) {
        cc = (xx[k][j+1][i].fx[s]-xx[k][j][i].fx[s])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[s]-xx[k][j-1][i].fx[s])/(thetaC[j]-thetaC[j-1]);
    }
    else {
        kc = (k+a3/2) % a3;

        if (s == 1) sgn = -1.0; else sgn = 1.0;

        cc = (sgn*xx[kc][Nthm][i].fx[s]-xx[k][j][i].fx[s])/(2.0*thetaC[j]);
        dd = (xx[k][j][i].fx[s]-xx[k][j-1][i].fx[s])/(thetaC[j]-thetaC[j-1]);
    }

    return minmod(cc, dd);
}


/* reconstructed electric field within the cell (i, j, k) */
inline double reconstructed_efd(Field ***uu, int i, int j, int k, int s, double rg, double thetag, double phig)
{
    double efd, slope_r = 0.0, slope_theta = 0.0, slope_phi = 0.0;
    
    if ((rg - rC[i]) != 0.0) slope_r = efd_limited_slope_r(uu, i, j, k, s)*(rg - rC[i]);
    if ((thetag - thetaC[j]) != 0.0) slope_theta = efd_limited_slope_theta(uu, i, j, k, s)*(thetag - thetaC[j]);
    if ((phig - phi[k]) != 0.0) slope_phi = limited_slope_phi(uu, i, j, k, s)*(phig - phi[k]);

    efd = uu[k][j][i].fx[s] + slope_r + slope_theta + slope_phi;

    return efd;
}


inline double limited_slope_Br_phi(Field ***xx, int i, int j, int k)
{
    int km, kp;
    double ee, ff;

    if (k == 0) km = Np; else km = k-1;
    if (k < Np) kp = k+1; else kp = 0;

    ee = (xx[kp][j][i].fx[23]-xx[k][j][i].fx[23])/dph;
    ff = (xx[k][j][i].fx[23]-xx[km][j][i].fx[23])/dph;

    return minmod(ee, ff);
}

inline double limited_slope_Btheta_r(Field ***xx, int i, int j, int k)
{
    double cc, dd, rrb;

    if (i < Nr) cc = (xx[k][j][i+1].fx[24]-xx[k][j][i].fx[24])/(rfavg[i+1]-rfavg[i]);
    else cc = (xx[k][j][i].fx[24]-xx[k][j][i-1].fx[24])/(rfavg[i]-rfavg[i-1]);
    if (i == 0) {
        rrb=(rfavgm1-rC[1])/(rC[2]-rC[1]);

        //y1=xx[k][j][1].fx[24];
        //y2=xx[k][j][2].fx[24];
        //xx[k][j][i-1].fx[24] = y1+rrb*(y2-y1);
        dd = (xx[k][j][i].fx[24]-(xx[k][j][1].fx[24]+rrb*(xx[k][j][2].fx[24]-xx[k][j][1].fx[24])))
             /(rfavg[i]-rfavgm1);
    }
    else dd = (xx[k][j][i].fx[24]-xx[k][j][i-1].fx[24])/(rfavg[i]-rfavg[i-1]);

    return minmod(cc, dd);
}

inline double limited_slope_Btheta_phi(Field ***xx, int i, int j, int k)
{
    int km, kp;
    double ee, ff;

    if (k == 0) km = Np; else km = k-1;
    if (k < Np) kp = k+1; else kp = 0;

    ee = (xx[kp][j][i].fx[24]-xx[k][j][i].fx[24])/dph;
    ff = (xx[k][j][i].fx[24]-xx[km][j][i].fx[24])/dph;

    return minmod(ee, ff);
}
#include <iostream> 
inline double limited_slope_Bphi_r(Field ***xx, int i, int j, int k)
{
    double cc, dd, rrb;

    if (i < Nr) cc = (xx[k][j][i+1].fx[25]-xx[k][j][i].fx[25])/(rfavg[i+1]-rfavg[i]);
    else cc = (xx[k][j][i].fx[25]-xx[k][j][i-1].fx[25])/(rfavg[i]-rfavg[i-1]);
    if (i == 0) {
        rrb=(rfavgm1-rC[1])/(rC[2]-rC[1]);

        //y1=xx[k][j][1].fx[25];
        //y2=xx[k][j][2].fx[25];
        //xx[k][j][i-1].fx[25] = y1+rrb*(y2-y1);
        dd = (xx[k][j][i].fx[25]-(xx[k][j][1].fx[25]+rrb*(xx[k][j][2].fx[25]-xx[k][j][1].fx[25])))
             /(rfavg[i]-rfavgm1);
    }
    else dd = (xx[k][j][i].fx[25]-xx[k][j][i-1].fx[25])/(rfavg[i]-rfavg[i-1]);

    return minmod(cc, dd);
}

inline double limited_slope_Bphi_theta(Field ***xx, int i, int j, int k)
{
    double cc, dd;

    if (j == 0) {
        cc = (xx[k][j+1][i].fx[25]-xx[k][j][i].fx[25])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[25]-xx[(k+a3/2) % a3][0][i].fx[25])/(2.0*thetaC[j]);
    }
    else if (j > 0 && j < Nthm) {
        cc = (xx[k][j+1][i].fx[25]-xx[k][j][i].fx[25])/(thetaC[j+1]-thetaC[j]);
        dd = (xx[k][j][i].fx[25]-xx[k][j-1][i].fx[25])/(thetaC[j]-thetaC[j-1]);
    }
    else {
        cc = (xx[(k+a3/2) % a3][Nthm][i].fx[25]-xx[k][j][i].fx[25])/(2.0*thetaC[j]);
        dd = (xx[k][j][i].fx[25]-xx[k][j-1][i].fx[25])/(thetaC[j]-thetaC[j-1]);
    }

    return minmod(cc, dd);
}

/*-----------------------------------------------------------------------------------------*/
/*** reconstruction point Br, Btheta, Bphi with surface averaged <Br>, <Btheta>, <Bphi> ****/
/*-----------------------------------------------------------------------------------------*/

/* reconstructed Br within a cell (i, j, k). Used to calculate Br at faces of the cell */
inline double reconstructed_Br(Field ***xx, int i, int j, int k, double rg, double thetag, double phig)
{
     int ip = i+1;
     double Brijk;

     Brijk = 0.5*(xx[k][j][ip].fx[23] + xx[k][j][i].fx[23])
            +( (xx[k][j][ip].fx[23] - xx[k][j][i].fx[23])*(rg - rr[i])
              +( limited_slope_Br_theta(xx, ip, j, k)*(rg-rh[i])
                +limited_slope_Br_theta(xx, i, j, k)*(rh[ip]-rg))*(thetag - thetaC[j])
              +( limited_slope_Br_phi(xx, ip, j, k)*(rg-rh[i])
                +limited_slope_Br_theta(xx, i, j, k)*(rh[ip]-rg))*(phig - phi[k]))/dr;

    return Brijk;
}

/* reconstructed Btheta within a cell (i, j, k). Used to calculate Btheta at faces of the cell */
inline double reconstructed_Btheta(Field ***xx, int i, int j, int k, double rg, double thetag, double phig)
{
    int jp=j+1;
    double Bthetaijk;

    Bthetaijk = 0.5*(xx[k][j][i].fx[24] + xx[k][j][i].fx[24])
            +( (xx[k][jp][i].fx[24] - xx[k][j][i].fx[24])*(thetag - theta[j])
              +( limited_slope_Btheta_r(xx, i, jp, k)*(thetag-thetah[j])
                +limited_slope_Btheta_r(xx, i, j, k)*(thetah[j+1]-thetag))*(rg - rfavg[i])
              +( limited_slope_Btheta_phi(xx, i, jp, k)*(thetag-thetah[j])
                +limited_slope_Btheta_phi(xx, i, j, k)*(thetah[j+1]-thetag))*(phig - phi[k]))/dth;

    return Bthetaijk;
}

/* reconstructed Bphi within a cell (i, j, k). Used to calculate Bphi at faces of the cell */
inline double reconstructed_Bphi(Field ***xx, int i, int j, int k, double rg, double thetag, double phig)
{
    int kp;
    double Bphiijk;

    if (k < Np) kp = k+1; else kp = 0;

    Bphiijk = 0.5*(xx[kp][j][i].fx[25] + xx[k][j][i].fx[25])
            +( (xx[kp][j][i].fx[25] - xx[k][j][i].fx[25])*(phig - phi[k])
              +( limited_slope_Bphi_r(xx, i, j, kp)*(phig-phih[k])
                +limited_slope_Bphi_r(xx, i, j, k)*(phih[k+1]-phig))*(rg - rfavg[i])
              +( limited_slope_Bphi_theta(xx, i, j, kp)*(phig-phih[k])
                +limited_slope_Bphi_theta(xx, i, j, k)*(phih[k+1]-phig))*(thetag - theta[j]))/dph;

    return Bphiijk;
}

#endif