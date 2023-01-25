#ifndef OPERATORS_H
#define OPERATORS_H

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
    double aa, bb;

    if (i == 0) { //dN/dr|0 = dN/dr|1
        aa = (xx[k][j][1].fx[s]-xx[k][j][i].fx[s])/dr;
        bb = aa;
    }
    else if (i > 0 and i < Nr) {
        aa = (xx[k][j][i+1].fx[s]-xx[k][j][i].fx[s])/dr;
        bb = (xx[k][j][i].fx[s]-xx[k][j][i-1].fx[s])/dr;
    }
    else {
        bb = (xx[k][j][i].fx[s]-xx[k][j][i-1].fx[s])/dr;
        aa=bb;
    }

    return minmod(aa, bb);
}

inline double limited_slope_theta(Field ***xx, int i, int j, int k, int s)
{
    double cc, dd, sgn;

    if (s==8 || s==11 || s==14 || s==28 || s==32 || s==35) sgn = -1.0;
    else sgn = 1.0;

    //when calculating gradient in theta, here treat theta as increase
    if (j == 1) {
        cc = (xx[k][j+1][i].fx[s]-xx[k][j][i].fx[s])/dth;
        dd = (xx[k][j][i].fx[s]-sgn*xx[(k+a3/2) % a3][j][i].fx[s])/dth;
    }
    else if (j > 1 && j < Nthm) {
        cc = (xx[k][j+1][i].fx[s]-xx[k][j][i].fx[s])/dth;
        dd = (xx[k][j][i].fx[s]-xx[k][j-1][i].fx[s])/dth;
    }
    else {
        cc = (sgn*xx[(k+a3/2) % a3][j][i].fx[s]-xx[k][j][i].fx[s])/dth;
        dd = (xx[k][j][i].fx[s]-xx[k][j-1][i].fx[s])/dth;
    }

    return minmod(cc, dd);
}

inline double limited_slope_phi(Field ***xx, int i, int j, int k, int s)
{
    double ee, ff;

    ee = (xx[k+1][j][i].fx[s]-xx[k][j][i].fx[s])/dph;
    ff = (xx[k][j][i].fx[s]-xx[k-1][j][i].fx[s])/dph;

    return minmod(ee, ff);
}

inline double difference_r(Field ***xx, int i, int j, int k, int s)
{
    double yL, yR, yip, yim, slope_ri;

    yL = xx[k][j][i-1].fx[s] + limited_slope_r(xx, i-1, j, k, s)*(rh[i]-rr[i-1]);

    slope_ri=limited_slope_r(xx, i, j, k, s);
    yR = xx[k][j][i].fx[s] + slope_ri*(rh[i]-rr[i]);

    yim = 0.5*(yL + yR);

    yL = xx[k][j][i].fx[s] + slope_ri*(rh[i+1]-rr[i]);
    yR = xx[k][j][i+1].fx[s] + limited_slope_r(xx, i+1, j, k, s)*(rh[i+1]-rr[i+1]);

    yip = 0.5*(yL + yR);

    double dydr = (yip - yim)/dr;

    return dydr;
}

inline double difference_theta(Field ***xx, int i, int j, int k, int s)
{
    double yL, yR, yjp, yjm, slope_thetaj;
    int    kc, jm, jp;

    if (j == 1) {
        kc = (k+a3/2) % a3; 
        if (s==8 || s==11 || s==14 || s==28 || s==32 || s==35)
            yL = -(xx[kc][j][i].fx[s] - limited_slope_theta(xx, i, j, kc, s)*(thetah[j]-theta[j]));
        else yL = xx[kc][j][i].fx[s] + limited_slope_theta(xx, i, j, kc, s)*(thetah[j]-theta[j]);
    }
    else {
        jm=j-1;
        yL = xx[k][jm][i].fx[s] + limited_slope_theta(xx, i, jm, k, s)*(thetah[j]-theta[jm]);
    }

    slope_thetaj=limited_slope_theta(xx, i, j, k, s);
    yR = xx[k][j][i].fx[s] + slope_thetaj*(thetah[j]-theta[j]);

    yjm = 0.5*(yL + yR);

    jp=j+1;
    yL = xx[k][j][i].fx[s] + slope_thetaj*(thetah[jp]-theta[j]);

    if (j == Nthm) {
        kc = (k+a3/2) % a3;
        if (s==8 || s==11 || s==14 || s==28 || s==32 || s==35)
            yR = -(xx[kc][j][i].fx[s] - limited_slope_theta(xx, i, j, kc, s)*(thetah[jp]-theta[j]));
        else yR = xx[kc][j][i].fx[s] + limited_slope_theta(xx, i, j, kc, s)*(thetah[jp]-theta[j]);
    }
    else {
        yR = xx[k][jp][i].fx[s] + limited_slope_theta(xx, i, jp, k, s)*(thetah[jp]-theta[jp]);
    }

    yjp = 0.5*(yL + yR);

    double dydtheta = (yjp - yjm)/dth;

    return dydtheta;
}

inline double difference_phi(Field ***xx, int i, int j, int k, int s)
{
    double yL, yR, ykp, ykm, slope_phk;
    int    km=k-1, kp=k+1;

    if (k==0) yL = xx[Np][j][i].fx[s] + limited_slope_phi(xx, i, j, Np, s)*(phih[Np+1]-phi[Np]);
    else yL = xx[km][j][i].fx[s] + limited_slope_phi(xx, i, j, km, s)*(phih[k]-phi[km]);

    slope_phk=limited_slope_phi(xx, i, j, k, s);
    yR = xx[k][j][i].fx[s] + slope_phk*(phih[k]-phi[k]);

    ykm = 0.5*(yL + yR);

    yL = xx[k][j][i].fx[s] + slope_phk*(phih[kp]-phi[k]);

    if (k == Np) yR = xx[0][j][i].fx[s] + limited_slope_phi(xx, i, j, 0, s)*(phih[0]-phi[0]);
    else yR = xx[kp][j][i].fx[s] + limited_slope_phi(xx, i, j, kp, s)*(phih[kp]-phi[kp]);

    ykp = 0.5*(yL + yR);

    double dydphi = (ykp - ykm)/dph;

    return dydphi;
}

inline vector3D gradient(Field ***xx, int i, int j, int k, int yj, int xi, int s)
{
    vector3D y;

    double dydr = difference_r(xx, i, j, k, s);     //d/d_{r}
    double dydt = difference_theta(xx, i, j, k, s); //d/d_{theta}
    double dydp = difference_phi(xx, i, j, k, s);   //d/d_{phi}

    y.r = dydr; y.t = dydt/rr[i]; y.p = dydp/rsin[yj][xi];

    return y;
}

inline double divergence(Field ***xx, int i, int j, int k, int yj, int xi, int s)
{
    double dyr_dr = difference_r(xx, i, j, k, s);       //d/d_{r}
    double dyt_dt = difference_theta(xx, i, j, k, s+1); //d/d_{theta}
    double dyp_dp = difference_phi(xx, i, j, k, s+2);   //d/d_{phi}

    return 2.0*xx[k][j][i].fx[s]/rr[i]+xx[k][j][i].fx[s+1]*cot_div_r[yj][xi]
          +dyr_dr + dyt_dt/rr[i] + dyp_dp/rsin[yj][xi]; 
}

#endif