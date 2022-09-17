#include "reconstruction.h"
#include <iostream>
#include <fstream>

inline void max_speed_r_face(Field ***xx, Field ***zz, int i, int j, int k, double rg, double thetag,
    double phig, double a_imjk[])
{
    double rhoim_mjk = 0.0, rhoim_pjk = 0.0, rhonim_mjk = 0.0, rhonim_pjk = 0.0;
    double Br0, Btheta0, BB2, Br, Btheta, Bphi, rhoeff_mu0, ratio;
    const double five3rd=5.0/3.0;

    int im = i-1;

    //background magnetic field at (rg, thetag, phig)
    ambient_mfd(rg, thetag, Br0, Btheta0);

    //mass density in kg cm^{-3}
    for (int s = 0; s < sl; s++) {
        rhoim_mjk += reconstructed(xx, im, j, k, s, rg, thetag, phig)*ms[s];
        rhoim_pjk += reconstructed(xx, i, j, k, s, rg, thetag, phig)*ms[s];

        rhonim_mjk += reconstructed(xx, im, j, k, 12+s, rg, thetag, phig)*ms[s];
        rhonim_pjk += reconstructed(xx, i, j, k, 12+s, rg, thetag, phig)*ms[s];
    }

/* evaluate inteface speed a^{+}_{r_{i-1/2}, theta_{face_avg}[j]}, phi_k */
    //pressure in Joule cm^{-3}
    double P = reconstructed(xx, im, j, k, 10, rg, thetag, phig)
              +reconstructed(xx, im, j, k, 11, rg, thetag, phig);

    // squared sound at (im-, thetaC[j], k) in (m/s)^2
    double Cs2 = five3rd*P/rhoim_mjk;

    Br = reconstructed_Br(xx, im, j, k, rg, thetag, phig) + Br0;
    Btheta = reconstructed_Btheta(xx, im, j, k, rg, thetag, phig) + Btheta0;
    Bphi = reconstructed_Bphi(xx, im, j, k, rg, thetag, phig);
    BB2 = Br*Br+Btheta*Btheta+Bphi*Bphi;

    //mass density in kg/cm^3. So multplying by 1e6 to convert to kg/m^3
    ratio=reconstructed(zz, im, j, k, 23, rg, thetag, phig);
    rhoeff_mu0=1.0e6*(rhoim_mjk+ratio*rhonim_mjk)*mu0;

    //Alfven speed square at (im-, thetaC[j], k) in (m.s)^2. Fast speed along r direction in m/s
    double CA2 = BB2/rhoeff_mu0;
    double Cfr=sqrt(0.5*(Cs2+CA2+sqrt((Cs2+CA2)*(Cs2+CA2)-4.0*Cs2*(Br*Br/rhoeff_mu0))));
    double ur = reconstructed(xx, im, j, k, 7, rg, thetag, phig)/rhoim_mjk;

    double aa = ur + Cfr;
    double bb = ur - Cfr;

    // squared sound & Alfven speeds at (im+, thetaC[j], k)
    P = reconstructed(xx, i, j, k, 10, rg, thetag, phig)
       +reconstructed(xx, i, j, k, 11, rg, thetag, phig);
    Cs2 = five3rd*P/rhoim_pjk;

    if (i<Nr) Br = reconstructed_Br(xx, i, j, k, rg, thetag, phig) + Br0;
    else Br = xx[k][j][i].fx[23] + Br0;
    Btheta = reconstructed_Btheta(xx, i, j, k, rg, thetag, phig) + Btheta0;
    Bphi = reconstructed_Bphi(xx, i, j, k, rg, thetag, phig);
    BB2 = Br*Br+Btheta*Btheta+Bphi*Bphi;

    ratio=reconstructed(zz, i, j, k, 23, rg, thetag, phig);
    rhoeff_mu0=1.0e6*(rhoim_pjk+ratio*rhonim_pjk)*mu0;

    CA2 = BB2/rhoeff_mu0;
    Cfr=sqrt(0.5*(Cs2+CA2+sqrt((Cs2+CA2)*(Cs2+CA2)-4.0*Cs2*(Br*Br/(rhoeff_mu0)))));
    ur = reconstructed(xx, i, j, k, 7, rg, thetag, phig)/rhoim_pjk;

    double cc= ur + Cfr;
    double dd= ur - Cfr;

    a_imjk[0] = max(aa, cc);
    a_imjk[0] = max(a_imjk[0], 0.0); //velocity in m/s

    a_imjk[1] = min(bb, dd);
    a_imjk[1] = min(a_imjk[1], 0.0); //velocity in m/s
}

inline void max_speed_theta_face(Field ***xx, Field ***zz, int i, int j, int k, double rg, double thetag,
    double phig, double b_ijmk[])
{
    double rhoijm_mk = 0.0, rhoijm_pk = 0.0, rhonijm_mk = 0.0, rhonijm_pk = 0.0;
    double Br0, Btheta0, BB2, Br, Btheta, Bphi, rhoeff_mu0, ratio;
    const double five3rd=5.0/3.0;

    int jm, jj, k_0_c, k_Nth_c;

    if (j == 1) {k_0_c=(k+a3/2) % a3; jm=j;} else {k_0_c=k; jm=j-1;}
    if (j < Nth) {k_Nth_c=k; jj=j;} else {k_Nth_c=(k+a3/2) % a3; jj=Nthm;}

    ambient_mfd(rg, thetag, Br0, Btheta0);

    for (int s = 0; s < sl; s++) {
        rhoijm_mk += reconstructed(xx, i, jm, k_0_c, s, rg, thetag, phig)*ms[s];
        rhoijm_pk += reconstructed(xx, i, jj, k_Nth_c, s, rg, thetag, phig)*ms[s];

        rhonijm_mk += reconstructed(xx, i, jm, k_0_c, 12+s, rg, thetag, phig)*ms[s];
        rhonijm_pk += reconstructed(xx, i, jj, k_Nth_c, 12+s, rg, thetag, phig)*ms[s];
    }

    // sound & Alfven speeds at (rfavg, jm-, k) in (m/s)^2
    double P = reconstructed(xx, i, jm, k_0_c, 10, rg, thetag, phig)
              +reconstructed(xx, i, jm, k_0_c, 11, rg, thetag, phig);
    double Cs2 = five3rd*P/rhoijm_mk;

    Br = reconstructed_Br(xx, i, jm, k_0_c, rg, thetag, phig) + Br0;
    Btheta = reconstructed_Btheta(xx, i, jm, k_0_c, rg, thetag, phig) + Btheta0;
    Bphi = reconstructed_Bphi(xx, i, jm, k_0_c, rg, thetag, phig);
    BB2 = Br*Br+Btheta*Btheta+Bphi*Bphi;

    ratio=reconstructed(zz, i, jm, k_0_c, 23, rg, thetag, phig);
    rhoeff_mu0=1.0e6*(rhoijm_mk+ratio*rhonijm_mk)*mu0;

    double CA2 = BB2/rhoeff_mu0;
    double Cftheta=sqrt(0.5*(Cs2+CA2+sqrt((Cs2+CA2)*(Cs2+CA2)-4.0*Cs2*(Btheta*Btheta/rhoeff_mu0))));
    double utheta = reconstructed(xx, i, jm, k_0_c, 8, rg, thetag, phig)/rhoijm_mk;

    double aa= utheta + Cftheta;
    double bb= utheta - Cftheta;

    // sound & Alfven speeds at (rfavg, jm+, k)
    P = reconstructed(xx, i, jj, k_Nth_c, 10, rg, thetag, phig)
       +reconstructed(xx, i, jj, k_Nth_c, 11, rg, thetag, phig);
    Cs2 = five3rd*P/rhoijm_pk;

    Br = reconstructed_Br(xx, i, jj, k_Nth_c, rg, thetag, phig) + Br0;
    Btheta = reconstructed_Btheta(xx, i, jj, k_Nth_c, rg, thetag, phig) + Btheta0;
    Bphi = reconstructed_Bphi(xx, i, jj, k_Nth_c, rg, thetag, phig);
    BB2 = Br*Br+Btheta*Btheta+Bphi*Bphi;

    ratio=reconstructed(zz, i, jj, k_Nth_c, 23, rg, thetag, phig);
    rhoeff_mu0=1.0e6*(rhoijm_pk+ratio*rhonijm_pk)*mu0;

    CA2 = BB2/rhoeff_mu0;
    Cftheta=sqrt(0.5*(Cs2+CA2+sqrt((Cs2+CA2)*(Cs2+CA2)-4.0*Cs2*(Btheta*Btheta/rhoeff_mu0))));
    utheta = reconstructed(xx, i, jj, k_Nth_c, 8, rg, thetag, phig)/rhoijm_pk;

    double cc= utheta + Cftheta;
    double dd= utheta - Cftheta;

    b_ijmk[0] = max(aa, cc);
    b_ijmk[0] = max(b_ijmk[0], 0.0); //velocity in m/s

    b_ijmk[1] = min(bb, dd);
    b_ijmk[1] = min(b_ijmk[1], 0.0);
}

// sound & Alfven speeds at (rfavg, j, km-)
inline void max_speed_phi_face(Field ***xx, Field ***zz, int i, int j, int k, double rg, double thetag,
    double phig, double c_ijkm[])
{
    double rhoijkm_m = 0.0, rhoijkm_p = 0.0, rhonijkm_m = 0.0, rhonijkm_p = 0.0;
    double Br0, Btheta0, BB2, Br, Btheta, Bphi, phig_m, rhoeff_mu0, ratio;
    const double five3rd=5.0/3.0;

    int km;
    if (k == 0) {km = Np; phig_m=pi2-phig;} else {km = k-1; phig_m=phig;}

    ambient_mfd(rg, thetag, Br0, Btheta0);

    for (int s = 0; s < sl; s++) {
        rhoijkm_m += reconstructed(xx, i, j, km, s, rg, thetag, phig_m)*ms[s];
        rhoijkm_p += reconstructed(xx, i, j, k, s, rg, thetag, phig)*ms[s];

        rhonijkm_m += reconstructed(xx, i, j, km, 12+s, rg, thetag, phig_m)*ms[s];
        rhonijkm_p += reconstructed(xx, i, j, k, 12+s, rg, thetag, phig)*ms[s];
    }

    double P = reconstructed(xx, i, j, km, 10, rg, thetag, phig_m)
              +reconstructed(xx, i, j, km, 11, rg, thetag, phig_m);
    double Cs2 = five3rd*P/rhoijkm_m;

    Br = reconstructed_Br(xx, i, j, km, rg, thetag, phig_m) + Br0;
    Btheta = reconstructed_Btheta(xx, i, j, km, rg, thetag, phig_m) + Btheta0;
    Bphi = reconstructed_Bphi(xx, i, j, km, rg, thetag, phig_m);
    BB2 = Br*Br+Btheta*Btheta+Bphi*Bphi;

    ratio=reconstructed(zz, i, j, km, 23, rg, thetag, phig_m);
    rhoeff_mu0=1.0e6*(rhoijkm_m+ratio*rhonijkm_m)*mu0;

    double CA2 = BB2/rhoeff_mu0;
    double Cfphi=sqrt(0.5*(Cs2+CA2+sqrt((Cs2+CA2)*(Cs2+CA2)-4.0*Cs2*(Bphi*Bphi/rhoeff_mu0))));
    double uphi = reconstructed(xx, i, j, km, 9, rg, thetag, phig_m)/rhoijkm_m;

    double aa= uphi + Cfphi;
    double bb= uphi - Cfphi;

    // sound & Alfven speeds at (rfavg, jm+, k)
    P = reconstructed(xx, i, j, k, 10, rg, thetag, phig)
       +reconstructed(xx, i, j, k, 11, rg, thetag, phig);
    Cs2 = five3rd*P/rhoijkm_p;

    Br = reconstructed_Br(xx, i, j, k, rg, thetag, phig) + Br0;
    Btheta = reconstructed_Btheta(xx, i, j, k, rg, thetag, phig) + Btheta0;
    Bphi = reconstructed_Bphi(xx, i, j, k, rg, thetag, phig);
    BB2 = Br*Br+Btheta*Btheta+Bphi*Bphi;

    ratio=reconstructed(zz, i, j, k, 23, rg, thetag, phig);
    rhoeff_mu0=1.0e6*(rhoijkm_p+ratio*rhonijkm_p)*mu0;

    CA2 = BB2/rhoeff_mu0;
    Cfphi=sqrt(0.5*(Cs2+CA2+sqrt((Cs2+CA2)*(Cs2+CA2)-4.0*Cs2*(Bphi*Bphi/rhoeff_mu0))));
    uphi = reconstructed(xx, i, j, k, 9, rg, thetag, phig)/rhoijkm_p;

    double cc= uphi + Cfphi;
    double dd= uphi - Cfphi;

    c_ijkm[0] = max(aa, cc);
    c_ijkm[0] = max(c_ijkm[0], 0.0); //velocity in m/s

    c_ijkm[1] = min(bb, dd);
    c_ijkm[1] = min(c_ijkm[1], 0.0);
}
