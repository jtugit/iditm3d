#include "reconstruction.h"
#include <iostream>

inline void neu_max_speed_r_face(Field ***xx, Field ***uu, int i, int j, int k, double an_imjk[])
{
    double rhonim_mjk = 0.0, rhonim_pjk = 0.0;
    const double five3rd=5.0/3.0;

    int im = i-1;

    for (int s = 0; s < sm; s++) {
        rhonim_mjk += reconstructed(xx, im, j, k, 12+s, rh[i], thetaC[j], phi[k])*ms[s];
        rhonim_pjk += reconstructed(xx, i, j, k, 12+s, rh[i], thetaC[j], phi[k])*ms[s];
    }

/* evaluate inteface speed a^{+}_{r_{i-1/2}, theta_{face_avg}[j], phi_k} */
    double Pn = reconstructed(xx, im, j, k, 22, rh[i], thetaC[j], phi[k]);

    // sound & Alfven speeds at (im-, thetaC[j], k)
    double Cn = sqrt(five3rd*Pn/rhonim_mjk);

    double unr = reconstructed(uu, im, j, k, 19, rh[i], thetaC[j], phi[k]);

    double aa = unr + Cn;
    double bb = unr - Cn;

    // sound & Alfven speeds at (im+, thetaC[j], k)
    Pn = reconstructed(xx, i, j, k, 22, rh[i], thetaC[j], phi[k]);
    Cn = sqrt(five3rd*Pn/rhonim_pjk);

    unr = reconstructed(uu, i, j, k, 19, rh[i], thetaC[j], phi[k]);

    double cc= unr + Cn;
    double dd= unr - Cn;

    an_imjk[0] = max(aa, cc);
    an_imjk[0] = max(an_imjk[0], 0.0);

    an_imjk[1] = min(bb, dd);
    an_imjk[1] = min(an_imjk[1], 0.0);   
}

inline void neu_max_speed_theta_face(Field ***xx, Field ***uu, int i, int j, int k, double bn_ijmk[])
{
    double rhonijm_mk = 0.0, rhonijm_pk = 0.0, sgn;
    const double five3rd=5.0/3.0;

    int jm, kc;

    if (j == 0) {jm = 0; kc = (k+a3/2) % a3; sgn=-1.0;}
    else {jm = j-1; kc = k; sgn=1.0;}

    for (int s = 0; s < sm; s++) {
        rhonijm_mk += reconstructed(xx, i, jm, kc, 12+s, rfavg[i], thetah[j], phi[kc])*ms[s];
        rhonijm_pk += reconstructed(xx, i, j, k, 12+s, rfavg[i], thetah[j], phi[k])*ms[s];
    }

/* evaluate inteface speed a^{+}_{r_{i-1/2}, theta_{face_avg}[j], phi_k} */
    double Pn = reconstructed(xx, i, jm, kc, 22, rfavg[i], thetah[j], phi[kc]);

    // sound & Alfven speeds at (im-, thetaC[j], k)
    double Cn = sqrt(five3rd*Pn/rhonijm_mk);

    double untheta = sgn*reconstructed(uu, i, jm, kc, 20, rfavg[i], thetah[j], phi[kc]);

    double aa = untheta + Cn;
    double bb = untheta - Cn;

    // sound & Alfven speeds at (im+, thetaC[j], k)
    Pn = reconstructed(xx, i, j, k, 22, rfavg[i], thetah[j], phi[k]);
    Cn = sqrt(five3rd*Pn/rhonijm_pk);

    untheta = reconstructed(xx, i, j, k, 20, rfavg[i], thetah[j], phi[k]);

    double cc= untheta + Cn;
    double dd= untheta - Cn;

    bn_ijmk[0] = max(aa, cc);
    bn_ijmk[0] = max(bn_ijmk[0], 0.0);

    bn_ijmk[1] = min(bb, dd);
    bn_ijmk[1] = min(bn_ijmk[1], 0.0);   
}

inline void neu_max_speed_phi_face(Field ***xx, Field ***uu, int i, int j, int k, double cn_ijkm[])
{
    double rhonijkm_m = 0.0, rhonijkm_p = 0.0;
    const double five3rd=5.0/3.0;

    int km, kprime;

    if (k == 0) {km = Np; kprime = Np;}
    else {km = k-1; kprime = k;}

    for (int s = 0; s < sm; s++) {
        rhonijkm_m += reconstructed(xx, i, j, km, 12+s, rfavg[i], theta[j], phih[kprime])*ms[s];
        rhonijkm_p += reconstructed(xx, i, j, k, 12+s, rfavg[i], theta[j], phih[k])*ms[s];
    }

/* evaluate inteface speed a^{+}_{r_{i-1/2}, theta_{face_avg}[j], phi_k} */
    double Pn = reconstructed(xx, i, j, km, 22, rfavg[i], thetah[j], phih[kprime]);

    // sound & Alfven speeds at (im-, thetaC[j], k)
    double Cn = sqrt(five3rd*Pn/rhonijkm_m);

    double unphi = reconstructed(uu, i, j, km, 21, rfavg[i], theta[j], phih[kprime]);

    double aa = unphi + Cn;
    double bb = unphi - Cn;

    // sound & Alfven speeds at (im+, thetaC[j], k)
    Pn = reconstructed(xx, i, j, k, 22, rfavg[i], theta[j], phih[k]);
    Cn = sqrt(five3rd*Pn/rhonijkm_p);

    unphi = reconstructed(xx, i, j, k, 21, rfavg[i], theta[j], phih[k]);

    double cc= unphi + Cn;
    double dd= unphi - Cn;

    cn_ijkm[0] = max(aa, cc);
    cn_ijkm[0] = max(cn_ijkm[0], 0.0);

    cn_ijkm[1] = min(bb, dd);
    cn_ijkm[1] = min(cn_ijkm[1], 0.0);   
}

/* ------ r-face fluxes at intefaces (i-1/2, j, k) of the cell i ---------------*/
inline void fluxes_r(Field ***xx, Field ***uu, Field ***vv, int i, int j, int k, double flux_r[])
{
    int im=i-1, s;

    double fluxim_mjk, fluxim_pjk;
    double Uim_mjk, Uim_pjk;
    double drh_rCm=rh[i]-rC[im], drh_rC=rh[i]-rC[i];

    //for plasma equations, terms associated with fluxes are partially treated explicitly
    for (s = 0; s < 12; s++) {
        //first calcuated fluxes on -/+ sides of i-1/2
        fluxim_mjk = vv[k][j][i].fx[s]+flux_limited_slope_r(vv, im, j, k, s)*drh_rCm;
        fluxim_pjk = vv[k][j][i].fx[s]+flux_limited_slope_r(vv, i, j, k, s)*drh_rC;

        flux_r[s] = 0.5*(fluxim_mjk + fluxim_pjk);
    }

//********* flux for n_q (q=0, 1, 2... sm) equations *****************
    double an_imjk[2];
    neu_max_speed_r_face(xx, uu, i, j, k, an_imjk);
    double an_imjk_multi=an_imjk[0]*an_imjk[1], an_imjk_subtr=an_imjk[0]-an_imjk[1];

    //for neutral equations terms associated with fluxes are completely treated explicitly
    for (s = 12; s < 23; s++) {
    //reconstructed conservative variables for neutrals at im- and im+ sides of the interface im=i-1/2
        Uim_mjk = xx[k][j][i].fx[s]+limited_slope_r(xx, im, j, k, s)*drh_rCm;
        Uim_pjk = xx[k][j][i].fx[s]+limited_slope_r(xx, i, j, k, s)*drh_rC;

    //**** reconstructed neutral fluxes at interface (i-1/2, j, k) of the cell i *****************
        //first calcuated fluxes on -/+ sides of i-1/2
        fluxim_mjk = vv[k][j][i].fx[s]+flux_limited_slope_r(vv, im, j, k, s)*drh_rCm;
        fluxim_pjk = vv[k][j][i].fx[s]+flux_limited_slope_r(vv, i, j, k, s)*drh_rC;

        flux_r[s] = ( (an_imjk[0]*fluxim_mjk - an_imjk[1]*fluxim_pjk)
                     -an_imjk_multi*(Uim_mjk-Uim_pjk))/an_imjk_subtr;
    }
}

/* ------ theta-face fluxes at intefaces (i, j-1/2, k) i ------*/
inline void fluxes_theta(Field ***xx, Field ***uu, Field ***ww, int i, int j, int k, double flux_theta[])
{
    int jm, kc, s;

    double fluxijm_mk, fluxijm_pk;
    double Uijm_mk, Uijm_pk;
    double drfavg_rC=rfavg[i]-rC[i], dth_thCm;
    double dth_thC=thetah[j]-thetaC[j];

    if (j == 0) {
        jm = 0; kc = (k+a3/2) % a3;
        dth_thCm=2.0*thetah[j];
    }
    else {
        jm = j-1; kc = k;
        dth_thCm=thetah[j]-thetaC[jm];
    }

    for (s = 0; s < 12; s++) {
        fluxijm_mk = ww[k][j][i].fx[s]+flux_limited_slope_r(ww, i, jm, kc, s)*drfavg_rC
                    +limited_slope_theta(ww, i, jm, kc, s)*dth_thCm;
        fluxijm_pk = ww[k][j][i].fx[s]+flux_limited_slope_r(ww, i, j, k, s)*drfavg_rC
                    +limited_slope_theta(ww, i, j, k, s)*dth_thC;

        flux_theta[s] = 0.5*(fluxijm_mk + fluxijm_pk);
    }

//********* flux for n_q (q=0, 1, 2... sm) equations *****************
    double bn_ijmk[2];
    neu_max_speed_theta_face(xx, uu, i, j, k, bn_ijmk);
    double bn_ijmk_multi = bn_ijmk[0]*bn_ijmk[1], bn_ijmk_subtr = bn_ijmk[0] - bn_ijmk[1];

    for (s = 12; s < 23; s++) {
        Uijm_mk = xx[k][j][i].fx[s]+limited_slope_r(xx, i, jm, kc, s)*drfavg_rC
                 +limited_slope_theta(xx, i, jm, kc, s)*dth_thCm;
        Uijm_pk = xx[k][j][i].fx[s]+limited_slope_r(xx, i, j, k, s)*drfavg_rC
                 +limited_slope_theta(xx, i, j, kc, s)*dth_thC;

        fluxijm_mk = ww[k][j][i].fx[s]+flux_limited_slope_r(ww, i, jm, kc, s)*drfavg_rC
                    +limited_slope_theta(ww, i, jm, kc, s)*dth_thCm;
        fluxijm_pk = ww[k][j][i].fx[s]+flux_limited_slope_r(ww, i, j, k, s)*drfavg_rC
                    +limited_slope_theta(ww, i, j, k, s)*dth_thC;

        flux_theta[s] = ( (bn_ijmk[0]*fluxijm_mk - bn_ijmk[1]*fluxijm_pk)
                         -bn_ijmk_multi*(Uijm_mk - Uijm_pk))/bn_ijmk_subtr;
    }
}

/* ------ phi-face fluxes at intefaces (i, j-1/2, k) i ------*/
inline void fluxes_phi(Field ***xx, Field ***uu, Field ***zz, int i, int j, int k, double flux_phi[])
{
    int km, s;

    double fluxijkm_m, fluxijkm_p;
    double Uijkm_m, Uijkm_p;
    double drfavg_rC=rfavg[i]-rC[i], dthh_thC=thetah[j]-thetaC[j], dphim_phikm, dphim_phi=phih[k]-phi[k];
                                                           //    phi[Np] phih[Np+1]  phi[0]
    if (k == 0) {km = Np; dphim_phikm = phih[Np+1]-phi[km];} // |      o      |          o   
    else {km = k-1; dphim_phikm = phih[k]-phi[km];}          //      k=Np phi_Np+1/2    k=0

    for (s = 0; s < 12; s++) {
        fluxijkm_m = zz[km][j][i].fx[s]+flux_limited_slope_r(zz, i, j, km, s)*drfavg_rC
                    +limited_slope_theta(zz, i, j, km, s)*dthh_thC
                    +limited_slope_phi(zz, i, j, km, s)*dphim_phikm;
        fluxijkm_p = zz[k][j][i].fx[s]+flux_limited_slope_r(zz, i, j, k, s)*drfavg_rC
                    +limited_slope_theta(zz, i, j, k, s)*dthh_thC
                    +limited_slope_phi(zz, i, j, k, s)*dphim_phi;

        flux_phi[s] = 0.5*(fluxijkm_m + fluxijkm_p);
    }

//********* flux for n_q (q=0, 1, 2... sm) equations *****************
    double cn_ijkm[2];
    neu_max_speed_theta_face(xx, uu, i, j, k, cn_ijkm);
    double cn_ijkm_multi = cn_ijkm[0]*cn_ijkm[1], cn_ijkm_subtr = cn_ijkm[0]-cn_ijkm[1];

    for (s = 12; s < 23; s++) {
        Uijkm_m = xx[k][j][i].fx[s]+limited_slope_r(xx, i, j, km, s)*drfavg_rC
                 +limited_slope_theta(xx, i, j, km, s)*dthh_thC
                 +limited_slope_phi(xx, i, j, km, s)*dphim_phikm;
        Uijkm_p = xx[k][j][i].fx[s]+limited_slope_r(xx, i, j, k, s)*drfavg_rC
                 +limited_slope_theta(xx, i, j, k, s)*dthh_thC
                 +limited_slope_phi(xx, i, j, k, s)*dphim_phi;

        fluxijkm_m = zz[k][j][i].fx[s]+flux_limited_slope_r(zz, i, j, km, s)*drfavg_rC
                    +limited_slope_theta(zz, i, j, km, s)*dthh_thC
                    +limited_slope_phi(zz, i, j, km, s)*dphim_phikm;
        fluxijkm_p = zz[k][j][i].fx[s]+flux_limited_slope_r(zz, i, j, k, s)*drfavg_rC
                    +limited_slope_theta(zz, i, j, k, s)*dthh_thC
                    +limited_slope_phi(zz, i, j, k, s)*dphim_phi;

        flux_phi[s] = ( (cn_ijkm[0]*fluxijkm_m - cn_ijkm[1]*fluxijkm_p)
                       -cn_ijkm_multi*(Uijkm_m - Uijkm_p))/cn_ijkm_subtr;
    }
}