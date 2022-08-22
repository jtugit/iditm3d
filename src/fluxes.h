#include "max_mhd_speeds.h"
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
    double Cn = five3rd*Pn/rhonim_mjk;

    double unr = reconstructed(uu, im, j, k, 19, rh[i], thetaC[j], phi[k]);

    double aa = unr + Cn;
    double bb = unr - Cn;

    // sound & Alfven speeds at (im+, thetaC[j], k)
    Pn = reconstructed(xx, i, j, k, 22, rh[i], thetaC[j], phi[k]);
    Cn = five3rd*Pn/rhonim_pjk;

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
    double Cn = five3rd*Pn/rhonijm_mk;

    double untheta = sgn*reconstructed(uu, i, jm, kc, 20, rfavg[i], thetah[j], phi[kc]);

    double aa = untheta + Cn;
    double bb = untheta - Cn;

    // sound & Alfven speeds at (im+, thetaC[j], k)
    Pn = reconstructed(xx, i, j, k, 22, rfavg[i], thetah[j], phi[k]);
    Cn = five3rd*Pn/rhonijm_pk;

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
    double Cn = five3rd*Pn/rhonijkm_m;

    double unphi = reconstructed(uu, i, j, km, 21, rfavg[i], theta[j], phih[kprime]);

    double aa = unphi + Cn;
    double bb = unphi - Cn;

    // sound & Alfven speeds at (im+, thetaC[j], k)
    Pn = reconstructed(xx, i, j, k, 22, rfavg[i], theta[j], phih[k]);
    Cn = five3rd*Pn/rhonijkm_p;

    unphi = reconstructed(xx, i, j, k, 21, rfavg[i], theta[j], phih[k]);

    double cc= unphi + Cn;
    double dd= unphi - Cn;

    cn_ijkm[0] = max(aa, cc);
    cn_ijkm[0] = max(cn_ijkm[0], 0.0);

    cn_ijkm[1] = min(bb, dd);
    cn_ijkm[1] = min(cn_ijkm[1], 0.0);   
}

/* ------ r-face fluxes at intefaces (i-1/2, j, k) of the cell i ---------------*/
inline void fluxes_r(Field ***xx, Field ***uu, Field ***vv, Field ***zz, int i, int j, int k, double flux_r[])
{
    int im=i-1, s;

    double fluxim_mjk, fluxim_pjk;
    double Uim_mjk[nvar], Uim_pjk[nvar];
    double a_imjk[2];

    max_speed_r_face(xx, zz, i, j, k, rh[i], thetaC[j], phi[k], a_imjk);
    double a_imjk_multi = a_imjk[0]*a_imjk[1], a_imjk_subtr = a_imjk[0]-a_imjk[1];
    if (a_imjk[0] >5000.0e3 || fabs(a_imjk[1]) > 5000.0e3) std::cout<<"a_imjk "<<a_imjk[0]<<" "<<i<<" "<<j<<" "<<k<<endl; 

    for (s = 0; s < 12; s++) {
    //reconstructed conservative variables for plasma at im- and im+ sides of the interface im=i-1/2
        Uim_mjk[s] = reconstructed(xx, im, j, k, s, rh[i], thetaC[j], phi[k]);
        Uim_pjk[s] = reconstructed(xx, i, j, k, s, rh[i], thetaC[j], phi[k]);

    //**** reconstructed ion fluxes at interface (i-1/2, j, k) of the cell i *****************
        //first calcuated fluxes on -/+ sides of i-1/2
        fluxim_mjk = reconstructed_flux(vv, im, j, k, s, rh[i], thetaC[j], phi[k]);
        fluxim_pjk = reconstructed_flux(vv, i, j, k, s, rh[i], thetaC[j], phi[k]);

        flux_r[s] = ( (a_imjk[0]*fluxim_mjk - a_imjk[1]*fluxim_pjk)
                     -a_imjk_multi*(Uim_mjk[s] - Uim_pjk[s]))/a_imjk_subtr;
    }

//********* flux for n_q (q=0, 1, 2... sm) equations *****************
    double an_imjk[2];
    neu_max_speed_r_face(xx, uu, i, j, k, an_imjk);
    double an_imjk_multi=an_imjk[0]*an_imjk[1], an_imjk_subtr=an_imjk[0]-an_imjk[1];

    for (s = 12; s < 23; s++) {
    //reconstructed conservative variables for neutrals at im- and im+ sides of the interface im=i-1/2
        Uim_mjk[s] = reconstructed(xx, im, j, k, s, rh[i], thetaC[j], phi[k]);
        Uim_pjk[s] = reconstructed(xx, i, j, k, s, rh[i], thetaC[j], phi[k]);

    //**** reconstructed neutral fluxes at interface (i-1/2, j, k) of the cell i *****************
        //first calcuated fluxes on -/+ sides of i-1/2
        fluxim_mjk = reconstructed_flux(vv, im, j, k, s, rh[i], thetaC[j], phi[k]);
        fluxim_pjk = reconstructed_flux(vv, i, j, k, s, rh[i], thetaC[j], phi[k]);

        flux_r[s] = ( (an_imjk[0]*fluxim_mjk - an_imjk[1]*fluxim_pjk)
                     -an_imjk_multi*(Uim_mjk[s]-Uim_pjk[s]))/an_imjk_subtr;
    }
}

/* ------ theta-face fluxes at intefaces (i, j-1/2, k) i ------*/
inline void fluxes_theta(Field ***xx, Field ***uu, Field ***ww, Field ***zz, int i, int j, int k, double flux_theta[])
{
    int jm, kc, s;

    double fluxijm_mk, fluxijm_pk;
    double Uijm_mk[nvar], Uijm_pk[nvar];
    double b_ijmk[2];

    max_speed_theta_face(xx, zz, i, j, k, rfavg[i], thetah[j], phi[k], b_ijmk);
    double b_ijmk_multi = b_ijmk[0]*b_ijmk[1], b_ijmk_subtr = b_ijmk[0] - b_ijmk[1];
    if (b_ijmk[0] >5000.0e3 || fabs(b_ijmk[1]) > 5000.0e3) std::cout<<"b_ijmk "<<i<<" "<<j<<" "<<k<<endl; 

    if (j == 0) {jm = 0; kc = (k+a3/2) % a3;}
    else {jm = j-1; kc = k;}

    for (s = 0; s < 12; s++) {
        Uijm_mk[s] = reconstructed(xx, i, jm, kc, s, rfavg[i], thetah[j], phi[kc]);
        if (j == 0 && (s == 8 || s == 20)) Uijm_mk[s] = -Uijm_mk[s];
        Uijm_pk[s] = reconstructed(xx, i, j, k, s, rfavg[i], thetah[j], phi[k]);

        fluxijm_mk = reconstructed_flux(ww, i, jm, kc, s, rfavg[i], thetah[j], phi[kc]);
        fluxijm_pk = reconstructed_flux(ww, i, j, k, s, rfavg[i], thetah[j], phi[k]);

        flux_theta[s] = ( (b_ijmk[0]*fluxijm_mk - b_ijmk[1]*fluxijm_pk)
                         -b_ijmk_multi*(Uijm_mk[s]-Uijm_pk[s]))/b_ijmk_subtr;
    }

//********* flux for n_q (q=0, 1, 2... sm) equations *****************
    double bn_ijmk[2];
    neu_max_speed_theta_face(xx, uu, i, j, k, bn_ijmk);
    double bn_ijmk_multi = bn_ijmk[0]*bn_ijmk[1], bn_ijmk_subtr = bn_ijmk[0] - bn_ijmk[1];

    for (s = 12; s < 23; s++) {
        Uijm_mk[s] = reconstructed(xx, i, jm, kc, s, rfavg[i], thetah[j], phi[kc]);
        if (j == 0 && s == 20) Uijm_mk[s] = -Uijm_mk[s];
        Uijm_pk[s] = reconstructed(xx, i, j, k, s, rfavg[i], thetah[j], phi[k]);

        fluxijm_mk = reconstructed(ww, i, jm, kc, s, rfavg[i], thetah[j], phi[kc]);
        fluxijm_pk = reconstructed(ww, i, j, k, s, rfavg[i], thetah[j], phi[k]);

        flux_theta[s] = ( (b_ijmk[0]*fluxijm_mk - b_ijmk[1]*fluxijm_pk)
                         -bn_ijmk_multi*(Uijm_mk[s] - Uijm_pk[s]))/bn_ijmk_subtr;
    }
}

/* ------ phi-face fluxes at intefaces (i, j-1/2, k) i ------*/
inline void fluxes_phi(Field ***xx, Field ***uu, Field ***zz, int i, int j, int k, double flux_phi[])
{
    int km, kprime, s;

    double fluxijkm_m, fluxijkm_p;
    double Uijkm_m[nvar], Uijkm_p[nvar];
    double c_ijkm[2];

    if (k == 0) {km = Np; kprime = a3;}
    else {km = k-1; kprime = k;}

    max_speed_phi_face(xx, zz, i, j, k, rfavg[i], theta[j], phih[k], c_ijkm);
    double c_ijkm_multi = c_ijkm[0]*c_ijkm[1], c_ijkm_subtr = c_ijkm[0]-c_ijkm[1];
    if (c_ijkm[0] >5000.0e3 || fabs(c_ijkm[1]) > 5000.0e3) std::cout<<"c_ijmk "<<i<<" "<<j<<" "<<k<<endl; 

    for (s = 0; s < 12; s++) {
        Uijkm_m[s] = reconstructed(xx, i, j, km, s, rfavg[i], theta[j], phih[kprime]);
        Uijkm_p[s] = reconstructed(xx, i, j, k, s, rfavg[i], theta[j], phih[k]);

        fluxijkm_m = reconstructed(zz, i, j, km, s, rfavg[i], theta[j], phih[kprime]);
        fluxijkm_p = reconstructed(zz, i, j, k, s, rfavg[i], theta[j], phih[k]);

        flux_phi[s] = ( (c_ijkm[0]*fluxijkm_m - c_ijkm[1]*fluxijkm_p)
                       -c_ijkm_multi*(Uijkm_m[s] - Uijkm_p[s]))/c_ijkm_subtr;
    }

//********* flux for n_q (q=0, 1, 2... sm) equations *****************
    double cn_ijkm[2];
    neu_max_speed_theta_face(xx, uu, i, j, k, cn_ijkm);
    double cn_ijkm_multi = cn_ijkm[0]*cn_ijkm[1], cn_ijkm_subtr = cn_ijkm[0]-cn_ijkm[1];

    for (s = 12; s < 23; s++) {
        Uijkm_m[s] = reconstructed(xx, i, j, km, s, rfavg[i], theta[j], phih[kprime]);
        Uijkm_p[s] = reconstructed(xx, i, j, k, s, rfavg[i], theta[j], phih[k]);

        fluxijkm_m = reconstructed(zz, i, j, km, s, rfavg[i], theta[j], phih[kprime]);
        fluxijkm_p = reconstructed(zz, i, j, k, s, rfavg[i], theta[j], phih[k]);

        flux_phi[s] = ( (c_ijkm[0]*fluxijkm_m - c_ijkm[1]*fluxijkm_p)
                       -cn_ijkm_multi*(Uijkm_m[s] - Uijkm_p[s]))/cn_ijkm_subtr;
    }
}