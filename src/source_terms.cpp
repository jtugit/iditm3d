#include "ele_cooling_rate.h"
#include "funcdef.h"

/** evaluate source terms excluding collision terms */
void source_terms(Field ***xx, Field ***uu, Field ***ww, Field ***zz, int xs, int i, int j, int k,
    int zk, int yj, int xi, double source[])
{
    int    s;
    const double two3rd=2.0/3.0;
    //double Ti_coll_coef = 0.0, inner_term = 0.0;
    //double Tn_coll_coef = 0.0, inner_term_n = 0.0;

    double ns_jC[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0}, nn_jC;
    double rhoi_jC = 0.0, rhoi_j = 0.0, ne_jC=0.0;
    double rhon_jC=0.0, rhon_j=0.0;
    const double ni00=ni_0/1.0e6, nn00=nn_0/1.0e6;

    for (s = 0; s < sl; s++) {
        //reconstructed ion mass density at (rfavg, thetaC, phi)
        ns_jC[s] = reconstructed(xx, i, j, k, s, rfavg[i], thetaC[j], phi[k]);
        rhoi_jC += ms[s]*ns_jC[s];
        ne_jC += ns_jC[s];

        //reconstructed ion mass density at (rfavg, theta, phi)
        rhoi_j += ms[s]*reconstructed(xx, i, j, k, s, rfavg[i], theta[j], phi[k]);

        //reconstructed rho_{n} at (rfavg, thetaC, phi_k)
        nn_jC=reconstructed(xx, i, j, k, 12+s, rfavg[i], thetaC[j], phi[k]);
        rhon_jC += ms[s]*nn_jC;

        //reconstructed rho_{n} at (rfavg, theta, phi_k)
        rhon_j += ms[s]*reconstructed(xx, i, j, k, 12+s, rfavg[i], theta[j], phi[k]);      

        //source terms for ion and neutral continuity equations
        source[s] = Ps[zk][yj][xi][s]/ni00; // - Ls[s]*ns_jC[s];
        source[12+s] = Ps[zk][yj][xi][7+s]/nn00; // - Ls[7+s]*nn_jC;
    }

    //reconstructed quantities at (rfavg, thetaC, phi_k)
    double uir_jC, uitheta_jC, uiphi_jC;
    uir_jC = reconstructed(uu, i, j, k, 7, rfavg[i], thetaC[j], phi[k]);
    uitheta_jC = reconstructed(uu, i, j, k, 8, rfavg[i], thetaC[j], phi[k]);
    uiphi_jC = reconstructed(uu, i, j, k, 9, rfavg[i], thetaC[j], phi[k]);

    double Br_jC, Btheta_jC, Bphi_jC;
    Br_jC = ww[k][j][i].fx[23]+reconstructed_Br(xx, i, j, k, rfavg[i], thetaC[j], phi[k]);
    Btheta_jC = ww[k][j][i].fx[24]+reconstructed_Btheta(xx, i, j, k, rfavg[i], thetaC[j], phi[k]);
    Bphi_jC = ww[k][j][i].fx[25]+reconstructed_Bphi(xx, i, j, k, rfavg[i], thetaC[j], phi[k]);

    double Pi_jC, Pe_jC;
    Pi_jC = reconstructed(xx, i, j, k, 10, rfavg[i], thetaC[j], phi[k]);
    Pe_jC = reconstructed(xx, i, j, k, 11, rfavg[i], thetaC[j], phi[k]);

    //reconstructed quantities at (rfavg, theta, phi_k)
    double uitheta_j, uiphi_j;
    uitheta_j = reconstructed(uu, i, j, k, 8, rfavg[i], theta[j], phi[k]);
    uiphi_j = reconstructed(uu, i, j, k, 9, rfavg[i], theta[j], phi[k]);

    double Br_j, Btheta_j, Bphi_j;
    Br_j = zz[k][j][i].fx[23]+reconstructed_Br(xx, i, j, k, rfavg[i], theta[j], phi[k]);
    Btheta_j = zz[k][j][i].fx[24]+reconstructed_Btheta(xx, i, j, k, rfavg[i], theta[j], phi[k]);
    Bphi_j = reconstructed_Bphi(xx, i, j, k, rfavg[i], theta[j], phi[k]);

    double Pi_j, Pe_j;
    Pi_j = reconstructed(xx, i, j, k, 10, rfavg[i], theta[j], phi[k]);
    Pe_j = reconstructed(xx, i, j, k, 11, rfavg[i], theta[j], phi[k]);

    //source terms for ion velocity
    source[7] = (rhoi_jC*(uitheta_jC*uitheta_jC + uiphi_jC*uiphi_jC)+2.0*(Pi_jC+Pe_jC)+Br_jC*Br_jC/ni_0)/rfavg[i]
               -rhoi_jC*gr[i-xs];
    source[8] = ( (rhoi_j*uiphi_j*uiphi_j+Pi_j+Pe_j+0.5/ni_0*(Br_j*Br_j+Btheta_j*Btheta_j-Bphi_j*Bphi_j))*cotth[j]
                 +(Btheta_jC*Br_jC/ni_0-rhoi_jC*uitheta_jC*uir_jC))/rfavg[i];
    source[9] = ( Bphi_jC*Br_jC/ni_0 - rhoi_jC*uiphi_jC*uir_jC
                 +(Bphi_j*Btheta_j/ni_0 - rhoi_j*uiphi_j*uitheta_j)*cotth[j])/rfavg[i];

    // for source terms of ion and electron energy equations
    double delta_uir, delta_uitheta, delta_uiphi;

    delta_uir = limited_slope_r(uu, i, j, k, 7);
    delta_uitheta = limited_slope_r(uu, i, j, k, 8);
    delta_uiphi = limited_slope_r(uu, i, j, k, 9);

    double unr_jC = reconstructed(uu, i, j, k, 19, rfavg[i], thetaC[j], phi[k]);
    double untheta_jC = reconstructed(uu, i, j, k, 20, rfavg[i], thetaC[j], phi[k]);
    double unphi_jC = reconstructed(uu, i, j, k, 21, rfavg[i], thetaC[j], phi[k]);

    source[10] =-two3rd*( Pi_jC*(delta_uir+(2.0*uir_jC+delta_uitheta+delta_uiphi/sinth[j])/rfavg[i])
                         +cotth[j]/rfavg[i]*Pi_j*uitheta_j);

    double delta_uer = limited_slope_r(uu, i, j, k, 3);
    double delta_uetheta = limited_slope_r(uu, i, j, k, 4);
    double delta_uephi = limited_slope_r(uu, i, j, k, 5);
    double uer_jC = reconstructed(uu, i, j, k, 3, rfavg[i], thetaC[j], phi[k]);
    double uetheta_j = reconstructed(uu, i, j, k, 4, rfavg[i], theta[j], phi[k]);

    //Qee must be in unit of Joule cm^{-3} s^{-1}
    source[11] =two3rd*(Qee[zk][yj][xi]-Pe_jC*(delta_uer+(2.0*uer_jC+delta_uetheta+delta_uephi/sinth[j])/rfavg[i])
                           -(Pe_j*uetheta_j)*cotth[j]/rfavg[i]);

    double untheta_j = reconstructed(uu, i, j, k, 20, rfavg[i], theta[j], phi[k]);
    double unphi_j = reconstructed(uu, i, j, k, 21, rfavg[i], theta[j], phi[k]);
    double Pn_jC = reconstructed(xx, i, j, k, 22, rfavg[i], thetaC[j], phi[k]);
    double Pn_j = reconstructed(xx, i, j, k, 22, rfavg[i], theta[j], phi[k]);

    source[19] = (rhon_jC*(untheta_jC*untheta_jC+unphi_jC*unphi_jC) + 2.0*Pn_jC)/rfavg[i]
                - rhon_jC*( 2.0*(rotat_t[zk][yj]*unphi_jC-rotat_p[zk][yj]*untheta_jC)
                           -cenf_r[zk][yj][xi] + gr[xi]);
    source[20] = ((rhon_j*unphi_j*unphi_j+Pn_j)*cotth[j] - rhon_jC*untheta_jC*unr_jC)/rfavg[i]
                -rhon_jC*(cenf_t[zk][yj][xi]+2.0*(rotat_p[zk][yj]*unr_jC-rotat_r[zk][yj]*unphi_jC));
    source[21] = -rhon_jC*(cenf_p[zk][yj][xi]+2.0*(rotat_r[zk][yj]*untheta_jC-rotat_t[zk][yj]*unr_jC))
                 -(rhon_jC*unphi_jC*unr_jC + rhon_j*unphi_j*untheta_j*cotth[j])/rfavg[i];

    double delta_unr = limited_slope_r(uu, i, j, k, 19);
    double delta_untheta = limited_slope_r(uu, i, j, k, 20);
    double delta_unphi = limited_slope_r(uu, i, j, k, 21);
    double Cn = neu_cooling_rate(xx, uu, i, j, k);

    source[22] = two3rd*( (Qeuv[zk][yj][xi]-Cn)
                         -Pn_jC*(delta_unr+(2.0*unr_jC+delta_untheta+delta_unphi/sinth[j])/rfavg[i])
                         -(Pn_j*untheta_j)*cotth[j]/rfavg[i]);
}