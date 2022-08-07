#include "max_mhd_speeds.h"

/*-------------------------------------------------------------------------------------------------------
 *  Cell vaolume averaged electric field components at (rC, thetaC, phi_k)
 *-------------------------------------------------------------------------------------------------------*/
inline void electric_field(Field ***xx, Field ***ww, Field ***uu, int i, int j, int k)
{
    double arcoef, brcoef;
    double atheta, btheta, ctheta;
    double aphi, bphi, cphi;
    double Br, Btheta, Bphi;

    int ip=i+1, jp=j+1, kp=k+1;

    arcoef=0.5*(xx[k][j][ip].fx[23]+xx[k][j][i].fx[23]);
    brcoef=(xx[k][j][ip].fx[23]-xx[k][j][i].fx[23])/dr;
    Br = ww[k][j][i].fx[23] + arcoef +brcoef*(rC[i]-rr[i]);

    atheta=0.5*(xx[k][jp][i].fx[24]+xx[k][j][i].fx[24]);
    btheta=( limited_slope_Btheta_r(xx, i, jp, k)*(thetaC[j]-thetah[j])
            +limited_slope_Btheta_r(xx, i, j, k)*(thetah[jp]-thetaC[j]))/dth;
    ctheta=(xx[k][jp][i].fx[24]-xx[k][j][i].fx[24])/dth;
    Btheta = ww[k][j][i].fx[24] + atheta +btheta*(rC[i]-rfavg[i]) + ctheta*(thetaC[j]-theta[j]);

    if (k < Np) kp = k+1; else kp = 0;

    aphi=0.5*(xx[kp][j][i].fx[25]+xx[k][j][i].fx[25]);
    bphi=( limited_slope_Bphi_r(xx, i, j, kp)*(phi[k]-phih[k])
          +limited_slope_Bphi_r(xx, i, j, k)*(phih[k+1]-phi[k]))/dph;
    cphi=( limited_slope_Bphi_theta(xx, i, j, kp)*(phi[k]-phih[k])
          +limited_slope_Bphi_theta(xx, i, j, k)*(phih[k+1]-phi[k]))/dph;
    Bphi = aphi + bphi*(rC[i]-rfavg[i]) + cphi*(thetaC[j]-theta[j]);

    double ene = e*uu[k][j][i].fx[6];

    // Er at (rC, thetaC, phi_k) (V/m)
    uu[k][j][i].fx[0]= uu[k][j][i].fx[5]*Btheta-uu[k][j][i].fx[4]*Bphi-limited_slope_r(xx,i,j,k,11)/ene;

    // Etheta at (rC, thetaC, phi_k)
    uu[k][j][i].fx[1]= uu[k][j][i].fx[3]*Bphi-uu[k][j][i].fx[5]*Br-limited_slope_theta(xx,i,j,k,11)/(rC[i]*ene);

    // Ephi at (rC, thetaC, phi_k)
    uu[k][j][i].fx[2]= uu[k][j][i].fx[4]*Br-uu[k][j][i].fx[3]*Btheta-limited_slope_phi(xx,i,j,k,11)/(rCsinC[j][i]*ene);
}

/*------------------------------------------------------------------------------
 *   Edge averaged electric field E*_r, E*_theta, E*_phi
 *------------------------------------------------------------------------------*/
inline void Estar(Field ***xx, Field ***uu, int i, int j, int k, Field ***vv)
{
    int im=i-1, jm, km, kc, kmc, kkc, kh;

    double a_imjk[2], b_ijmk[2], c_ijkm[2];
    double a_imj1k[2], a_imjk1[2], b_i1jmk[2], b_ijmk1[2], c_i1jkm[2], c_ij1km[2];

    double a_imjmk[2], a_imjkm[2];
    double b_imjmk[2], b_ijmkm[2];
    double c_imjkm[2], c_ijmkm[2];

    double Br_imjkmp, Br_imjkmm, Br_imjmpk, Br_imjmmk;
    double Btheta_ijmkmp, Btheta_ijmkmm, Btheta_impjmk, Btheta_immjmk;
    double Bphi_ijmpkm, Bphi_ijmmkm, Bphi_impjkm, Bphi_immjkm;

    double Er_ijmpkmp, Er_ijmmkmp, Er_ijmpkmm, Er_ijmmkmm;
    double Etheta_impjkmp, Etheta_immjkmp, Etheta_impjkmm, Etheta_immjkmm;
    double Ephi_impjmpk, Ephi_immjmpk, Ephi_impjmmk, Ephi_immjmmk;

    if (j == 0) {kc=(k+a3/2) % a3; jm = 0;}
    else{kc = k; jm = j-1;}

    if (k == 0) {km = Np; kh = Np+1;}
    else {km = k-1; kh = k;}

    if (j < Nth) {
        max_speed_r_face(xx, i, j, k, rh[i], theta[j], phi[k], a_imjk);
        max_speed_r_face(xx, i, jm, kc, rh[i], theta[jm], phi[kc], a_imj1k);
        max_speed_r_face(xx, i, j, km, rh[i], theta[j], phi[km], a_imjk1);

        max_speed_theta_face(xx, i, j, k, rr[i], thetah[j], phi[k], b_ijmk);
        max_speed_theta_face(xx, im, j, k, rr[im], thetah[j], phi[k], b_i1jmk);
        max_speed_theta_face(xx, i, j, km, rr[i], thetah[j], phi[km], b_ijmk1);

        max_speed_phi_face(xx, i, j, k, rr[i], theta[j], phih[k], c_ijkm);
        max_speed_phi_face(xx, im, j, k, rr[im], theta[j], phih[k], c_i1jkm);
        max_speed_phi_face(xx, i, jm, kc, rr[i], theta[jm], phih[kc], c_ij1km);

        a_imjkm[0]=max(a_imjk[0], a_imjk1[0]);
        a_imjkm[1]=min(a_imjk[1], a_imjk1[1]);

        c_imjkm[0]=max(c_ijkm[0], c_i1jkm[0]);
        c_imjkm[1]=min(c_ijkm[1], c_i1jkm[1]);

        //reconstructed Br, and Bphi
        Br_imjkmp=reconstructed_Br(xx, i, j, k, rh[i], theta[j], phih[k]);
        Br_imjkmm=reconstructed_Br(xx, i, j, km, rh[i], theta[j], phih[kh]);
        Br_imjmpk=reconstructed_Br(xx, i, j, k, rh[i], thetah[j], phi[k]);
        Br_imjmmk=reconstructed_Br(xx, i, jm, kc, rh[i], thetah[j], phi[kc]);

        Btheta_ijmkmp=reconstructed_Btheta(xx, i, j, k, rr[i], thetah[j], phih[k]);
        Btheta_ijmkmm=reconstructed_Btheta(xx, i, j, km, rr[i], thetah[j], phih[kh]);
        Btheta_impjmk=reconstructed_Bphi(xx, i, j, k, rh[i], thetah[j], phi[k]);
        Btheta_immjmk=reconstructed_Bphi(xx, im, j, k, rh[i], thetah[j], phi[k]);

        Bphi_ijmpkm=reconstructed_Bphi(xx, i, j, k, rr[i], thetah[j], phih[k]);
        Bphi_ijmmkm=reconstructed_Bphi(xx, i, jm, kc, rr[i], thetah[j], phih[kc]);
        Bphi_impjkm=reconstructed_Bphi(xx, i, j, k, rh[i], theta[j], phih[k]);
        Bphi_immjkm=reconstructed_Bphi(xx, im, j, k, rh[i], theta[j], phih[k]);

        //reconstructed Er
        if (j == 0) kmc = (km+a3/2) % a3;
        else kmc=km;

        Er_ijmpkmp=reconstructed_efd(uu, i, j, k, 0, rr[i], thetah[j], phih[k]);
        Er_ijmmkmp=reconstructed_efd(uu, i, jm, kc, 0, rr[i], thetah[j], phih[kc]);
        Er_ijmpkmm=reconstructed_efd(uu, i, j, km, 0, rr[i], thetah[j], phih[kh]);
        Er_ijmmkmm=reconstructed_efd(uu, i, jm, kmc, 0, rr[i], thetah[j], phih[kmc]);

        Etheta_impjkmp=reconstructed_efd(uu, i, j, k, 1, rh[i], theta[j], phih[k]);
        Etheta_immjkmp=reconstructed_efd(uu, im, j, k, 1, rh[i], theta[j], phih[k]);
        Etheta_impjkmm=reconstructed_efd(uu, i, j, km, 1, rh[i], theta[j], phih[kh]);
        Etheta_immjkmm=reconstructed_efd(uu, im, j, km, 1, rh[i], theta[j], phih[kh]);

        Ephi_impjmpk=reconstructed_efd(uu, i, j, k, 2, rh[i], thetah[j], phi[k]);
        Ephi_immjmpk=reconstructed_efd(uu, im, j, k, 2, rh[i], thetah[j], phi[k]);
        Ephi_impjmmk=reconstructed_efd(uu, i, jm, kc, 2, rh[i], thetah[j], phi[kc]);
        Ephi_immjmmk=reconstructed_efd(uu, im, jm, kc, 2, rh[i], thetah[j], phi[kc]);

        vv[k][j][i].fx[24]= a_imjkm[1]*c_imjkm[0]/(a_imjkm[0]-a_imjkm[1])*(Bphi_impjkm-Bphi_immjkm)
                           -c_imjkm[1]*c_imjkm[0]/(c_imjkm[0]-c_imjkm[1])*(Br_imjkmp-Br_imjkmm)
                           -( a_imjkm[0]*c_imjkm[1]*Etheta_immjkmp-a_imjkm[1]*c_imjkm[1]*Etheta_impjkmp
                             -a_imjkm[0]*c_imjkm[0]*Etheta_immjkmm+a_imjkm[1]*c_imjkm[0]*Etheta_impjkmm)
                            /((a_imjkm[0]-a_imjkm[1])*(c_imjkm[0]-c_imjkm[1]));
    }
    else {
        kkc=(k+a3/2) % a3; kmc = (km+a3/2) % a3;

        max_speed_r_face(xx, i, Nthm, kkc, rh[i], theta[Nthm], phi[kkc], a_imjk);
        max_speed_r_face(xx, i, jm, k, rh[i], theta[jm], phi[k], a_imj1k);
        //a_imjk1 at j=Nth is not needed

        max_speed_theta_face(xx, i, j, k, rr[i], thetah[j], phi[k], b_ijmk);
        max_speed_theta_face(xx, im, j, k, rr[im], thetah[j], phi[k], b_i1jmk);
        max_speed_theta_face(xx, i, j, km, rr[i], thetah[j], phi[km], b_ijmk1);

        max_speed_phi_face(xx, i, Nthm, kkc, rr[i], theta[Nthm], phih[kkc], c_ijkm);
        //c_i1jkm at j=Nth is not needed;
        max_speed_phi_face(xx, i, jm, k, rr[i], theta[jm], phih[k], c_ij1km);

        //We don't need Br_imjkmp and Br_imjkmm at j=Nth
        Br_imjmpk=reconstructed_Br(xx, i, Nthm, kkc, rh[i], thetah[j], phi[kkc]);
        Br_imjmmk=reconstructed_Br(xx, i, jm, k, rh[i], thetah[j], phi[k]);

        Btheta_ijmkmp=reconstructed_Btheta(xx, i, Nthm, kkc, rr[i], thetah[j], phih[kkc]);
        Btheta_ijmkmm=reconstructed_Btheta(xx, i, Nthm, kmc, rr[i], thetah[j], phih[kmc]);
        Btheta_impjmk=reconstructed_Bphi(xx, i, Nthm, kkc, rh[i], thetah[j], phi[kkc]);
        Btheta_immjmk=reconstructed_Bphi(xx, im, Nthm, kkc, rh[i], thetah[j], phi[kkc]);

        //We don't need Bphi_impjkm and Bphi_immjkm at j=Nth
        Bphi_ijmpkm=reconstructed_Bphi(xx, i, Nthm, kkc, rr[i], thetah[j], phih[kkc]);
        Bphi_ijmmkm=reconstructed_Bphi(xx, i, jm, k, rr[i], thetah[j], phih[k]);

        Er_ijmpkmp=reconstructed_efd(uu, i, Nthm, kkc, 0, rr[i], thetah[j], phih[kkc]);
        Er_ijmmkmp=reconstructed_efd(uu, i, jm, k, 0, rr[i], thetah[j], phih[k]);
        Er_ijmpkmm=reconstructed_efd(uu, i, Nthm, kmc, 0, rr[i], thetah[j], phih[kmc]);
        Er_ijmmkmm=reconstructed_efd(uu, i, jm, km, 0, rr[i], thetah[j], phih[kh]);

        //All E^start_theta are not needed at j=Nth
        Ephi_impjmpk=reconstructed_efd(uu, i, Nthm, kkc, 2, rh[i], thetah[j], phi[kkc]);
        Ephi_immjmpk=reconstructed_efd(uu, im, Nthm, kkc, 2, rh[i], thetah[j], phi[kkc]);
        Ephi_impjmmk=reconstructed_efd(uu, i, jm, k, 2, rh[i], thetah[j], phi[k]);
        Ephi_immjmmk=reconstructed_efd(uu, im, jm, k, 2, rh[i], thetah[j], phi[k]);
    }

    a_imjmk[0]=max(a_imjk[0], a_imj1k[0]);
    a_imjmk[1]=min(a_imjk[1], a_imj1k[1]);

    b_imjmk[0]=max(b_ijmk[0], b_i1jmk[0]);
    b_imjmk[1]=min(b_ijmk[1], b_i1jmk[1]);
    b_ijmkm[0]=max(b_ijmk[0], b_ijmk1[0]);
    b_ijmkm[1]=min(b_ijmk[1], b_ijmk1[1]);

    c_ijmkm[0]=max(c_ijkm[0], c_ij1km[0]);
    c_ijmkm[1]=min(c_ijkm[1], c_ij1km[1]);

    vv[k][j][i].fx[23]= c_ijmkm[1]*c_ijmkm[0]/(c_ijmkm[0]-c_ijmkm[1])*(Btheta_ijmkmp-Btheta_ijmkmm)
                       -b_ijmkm[1]*b_ijmkm[0]/(b_ijmkm[0]-b_ijmkm[1])*(Bphi_ijmpkm-Bphi_ijmmkm)
                       -( b_ijmkm[0]*c_ijmkm[1]*Er_ijmmkmp-b_ijmkm[1]*c_ijmkm[1]*Er_ijmpkmp
                         -b_ijmkm[0]*c_ijmkm[0]*Er_ijmmkmm+b_ijmkm[1]*c_ijmkm[0]*Er_ijmpkmm)
                        /((b_ijmkm[0]-b_ijmkm[1])*(c_ijmkm[0]-c_ijmkm[1]));

    vv[k][j][i].fx[25]= b_imjmk[1]*b_imjmk[0]/(b_imjmk[0]-b_imjmk[1])*(Br_imjmpk-Br_imjmmk)
                       -a_imjmk[1]*a_imjmk[0]/(a_imjmk[0]-a_imjmk[1])*(Btheta_impjmk-Btheta_immjmk)
                       -( a_imjmk[0]*b_imjmk[1]*Ephi_immjmpk-a_imjmk[1]*b_imjmk[1]*Ephi_impjmpk
                         -a_imjmk[0]*b_imjmk[0]*Ephi_immjmmk+a_imjmk[1]*b_imjmk[0]*Ephi_impjmmk)
                        /((a_imjmk[0]-a_imjmk[1])*(b_imjmk[0]-b_imjmk[1]));
}
