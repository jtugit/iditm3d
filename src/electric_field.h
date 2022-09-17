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

    if (i==Nr) {
        arcoef=0.5*(xx[k][j][i].fx[23]+xx[k][j][i-1].fx[23]);
        brcoef=(xx[k][j][i].fx[23]-xx[k][j][i-1].fx[23])/dr;
    }
    else {
        arcoef=0.5*(xx[k][j][ip].fx[23]+xx[k][j][i].fx[23]);
        brcoef=(xx[k][j][ip].fx[23]-xx[k][j][i].fx[23])/dr;
    }
    Br = ww[k][j][i].fx[23] + arcoef +brcoef*(rC[i]-rr[i]);

    atheta=0.5*(xx[k][jp][i].fx[24]+xx[k][j][i].fx[24]);
    btheta=( limited_slope_Btheta_r(xx, i, jp, k)*(thetaC[j]-thetah[j])
            +limited_slope_Btheta_r(xx, i, j, k)*(thetah[jp]-thetaC[j]))/dth;
    ctheta=(xx[k][jp][i].fx[24]-xx[k][j][i].fx[24])/dth;
    Btheta = ww[k][j][i].fx[24] + atheta +btheta*(rC[i]-rfavg[i]) + ctheta*(thetaC[j]-theta[j]);

    aphi=0.5*(xx[kp][j][i].fx[25]+xx[k][j][i].fx[25]);
    bphi=( limited_slope_Bphi_r(xx, i, j, kp)*(phi[k]-phih[k])
          +limited_slope_Bphi_r(xx, i, j, k)*(phih[k+1]-phi[k]))/dph;
    cphi=( limited_slope_Bphi_theta(xx, i, j, kp)*(phi[k]-phih[k])
          +limited_slope_Bphi_theta(xx, i, j, k)*(phih[k+1]-phi[k]))/dph;
    Bphi = aphi + bphi*(rC[i]-rfavg[i]) + cphi*(thetaC[j]-theta[j]);

    double ene = e*uu[k][j][i].fx[17];

    // Er at (rC, thetaC, phi_k) (V/m)
    uu[k][j][i].fx[0]= uu[k][j][i].fx[9]*Btheta-uu[k][j][i].fx[8]*Bphi-limited_slope_r(xx,i,j,k,11)/ene;

    // Etheta at (rC, thetaC, phi_k)
    uu[k][j][i].fx[1]= uu[k][j][i].fx[7]*Bphi-uu[k][j][i].fx[9]*Br-limited_slope_theta(xx,i,j,k,11)/(rC[i]*ene);

    // Ephi at (rC, thetaC, phi_k)
    uu[k][j][i].fx[2]= uu[k][j][i].fx[8]*Br-uu[k][j][i].fx[7]*Btheta-limited_slope_phi(xx,i,j,k,11)/(rCsinC[j][i]*ene);
}

/*------------------------------------------------------------------------------
 *   Edge averaged electric field E*_r, E*_theta, E*_phi
 *------------------------------------------------------------------------------*/
inline void Estar(Field ***xx, Field ***uu, int i, int j, int k, Field ***vv)
{
    int im=i-1, jm, km=k-1, kc, kh, kcm;

    double a_imjmk[2], a_imjkm[2];
    double b_imjmk[2], b_ijmkm[2];
    double c_imjkm[2], c_ijmkm[2];

    double Br_imjkmp, Br_imjkmm, Br_imjmpk, Br_imjmmk;
    double Btheta_ijmkmp, Btheta_ijmkmm, Btheta_impjmk, Btheta_immjmk;
    double Bphi_ijmpkm, Bphi_ijmmkm, Bphi_impjkm, Bphi_immjkm;

    double Er_ijmpkmp, Er_ijmmkmp, Er_ijmpkmm, Er_ijmmkmm;
    double Etheta_impjkmp, Etheta_immjkmp, Etheta_impjkmm, Etheta_immjkmm;
    double Ephi_impjmpk, Ephi_immjmpk, Ephi_impjmmk, Ephi_immjmmk;

    if (k == 0) {km=Np; kh=km+1;} 
    else {km=k-1; kh=k;}

    if (j < Nth) {
        if (j == 1) {
            kc=(k+a3/2) % a3; jm = j; kcm = (k-1+a3/2) % a3;
        }
        else{
            kc = k; jm = j-1; kcm=k-1;
        }

        a_imjmk[0]=max(uu[k][j][i].fx[3], uu[kc][jm][i].fx[3]);
        a_imjmk[1]=min(uu[k][j][i].fx[4], uu[kc][jm][i].fx[4]);
        a_imjkm[0]=max(uu[k][j][i].fx[3], uu[km][j][i].fx[3]);
        a_imjkm[1]=min(uu[k][j][i].fx[4], uu[km][j][i].fx[4]);

        b_imjmk[0]=max(uu[k][j][i].fx[5], uu[k][j][im].fx[5]);
        b_imjmk[1]=min(uu[k][j][i].fx[6], uu[k][j][im].fx[6]);
        b_ijmkm[0]=max(uu[k][j][i].fx[5], uu[km][j][i].fx[5]);
        b_ijmkm[1]=min(uu[k][j][i].fx[6], uu[km][j][i].fx[6]);

        c_imjkm[0]=max(uu[k][j][i].fx[15], uu[k][j][im].fx[15]);
        c_imjkm[1]=min(uu[k][j][i].fx[16], uu[k][j][im].fx[16]);
        c_ijmkm[0]=max(uu[k][j][i].fx[15], uu[kc][jm][i].fx[15]);
        c_ijmkm[1]=min(uu[k][j][i].fx[16], uu[kc][jm][i].fx[16]);

        //reconstructed Br, and Bphi
        Br_imjkmp=reconstructed_Br(xx, i, j, k, rh[i], theta[j], phih[k]);
        Br_imjkmm=reconstructed_Br(xx, i, j, km, rh[i], theta[j], phih[kh]);
        Br_imjmpk=reconstructed_Br(xx, i, j, k, rh[i], thetah[j], phi[k]);
        Br_imjmmk=reconstructed_Br(xx, i, jm, kc, rh[i], thetah[j], phi[kc]);

        Btheta_ijmkmp=reconstructed_Btheta(xx, i, j, k, rr[i], thetah[j], phih[k]);
        Btheta_ijmkmm=reconstructed_Btheta(xx, i, j, km, rr[i], thetah[j], phih[kh]);
        Btheta_impjmk=reconstructed_Btheta(xx, i, j, k, rh[i], thetah[j], phi[k]);
        Btheta_immjmk=reconstructed_Btheta(xx, im, j, k, rh[i], thetah[j], phi[k]);

        Bphi_ijmpkm=reconstructed_Bphi(xx, i, j, k, rr[i], thetah[j], phih[k]);
        Bphi_ijmmkm=reconstructed_Bphi(xx, i, jm, kc, rr[i], thetah[j], phih[kc]);
        Bphi_impjkm=reconstructed_Bphi(xx, i, j, k, rh[i], theta[j], phih[k]);
        Bphi_immjkm=reconstructed_Bphi(xx, im, j, k, rh[i], theta[j], phih[k]);

        Er_ijmpkmp=reconstructed_efd(uu, i, j, k, 0, rr[i], thetah[j], phih[k]);
        Er_ijmmkmp=reconstructed_efd(uu, i, jm, kc, 0, rr[i], thetah[j], phih[kc]);
        Er_ijmpkmm=reconstructed_efd(uu, i, j, km, 0, rr[i], thetah[j], phih[kh]);
        Er_ijmmkmm=reconstructed_efd(uu, i, jm, kcm, 0, rr[i], thetah[j], phih[kcm+1]);

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
        jm=j-1; 
        int kcp=(k+a3/2) % a3;
        kcm=(k-1+a3/2) % a3;

        a_imjmk[0]=max(uu[kcp][jm][i].fx[3], uu[k][jm][i].fx[3]);
        a_imjmk[1]=min(uu[kcp][jm][i].fx[4], uu[k][jm][i].fx[4]);

        b_imjmk[0]=max(uu[kcp][jm][i].fx[5], uu[kcp][jm][im].fx[5]);
        b_imjmk[1]=min(uu[kcp][jm][i].fx[6], uu[kcp][jm][im].fx[6]);
        b_ijmkm[0]=max(uu[kcp][jm][i].fx[5], uu[kcm][jm][i].fx[5]);
        b_ijmkm[1]=min(uu[kcp][jm][i].fx[6], uu[kcm][jm][i].fx[6]);

        c_ijmkm[0]=max(uu[kcp][jm][i].fx[15], uu[k][jm][i].fx[15]);
        c_ijmkm[1]=min(uu[kcp][jm][i].fx[16], uu[k][jm][i].fx[16]);

        //reconstructed Br, and Bphi
        Br_imjmpk=reconstructed_Br(xx, i, jm, kcp, rh[i], thetah[j], phi[kcp]);
        Br_imjmmk=reconstructed_Br(xx, i, jm, k, rh[i], thetah[j], phi[k]);

        Btheta_ijmkmp=reconstructed_Btheta(xx, i, j, k, rr[i], thetah[j], phih[k]);
        Btheta_ijmkmm=reconstructed_Btheta(xx, i, j, km, rr[i], thetah[j], phih[kh]);
        Btheta_impjmk=reconstructed_Btheta(xx, i, j, k, rh[i], thetah[j], phi[k]);
        Btheta_immjmk=reconstructed_Btheta(xx, im, j, k, rh[i], thetah[j], phi[k]);

        Bphi_ijmpkm=reconstructed_Bphi(xx, i, jm, kcp, rr[i], thetah[j], phih[kcp]);
        Bphi_ijmmkm=reconstructed_Bphi(xx, i, jm, k, rr[i], thetah[j], phih[k]);

        Er_ijmpkmp=reconstructed_efd(uu, i, jm, kcp, 0, rr[i], thetah[j], phih[kcp]);
        Er_ijmmkmp=reconstructed_efd(uu, i, jm, k, 0, rr[i], thetah[j], phih[k]);
        Er_ijmpkmm=reconstructed_efd(uu, i, jm, kcm, 0, rr[i], thetah[j], phih[kcp]);
        Er_ijmmkmm=reconstructed_efd(uu, i, jm, km, 0, rr[i], thetah[j], phih[kh]);

        Ephi_impjmpk=reconstructed_efd(uu, i, jm, kcp, 2, rh[i], thetah[j], phi[kcp]);
        Ephi_immjmpk=reconstructed_efd(uu, im, jm, kcp, 2, rh[i], thetah[j], phi[kcp]);
        Ephi_impjmmk=reconstructed_efd(uu, i, jm, k, 2, rh[i], thetah[j], phi[k]);
        Ephi_immjmmk=reconstructed_efd(uu, im, jm, k, 2, rh[i], thetah[j], phi[k]);
    }

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
