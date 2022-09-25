#include "reconstruction.h"

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
inline void Estar(Field ***uu, int i, int j, int k, Field ***vv)
{
    int im=i-1, jm, km=k-1, kc, kh, kcm;

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
            kc = k; jm = j-1;
            if (k == 0) kcm=Np; else kcm=k-1;
        }

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

        vv[k][j][i].fx[24]= 0.25*(Etheta_immjkmp+Etheta_impjkmp+Etheta_immjkmm+Etheta_impjkmm);
    }
    else {
        jm=j-1; 
        int kcp=(k+a3/2) % a3;
        kcm=(k-1+a3/2) % a3;

        Er_ijmpkmp=reconstructed_efd(uu, i, jm, kcp, 0, rr[i], thetah[j], phih[kcp]);
        Er_ijmmkmp=reconstructed_efd(uu, i, jm, k, 0, rr[i], thetah[j], phih[k]);
        Er_ijmpkmm=reconstructed_efd(uu, i, jm, kcm, 0, rr[i], thetah[j], phih[kcp]);
        Er_ijmmkmm=reconstructed_efd(uu, i, jm, km, 0, rr[i], thetah[j], phih[kh]);

        Ephi_impjmpk=reconstructed_efd(uu, i, jm, kcp, 2, rh[i], thetah[j], phi[kcp]);
        Ephi_immjmpk=reconstructed_efd(uu, im, jm, kcp, 2, rh[i], thetah[j], phi[kcp]);
        Ephi_impjmmk=reconstructed_efd(uu, i, jm, k, 2, rh[i], thetah[j], phi[k]);
        Ephi_immjmmk=reconstructed_efd(uu, im, jm, k, 2, rh[i], thetah[j], phi[k]);
    }

    vv[k][j][i].fx[23]= 0.25*(Er_ijmmkmp+Er_ijmpkmp+Er_ijmmkmm+Er_ijmpkmm);

    vv[k][j][i].fx[25]= 0.25*(Ephi_immjmpk+Ephi_impjmpk+Ephi_immjmmk+Ephi_impjmmk);
}
