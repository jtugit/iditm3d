#include "param.h"
#include "operators.h"

inline void E_gradPe(Field ***xx, Field ***uu, int i, int j, int k, int yj, int zk, double Ec_gradPe[]) 
{
    double dTe_dr, dne_dr, dTe_dth, dne_dth, dTe_dph, dne_dph;
    double ne, Te, dpe_dr, dpe_dth, dpe_dph, inv_ene;

    dTe_dr=difference_r(xx, i, j, k, 19);       //dTe/dr
    dTe_dth=difference_theta(xx, i, j, k, 19);  //dTe/dtheta
    dTe_dph=difference_phi(xx, i, j, k, 19);    //dTe/dphi

    dne_dr=difference_r(uu, i, j, k, 6);
    dne_dth=difference_theta(uu, i, j, k, 6);
    dne_dph=difference_phi(uu, i, j, k, 6);

    Te=xx[k][j][i].fx[19];
    ne=uu[k][j][i].fx[6];

    dpe_dr =ne*dTe_dr+Te*dne_dr;
    dpe_dth=(ne*dTe_dth+Te*dne_dth);
    dpe_dph=(ne*dTe_dph+Te*dne_dph);

    //negative of electric field part due to electron pressure gradient in Cartesian coordinates
    inv_ene=1.0/(e*ne);
    Ec_gradPe[0] =(K11[zk][yj]*dpe_dr+K12[zk][yj]*dpe_dth+K13[zk]*dpe_dph)*inv_ene;
    Ec_gradPe[1] =(K21[zk][yj]*dpe_dr+K22[zk][yj]*dpe_dth+K23[zk]*dpe_dph)*inv_ene;
    Ec_gradPe[2] =(K31[yj]*dpe_dr+K32[yj]*dpe_dth)*inv_ene;
}

/*-------------------------------------------------------------------------------------------------------
 *  Cartesian components of electric field: -(u_e x B) part
 *-------------------------------------------------------------------------------------------------------*/
inline void electric_field_vxB(Field ***xx, Field ***uu, double ue[], int i, int j, int k, int xi,
    int yj, int zk, double Ec_VxB[])
{
    double Bx, By, Bz, uex, uey, uez;

    //Cartesian Bx, By, and Bz in terms of special sperhical components of B
    Bx=( Jiv11[zk][yj]*xx[k][j][i].fx[31]+Jiv12[zk][yj][xi]*xx[k][j][i].fx[32]
         +Jiv13[zk][yj][xi]*xx[k][j][i].fx[33])/r2sintheta[yj][xi] + uu[k][j][i].fx[3];
    By=( Jiv21[zk][yj]*xx[k][j][i].fx[31]+Jiv22[zk][yj][xi]*xx[k][j][i].fx[32]
         +Jiv23[zk][yj][xi]*xx[k][j][i].fx[33])/r2sintheta[yj][xi] + uu[k][j][i].fx[4];
    Bz=(Jiv31[yj]*xx[k][j][i].fx[31]+Jiv32[yj][xi]*xx[k][j][i].fx[32])/r2sintheta[yj][xi]+uu[k][j][i].fx[5];

    //uex, uey, uez in terms of uer, and ue_theta, and ue_phi
    uex = (K11[zk][yj]*ue[0]+K12[zk][yj]*ue[1]+K13[zk]*ue[2]);
    uey = (K21[zk][yj]*ue[0]+K22[zk][yj]*ue[1]+K23[zk]*ue[2]);
    uez = (K31[yj]*ue[0]+K32[yj]*ue[1]);

    //Cartesian components of the electric field
    Ec_VxB[0]=By*uez - Bz*uey; Ec_VxB[1]=Bz*uex - Bx*uez; Ec_VxB[2]=Bx*uey - By*uex;
}

inline void currents(Field ***localuu, Field ***uu, int i, int j, int k, int yj, int xi)
{
    double dBphi_dtheta = difference_theta(localuu, i, j, k, 26);
    double dBtheta_dphi = difference_phi(localuu, i, j, k, 25);
    double dBr_dphi = difference_phi(localuu, i, j, k, 24);
    double dBphi_dr = difference_r(localuu, i, j, k, 26);
    double dBtheta_dr = difference_r(localuu, i, j, k, 25);
    double dBr_dtheta = difference_theta(localuu, i, j, k, 24);

    uu[k][j][i].fx[16]=cot_div_r[yj][xi]*localuu[k][j][i].fx[26]+dBphi_dtheta/rr[i]-dBtheta_dphi/rsin[yj][xi];
    uu[k][j][i].fx[17]=dBr_dphi/rsin[yj][xi] - dBphi_dr - localuu[k][j][i].fx[26]/rr[i];
    uu[k][j][i].fx[18]=localuu[k][j][i].fx[25]/rr[i]+dBtheta_dr - dBr_dtheta/rr[i];
}
