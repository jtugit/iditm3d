#include "param.h"
#include "operators.h"

inline void E_gradPe(Field ***xx, Field ***uu, int i, int j, int k, int xi, int yj, int zk, 
    int xm, int ym, double Ec_gradPe[]) 
{
    double dTe_dr, dne_dr, dTe_dth, dne_dth, dTe_dph, dne_dph;
    double ne, Te, dpe_dr, dpe_dth, dpe_dph, inv_ene;
    uint64_t kj, ji, kji;

    kj = (uint64_t)(zk*ym+yj); ji = (uint64_t)(yj*xm+xi); kji=(uint64_t)(zk*ym*xm+yj*xm+xi);

    dTe_dr=difference_r(xx, i, j, k, 19);
    dTe_dth=difference_theta(xx, i, j, k, 19);
    dTe_dph=difference_phi(xx, i, j, k, 19);

    dne_dr=difference_r(uu, i, j, k, 6);
    dne_dth=difference_theta(uu, i, j, k, 6);
    dne_dph=difference_phi(uu, i, j, k, 6);

    Te=xx[k][j][i].fx[19];
    ne=uu[k][j][i].fx[6];

    dpe_dr =ne*dTe_dr+Te*dne_dr;
    dpe_dth=(ne*dTe_dth+Te*dne_dth);
    dpe_dph=(ne*dTe_dph+Te*dne_dph);

    //negative of electric field part due to electron pressure gradient in Cartesian coordinates
    inv_ene=1.0/(e*uu[k][j][i].fx[6]);
    Ec_gradPe[0] =(Jmat.J11[kj]*dpe_dr+Jmat.J21[kji]*dpe_dth+Jmat.J31[kji]*dpe_dph)*inv_ene;
    Ec_gradPe[1] =(Jmat.J12[kj]*dpe_dr+Jmat.J22[kji]*dpe_dth+Jmat.J32[kji]*dpe_dph)*inv_ene;
    Ec_gradPe[2] =(Jmat.J13[(uint64_t)yj]*dpe_dr+Jmat.J23[ji] *dpe_dth)*inv_ene;
}

/*-------------------------------------------------------------------------------------------------------
 *  Cartesian components of electric field: -(u_e x B) part
 *-------------------------------------------------------------------------------------------------------*/
inline void electric_field_vxB(Field ***xx, Field ***uu, double ue[], int i, int j, int k, int xi,
    int yj, int zk, int xm, int ym, double Ec_VxB[])
{
    uint64_t kj = (uint64_t)(zk*ym+yj), kji=(uint64_t)(zk*ym*xm+yj*xm+xi), ji=(uint64_t)(yj*xm+xi);
    double Bx, By, Bz, uex, uey, uez;

    //Cartesian Bx, By, and Bz in terms of special sperhical components of B
    Bx = ( Jinv.Jiv11[kj]*xx[k][j][i].fx[31]+Jinv.Jiv12[kji]*xx[k][j][i].fx[32]
          +Jinv.Jiv13[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi] + uu[k][j][i].fx[3];
    By = ( Jinv.Jiv21[kj]*xx[k][j][i].fx[31]+Jinv.Jiv22[kji]*xx[k][j][i].fx[32]
          +Jinv.Jiv23[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi] + uu[k][j][i].fx[4];
    Bz = (Jinv.Jiv31[(uint64_t)yj]*xx[k][j][i].fx[31]+Jinv.Jiv32[ji]*xx[k][j][i].fx[32])/r2sintheta[yj][xi]
        + uu[k][j][i].fx[5];

    //uex, uey, uez in terms of uer, and ue_theta, and ue_phi
    uex = (Kmat.K11[kj]*ue[0]+Kmat.K12[kj]*ue[1]+Kmat.K13[(uint64_t)zk]*ue[2]);
    uey = (Kmat.K21[kj]*ue[0]+Kmat.K22[kj]*ue[1]+Kmat.K23[(uint64_t)zk]*ue[2]);
    uez = (Kmat.K31[(uint64_t)yj]*ue[0]+Kmat.K32[(uint64_t)yj]*ue[1]);

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

    double rsindph=rsin[yj][xi]*dph;

    uu[k][j][i].fx[16]=cot_div_r[yj][xi]*localuu[k][j][i].fx[26]+dBphi_dtheta/rr[i]-dBtheta_dphi/rsindph;
    uu[k][j][i].fx[17]=dBr_dphi/rsindph - dBphi_dr - localuu[k][j][i].fx[26]/rr[i];
    uu[k][j][i].fx[18]=localuu[k][j][i].fx[25]/rr[i]+dBtheta_dr - dBr_dtheta/rr[i];
}
