#include "param.h"
#include "operators.h"

inline vector3D E_gradPe(Field ***xx, Field ***uu, int i, int j, int k, int xi, int yj, int zk, int xm, int ym, int zm) 
{
    double dTe_dr, dne_dr, dTe_dth, dne_dth, dTe_dph, dne_dph;
    double ne, Te, dpe_dr, dpe_dth, dpe_dph, inv_ene;
    double Epart_x, Epart_y, Epart_z, kb_ene;
    int    kj, ji, kji;
    vector3D efd_gradPe;

    kj = zk*ym+yj;

    ji = yj*xm+xi; kji=zk*ym*xm+yj*xm+xi;

    dTe_dr=difference_r(xx, i, j, k, 19);
    dTe_dth=difference_theta(xx, i, j, k, 19);
    dTe_dph=difference_phi(xx, i, j, k, 19);

    dne_dr=difference_r(uu, i, j, k, 6);
    dne_dth=difference_theta(uu, i, j, k, 6);
    dne_dph=difference_phi(uu, i, j, k, 6);

    Te=xx[k][j][i].fx[19];
    ne=uu[k][j][i].fx[6];

    dpe_dr =ne*dTe_dr+Te*dne_dr;
    dpe_dth=(ne*dTe_dth+Te*dne_dth)/rr[i];
    dpe_dph=(ne*dTe_dph+Te*dne_dph)/rsin[yj][xi];

    //part of electric field due to electron pressure gradient in Cartesian coordinates
    Epart_x =-(Jmat.J11[kj]*dpe_dr+Jmat.J21[kji]*dpe_dth+Jmat.J31[kji]*dpe_dph)/inv_ene;
    Epart_y =-(Jmat.J12[kj]*dpe_dr+Jmat.J22[kji]*dpe_dth+Jmat.J32[kji]*dpe_dph)/inv_ene;
    Epart_z =-(Jmat.J13[yj]*dpe_dr+Jmat.J23[ji] *dpe_dth)/inv_ene;

    //convert this part of electric field into special spherical components
    efd_gradPe.r=Jinv.Jiv11[kj] *Epart_x+Jinv.Jiv21[kj] *Epart_y+Jinv.Jiv31[yj]*Epart_z;
    efd_gradPe.t=Jinv.Jiv12[kji]*Epart_x+Jinv.Jiv22[kji]*Epart_y+Jinv.Jiv32[ji]*Epart_z;
    efd_gradPe.p=Jinv.Jiv13[kji]*Epart_x+Jinv.Jiv23[kji]*Epart_y;

    return efd_gradPe;
}
/*-------------------------------------------------------------------------------------------------------
 *  special spherical components of electric field: -(u_e x B) part
 *-------------------------------------------------------------------------------------------------------*/
inline void electric_field_vxB(Field ***xx, Field ***uu, int i, int j, int k, int xi, int yj, int zk, int xm, int ym, int zm)
{
    int kj = zk*ym+yj, kji=zk*ym*xm+yj*xm+xi, ji=yj*xm+xi;

    double Bx, By, Bz, uex, uey, uez;

    //Cartesian Bx, By, and Bz in terms of special sperhical components of B
    Bx = ( Jinv.Jiv11[kj]*xx[k][j][i].fx[31]+Jinv.Jiv12[kji]*xx[k][j][i].fx[32]
          +Jinv.Jiv13[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi] + uu[k][j][i].fx[25];
    By = ( Jinv.Jiv21[kj]*xx[k][j][i].fx[31]+Jinv.Jiv22[kji]*xx[k][j][i].fx[32]
          +Jinv.Jiv23[kji]*xx[k][j][i].fx[33])/r2sintheta[yj][xi] + uu[k][j][i].fx[26];
    Bz = (Jinv.Jiv31[yj]*xx[k][j][i].fx[31]+Jinv.Jiv32[ji]*xx[k][j][i].fx[32])/r2sintheta[yj][xi]
        + uu[k][j][i].fx[27];

    //uex, uey, uez in terms of uer, and ue_theta, and ue_phi
    uex = (Kmat.K11[kj]*uu[k][j][i].fx[7]+Kmat.K12[kj]*uu[k][j][i].fx[8]+Kmat.K13[zk]*uu[k][j][i].fx[9]);
    uey = (Kmat.K21[kj]*uu[k][j][i].fx[7]+Kmat.K22[kj]*uu[k][j][i].fx[8]+Kmat.K23[zk]*uu[k][j][i].fx[9]);
    uez = (Kmat.K31[yj]*uu[k][j][i].fx[7]+Kmat.K32[yj]*uu[k][j][i].fx[8]);

    //Cartesian components of the electric field
    double Epartx=By*uez - Bz*uey, Eparty=Bz*uex - Bx*uez, Epartz=Bx*uey - By*uex;

    //convert Cartesian components EFD to special spherical components
    uu[k][j][i].fx[3] = Jinv.Jiv11[kj] *Epartx+Jinv.Jiv21[kj] *Eparty+Jinv.Jiv31[yj]*Epartz;
    uu[k][j][i].fx[4] = Jinv.Jiv12[kji]*Epartx+Jinv.Jiv22[kji]*Eparty+Jinv.Jiv32[ji]*Epartz;
    uu[k][j][i].fx[5] = Jinv.Jiv13[kji]*Epartx+Jinv.Jiv23[kji]*Eparty;
}

inline void currents(Field ***xx, Field ***uu, int i, int j, int k, int yj, int xi)
{
    double dBphi_dtheta = difference_theta(xx, i, j, k, 33);
    double dBtheta_dphi = difference_phi(xx, i, j, k, 32);
    double dBr_dphi = difference_phi(xx, i, j, k, 31);
    double dBphi_dr = difference_r(xx, i, j, k, 33);
    double dBtheta_dr = difference_r(xx, i, j, k, 32);
    double dBr_dtheta = difference_theta(xx, i, j, k, 31);

    uu[k][j][i].fx[20]=(cot_div_r[yj][xi]*xx[k][j][i].fx[33]+dBphi_dtheta/rr[i] - dBtheta_dphi/rsin[yj][xi])/mu0;
    uu[k][j][i].fx[21]=(dBr_dphi/rsin[yj][xi] - dBphi_dr - xx[k][j][i].fx[31]/rr[i])/mu0;
    uu[k][j][i].fx[22]=(xx[k][j][i].fx[32]/rr[i]+dBtheta_dr - dBr_dtheta/rr[i])/mu0;
}
