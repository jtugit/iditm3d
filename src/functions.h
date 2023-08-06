#include "funcdef.h"
#include "operators.h"
#include "electric_field.h"

inline double functions(Field ***xx, Field ***xn, Field ***uu, int i, int j, int k, int xi, int yj, int zk, int s)
{
    double nsmore, ne, uir, uith, uiph, ui[3], ff;
    double Ec_VxB[3], Ec_gradPe[3], rhos[3], rhossum_nusq_rhon[3], rhon;
    double   usr[3], usth[3], usph[3], unr, unth, unph;
    double   Br, Bt, Bp, TiO, TiH, TiHe, Te, Tn, neme, ui_un_sq[3], uOi_uHi_sq, uOi_uHei_sq, uHi_uHei_sq;
    const double one3rd=1.0/3.0;

    ne = uu[k][j][i].fx[6];
    for (int ss = 0; ss < 3; ss++) rhos[ss] = xn[k][j][i].fx[ss]*ams[ss];

    usr[0]=xx[k][j][i].fx[7]; usth[0]=xx[k][j][i].fx[8]; usph[0]=xx[k][j][i].fx[9];
    usr[1]=xx[k][j][i].fx[10]; usth[1]=xx[k][j][i].fx[11]; usph[1]=xx[k][j][i].fx[12];
    usr[2]=xx[k][j][i].fx[13]; usth[2]=xx[k][j][i].fx[14]; usph[2]=xx[k][j][i].fx[15];
    unr=xx[k][j][i].fx[27]; unth=xx[k][j][i].fx[28]; unph=xx[k][j][i].fx[29];

    nsmore=xn[k][j][i].fx[0]+xn[k][j][i].fx[3]+xn[k][j][i].fx[4]+xn[k][j][i].fx[5]+xn[k][j][i].fx[6];
    uir=(nsmore*usr[0]+xn[k][j][i].fx[1]*usr[1]+xn[k][j][i].fx[2]*usr[2])/ne;
    uith=(nsmore*usth[0]+xn[k][j][i].fx[1]*usth[1]+xn[k][j][i].fx[2]*usth[2])/ne;
    uiph=(nsmore*usph[0]+xn[k][j][i].fx[1]*usph[1]+xn[k][j][i].fx[2]*usph[2])/ne;
 
    //background magnetic field in normal spherical components
    Br=uu[k][j][i].fx[0]; Bt=uu[k][j][i].fx[1]; Bp=uu[k][j][i].fx[2];

    TiO=xx[k][j][i].fx[16]; TiH=xx[k][j][i].fx[17]; TiHe=xx[k][j][i].fx[18];
    Te=xx[k][j][i].fx[19]; Tn=xx[k][j][i].fx[30];
    unr=xx[k][j][i].fx[27]; unth=unr=xx[k][j][i].fx[27]; unph=unr=xx[k][j][i].fx[27];

    neme = uu[k][j][i].fx[6]*ame;

    ui_un_sq[0]=(usr[0]-unr)*(usr[0]-unr)+(usth[0]-unth)*(usth[0]-unth)+(usph[0]-unph)*(usph[0]-unph);
    ui_un_sq[1]=(usr[1]-unr)*(usr[1]-unr)+(usth[1]-unth)*(usth[1]-unth)+(usph[1]-unph)*(usph[1]-unph);
    ui_un_sq[2]=(usr[2]-unr)*(usr[2]-unr)+(usth[2]-unth)*(usth[2]-unth)+(usph[2]-unph)*(usph[2]-unph);
    uOi_uHi_sq =(usr[0]-usr[1])*(usr[0]-usr[1])+(usth[0]-usth[1])*(usth[0]-usth[1])+(usph[0]-usph[1])*(usph[0]-usph[1]);
    uOi_uHei_sq=(usr[0]-usr[2])*(usr[0]-usr[2])+(usth[0]-usth[2])*(usth[0]-usth[2])+(usph[0]-usph[2])*(usph[0]-usph[2]);
    uHi_uHei_sq=(usr[1]-usr[2])*(usr[1]-usr[2])+(usth[1]-usth[2])*(usth[1]-usth[2])+(usph[1]-usph[2])*(usph[1]-usph[2]);

    rhon = 0.0;
    for (int ss = 20; ss < 27; ss++) rhon += xn[k][j][i].fx[ss]*ams[ss-20];

    rhossum_nusq_rhon[0]=( rhos[0]*nust[zk][yj][xi][11]+rhos[3]*nust[zk][yj][xi][30]
                          +rhos[4]*nust[zk][yj][xi][32]+rhos[5]*nust[zk][yj][xi][34]
                          +rhos[6]*nust[zk][yj][xi][36])/rhon;
    rhossum_nusq_rhon[1]=rhos[1]*nust[zk][yj][xi][19]/rhon;
    rhossum_nusq_rhon[2]=rhos[2]*nust[zk][yj][xi][27]/rhon;

//------------ for O+ ---------------------------------------------------------------------------------------
    //----------- uO+_r equation
    if (s == 7)
        ff= xx[k][j][i].fx[7] - xn[k][j][i].fx[7]
           +dt_half*( qms[0]*((uith-usth[0])*Bp-(uiph-usph[0])*Bt)
                     +nust[zk][yj][xi][11]*(usr[0]-unr)+nust[zk][yj][xi][7]*(usr[0]-usr[1])
                     +nust[zk][yj][xi][8]*(usr[0]-usr[2]));   //implicit part C-N scheme

    //----------- uO+_{theta} equation
    else if (s == 8)
        ff= xx[k][j][i].fx[8] - xn[k][j][i].fx[8]
           +dt_half*( qms[0]*((uiph-usph[0])*Br-(uir-usr[0])*Bp)
                     +nust[zk][yj][xi][11]*(usth[0]-unth)+nust[zk][yj][xi][7]*(usth[0]-usth[1])
                     +nust[zk][yj][xi][8]*(usth[0]-usth[2]));

    //----------- uO+_{phi} equation
    else if (s == 9)
        ff= xx[k][j][i].fx[9] - xn[k][j][i].fx[9]
           +dt_half*( qms[0]*((uir-usr[0])*Bt-(uith-usth[0])*Br)
                     +nust[zk][yj][xi][11]*(usph[0]-unph)+nust[zk][yj][xi][7]*(usph[0]-usph[1])
                     +nust[zk][yj][xi][8]*(usph[0]-usph[2]));

    //----------- uH+_r equation
    else if (s == 10)
        ff= xx[k][j][i].fx[10] - xn[k][j][i].fx[10]
           +dt_half*( qms[1]*((uith-usth[1])*Bp-(uiph-usph[1])*Bt)
                     +nust[zk][yj][xi][19]*(usr[1]-unr)+nust[zk][yj][xi][14]*(usr[1]-usr[0])
                     +nust[zk][yj][xi][15]*(usr[1]-usr[2]));   //implicit part C-N scheme

    //----------- uH+_{theta} equation
    else if (s == 11)
        ff= xx[k][j][i].fx[11] - xn[k][j][i].fx[11]
           +dt_half*( qms[1]*((uiph-usph[1])*Br-(uir-usr[1])*Bp)
                     +nust[zk][yj][xi][19]*(usth[1]-unth)+nust[zk][yj][xi][14]*(usth[1]-usth[0])
                     +nust[zk][yj][xi][15]*(usth[1]-usth[2]));

    //----------- uH+_{phi} equation
    else if (s == 12)
        ff= xx[k][j][i].fx[12] - xn[k][j][i].fx[12]
           +dt_half*( qms[1]*((uir-usr[1])*Bt-(uith-usth[1])*Br)
                     +nust[zk][yj][xi][19]*(usph[1]-unph)+nust[zk][yj][xi][14]*(usph[1]-usph[0])
                     +nust[zk][yj][xi][15]*(usph[1]-usph[2]));

//------------ for He+ ---------------------------------------------------------------------------------------
    //----------- uHe+_r equation
    else if (s == 13)
        ff= xx[k][j][i].fx[13] - xn[k][j][i].fx[13]
           +dt_half*( qms[2]*((uith-usth[2])*Bp-(uiph-usph[2])*Bt)
                     +nust[zk][yj][xi][27]*(usr[2]-unr)+nust[zk][yj][xi][22]*(usr[2]-usr[0])
                     +nust[zk][yj][xi][23]*(usr[2]-usr[1]));   //implicit part C-N scheme

    //----------- uHe+_{theta} equation
    else if (s == 14)
        ff= xx[k][j][i].fx[14] - xn[k][j][i].fx[14]
           +dt_half*( qms[2]*((uiph-usph[2])*Br-(uir-usr[2])*Bp)+nust[zk][yj][xi][22]*(usth[2]-usth[0])
                     +nust[zk][yj][xi][23]*(usth[2]-usth[1]));

    //----------- uHe+_{phi} equation
    else if (s == 15)
        ff= xx[k][j][i].fx[15] - xn[k][j][i].fx[15]
           +dt_half*( qms[2]*((uir-usr[2])*Bt-(uith-usth[2])*Br)
                     +nust[zk][yj][xi][27]*(usph[2]-unph)+nust[zk][yj][xi][22]*(usph[2]-usph[0])
                     +nust[zk][yj][xi][23]*(usph[2]-usph[1]));
  
    else if (s == 16)
        //----------- TO+ equation
        ff= xx[k][j][i].fx[16] - xn[k][j][i].fx[16]
           +dt*( ams[0]*( nust[zk][yj][xi][12]*(TiO-Tn)+nust[zk][yj][xi][9]*(TiO-TiH)
                         +nust[zk][yj][xi][10]*(TiO-TiHe))
                +neme/rhos[0]*nust[zk][yj][xi][0]*(TiO-Te)
                -one3rd*ams[0]*( nust[zk][yj][xi][13]*ui_un_sq[0]
                                +ams[1]*nust[zk][yj][xi][9] *uOi_uHi_sq
                                +ams[2]*nust[zk][yj][xi][10]*uOi_uHei_sq)); //implicit part

    else if (s == 17)
       //----------- TH+ equation
        ff= xx[k][j][i].fx[17] - xn[k][j][i].fx[17]
           +dt*( ams[1]*( nust[zk][yj][xi][20]*(TiH-Tn)
                +nust[zk][yj][xi][16]*(TiH-TiO)+nust[zk][yj][xi][17]*(TiH-TiHe))
           +neme/rhos[1]*nust[zk][yj][xi][1]*(TiH-Te)
           -one3rd*ams[1]*( nust[zk][yj][xi][21]*ui_un_sq[1]
                           +nust[zk][yj][xi][18]*uOi_uHi_sq
                           +ams[2]*nust[zk][yj][xi][17]*uHi_uHei_sq)); //implicit part

    else if (s == 18)
        //----------- THe+ equation
        ff= xx[k][j][i].fx[18] - xn[k][j][i].fx[18]
           +dt*( ams[2]*( nust[zk][yj][xi][28]*(TiHe-Tn)+nust[zk][yj][xi][24]*(TiHe-TiO)
                         +nust[zk][yj][xi][25]*(TiHe-TiH))
                +neme/rhos[2]*nust[zk][yj][xi][2]*(TiHe-Te)
                -one3rd*ams[2]*( nust[zk][yj][xi][29]*ui_un_sq[2]
                                +nust[zk][yj][xi][26]*uOi_uHei_sq
                                +ams[1]*nust[zk][yj][xi][25]*uHi_uHei_sq)); //implicit part

    else if (s == 19)
//----------- Te equation --------------------------------------------------------------------------------
        ff= xx[k][j][i].fx[19] - xn[k][j][i].fx[19]
           +dt*( ame*( nust[zk][yj][xi][4]*(Te-TiO)+nust[zk][yj][xi][1]/ams[1]*(Te-TiH)
                      +nust[zk][yj][xi][2]/ams[2]*(Te-TiHe)+nust[zk][yj][xi][6]*(Te-Tn))
                -one3rd*ame*nust[zk][yj][xi][5]*( (uir-unr)*(uir-unr)+(uith-unth)*(uith-unth)
                                                 +(uiph-unph)*(uiph-unph)));
    else if (s == 27)
        //----------- un_{r} equation
        ff= xx[k][j][i].fx[27] - xn[k][j][i].fx[27]
           +dt_half*( rhossum_nusq_rhon[0]*(unr-usr[0])+rhossum_nusq_rhon[1]*(unr-usr[1])
                     +rhossum_nusq_rhon[2]*(unr-usr[2]));  //implicit part

    else if (s == 28)
        //----------- un_{theta} equation
        ff= xx[k][j][i].fx[28] - xn[k][j][i].fx[28]
           +dt_half*( rhossum_nusq_rhon[0]*(unth-usth[0])+rhossum_nusq_rhon[1]*(unth-usth[1])
                     +rhossum_nusq_rhon[2]*(unth-usth[2]));  //implicit part

    else if (s == 29)
        //----------- un_{phi} equation
        ff= xx[k][j][i].fx[29] - xn[k][j][i].fx[29]
           +dt_half*( rhossum_nusq_rhon[0]*(unph-usph[0])+rhossum_nusq_rhon[1]*(unph-usph[1])
                     +rhossum_nusq_rhon[2]*(unph-usph[2]));  //implicit part

    //----------- electric field equation
    else if (s >= 34 and s <= 36) {
        ui[0]=uir; ui[1]=uith; ui[2]=uiph;
        electric_field_vxB(xn, uu, ui, i, j, k, xi, yj, zk, Ec_VxB);
        E_gradPe(xn, uu, i, j, k, yj, zk, Ec_gradPe);

        double Ex = Ec_VxB[0]-Ec_gradPe[0], Ey = Ec_VxB[1]-Ec_gradPe[1], Ez = Ec_VxB[2]-Ec_gradPe[2];

        if (s ==  34) ff= xx[k][j][i].fx[34]-(Jiv11[zk][yj]*Ex+Jiv21[zk][yj]*Ey+Jiv31[yj]*Ez);

        else if (s ==  35) ff= xx[k][j][i].fx[35]-(Jiv12[zk][yj][xi]*Ex+Jiv22[zk][yj][xi]*Ey+Jiv32[yj][xi]*Ez);

        else if (s ==  36) ff= xx[k][j][i].fx[36]-(Jiv13[zk][yj][xi]*Ex+Jiv23[zk][yj][xi]*Ey);
    }

    return ff;
}
