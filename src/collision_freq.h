#include <cmath>

#include "param.h"

inline void collision_freq(Field ***xx, int i, int j, int k, int xi, int yj, int zk)
{
    int    s;
    double ni[7], nn[7], ne=0.0, Nn=0.0, Te, Te12;
    double Ti, Tn, Ti12, Ti32, logTi, Tr, logTr, Tr12, logTrOi, Tr12Oi;
    const double n00=n0*1.0e-6;

    for (s = 0; s < sl; s++) {
        ni[s] = exp(xx[k][j][i].fx[s])*n00;   //density in cm^{-3}
        ne += ni[s];

        nn[s] = exp(xx[k][j][i].fx[120+s])*n00;
        Nn += nn[s];
    }

/*--------------------------------------------------*/
/*---------------- collision frequencies -----------*/
/*--------------------------------------------------*/
    //electron Coulomb collision frequencies
    Te=xx[k][j][i].fx[19]*T0;
    Te12=sqrt(Te);

    //electron - ions collisions
    for (s = 0; s < sl; s++) nust[zk][yj][xi][s]=coe[0]*ni[s];

    /* electron - neutral collision frequencies */
    nust[zk][yj][xi][7] =coe[1]*nn[0]*(1.0+5.7e-4*Te)*Te12;    //e - O
    nust[zk][yj][xi][8] =coe[2]*nn[1]*(1.0-1.35e-4*Te)*Te12;    //e - H
    nust[zk][yj][xi][9] =coe[3]*nn[2]*Te12;                    //e - He
    nust[zk][yj][xi][10]=coe[4]*nn[3]*(1.0+3.6e-2*Te12)*Te12; //e - O2 
    nust[zk][yj][xi][11]=coe[5]*nn[4]*(1.0-1.21e-4*Te)*Te;    //e - N2
    nust[zk][yj][xi][12]=coe[6]*nn[5];                         //e - NO
    nust[zk][yj][xi][13]=coe[7]*nn[6];                         //e - N

/* ---- O+ collision frequencies ---------------------------------------------------------*/
    Ti=xx[k][j][i].fx[16]*T0; Tn=xx[k][j][i].fx[30]*T0;
    Ti12=sqrt(Ti); Ti32=Ti*Ti12; logTi=log10(Ti);
    Tr=0.5*(Ti+Tn); logTrOi=log10(Tr); Tr12Oi=sqrt(Tr);

/* O+ Coulomb collision frequencies */
    nust[zk][yj][xi][14]=ne*ame/(ni[0]*ams[0])*nust[zk][yj][xi][0]; // O+ - e
    nust[zk][yj][xi][15]=coiO[0]*ni[1]/Ti32;     // O+ - H+
    nust[zk][yj][xi][16]=coiO[1]*ni[2]/Ti32;     // O+ - He+
    nust[zk][yj][xi][17]=coiO[2]*ni[3]/Ti32;     // O+ - O2+
    nust[zk][yj][xi][18]=coiO[3]*ni[4]/Ti32;     // O+ - N2+
    nust[zk][yj][xi][19]=coiO[4]*ni[5]/Ti32;     // O+ - NO+
    nust[zk][yj][xi][20]=coiO[5]*ni[6]/Ti32;     // O+ - N+

    /* O+ - neutral collision frequencies */
    nust[zk][yj][xi][21]=coiO[6]*nn[0]*Tr12Oi*pow(1.0-0.064*logTrOi, 2.0);  // O+ - O
    nust[zk][yj][xi][22]=coiO[7]*nn[1]*Ti12*pow(1.0-0.047*logTi, 2.0);  // O+ - H
    nust[zk][yj][xi][23]=coiO[8]*nn[2];        // O+ - He
    nust[zk][yj][xi][24]=coiO[9]*nn[3];        // O+ - O2
    nust[zk][yj][xi][25]=coiO[10]*nn[4];        // O+ - N2
    nust[zk][yj][xi][26]=coiO[11]*nn[5];        // O+ - NO
    nust[zk][yj][xi][27]=coiO[12]*nn[6];        // O+ - N

/* ---- H+ collision frequencies ---------------------------------------------------------*/
    Ti=xx[k][j][i].fx[17]*T0;
    Ti12=sqrt(Ti); Ti32=Ti*Ti12; logTi=log10(Ti);
    logTr=log10(0.5*(Ti+Tn)); Tr12=sqrt(Tr);

    /* H+ Coulomb collision frequencies */
    nust[zk][yj][xi][28]=ne*ame/(ni[1]*ams[1])*nust[zk][yj][xi][1]; // H+ - e
    nust[zk][yj][xi][29]=coiH[0]*ni[0]/Ti32;     // H+ - O+
    nust[zk][yj][xi][30]=coiH[1]*ni[2]/Ti32;     // H+ - He+
    nust[zk][yj][xi][31]=coiH[2]*ni[3]/Ti32;     // H+ - O2+
    nust[zk][yj][xi][32]=coiH[3]*ni[4]/Ti32;     // H+ - N2+
    nust[zk][yj][xi][33]=coiH[4]*ni[5]/Ti32;     // H+ - NO+
    nust[zk][yj][xi][34]=coiH[5]*ni[6]/Ti32;     // H+ - N+

/* H+ - neutral collision frequencies */
    nust[zk][yj][xi][35]= coiH[6]*nn[0]*Ti12*pow(1.0-0.047*logTi, 2.0);  // H+ - O
    nust[zk][yj][xi][36]= coiH[7]*nn[1]*Tr12*pow(1.0-0.083*logTr, 2.0);  // H+ - H
    nust[zk][yj][xi][37]= coiH[8]*nn[2];       // H+ - He
    nust[zk][yj][xi][38]= coiH[9]*nn[3];       // H+ - O2
    nust[zk][yj][xi][39]= coiH[10]*nn[4];       // H+ - N2
    nust[zk][yj][xi][40]= coiH[11]*nn[5];      // H+ - NO
    nust[zk][yj][xi][41]= coiH[12]*nn[6];       // H+ - N

/* ---- He+ collision frequencies ---------------------------------------------------------*/
    Ti=xx[k][j][i].fx[18]*T0;
    logTr=log10(0.5*(Ti+Tn)); Tr12=sqrt(Tr);

    /* He+ Coulomb collision frequencies */
    nust[zk][yj][xi][42]=ne*ame/(ni[2]*ams[2])*nust[zk][yj][xi][2];   // He+ - e
    nust[zk][yj][xi][43]=coiHe[0]*ni[0]/Ti32;      // He+ - O+
    nust[zk][yj][xi][44]=coiHe[1]*ni[1]/Ti32;      // He+ - H+
    nust[zk][yj][xi][45]=coiHe[2]*ni[3]/Ti32;      // He+ - O2+
    nust[zk][yj][xi][46]=coiHe[3]*ni[4]/Ti32;      // He+ - N2+
    nust[zk][yj][xi][47]=coiHe[4]*ni[5]/Ti32;      // He+ - NO+
    nust[zk][yj][xi][48]=coiHe[5]*ni[6]/Ti32;      // He+ - N+

/* He+ - neutral collision frequencies */
    nust[zk][yj][xi][49]= coiHe[6]*nn[0];         // He+ - O
    nust[zk][yj][xi][50]= coiHe[7]*nn[1];        // He+ - H
    nust[zk][yj][xi][51]= coiHe[8]*nn[2]*Tr12*pow(1.0-0.083*logTr, 2.0); // He+ - He
    nust[zk][yj][xi][52]= coiHe[9]*nn[3];         // He+ - O2
    nust[zk][yj][xi][53]= coiHe[10]*nn[4];         // He+ - N2
    nust[zk][yj][xi][54]= coiHe[11]*nn[5];        // He+ - NO
    nust[zk][yj][xi][55]= coiHe[12]*nn[6];         // He+ - N

/* ---- O2+ collision frequencies --------------------------------------------------------*/
    /* O2+ Coulomb collision frequencies */
    nust[zk][yj][xi][57]=1.30e-1*ni[0]/Ti32;    // O2+ - O+
    nust[zk][yj][xi][58]=3.90e-2*ni[1]/Ti32;    // O2+ - H+
    nust[zk][yj][xi][59]=7.50e-2*ni[2]/Ti32;    // O2+ - He+
    nust[zk][yj][xi][60]=1.50e-1*ni[4]/Ti32;    // O2+ - N2+
    nust[zk][yj][xi][61]=1.60e-1*ni[5]/Ti32;    // O2+ - NO+
    nust[zk][yj][xi][62]=1.20e-1*ni[6]/Ti32;    // O2+ - N+

    /* O2+ - neutral collision frequencies */
    nust[zk][yj][xi][63]= 2.31e-10*nn[0];              // O2+ - O
    nust[zk][yj][xi][64]= 6.50e-11*nn[1];              // O2+ - H
    nust[zk][yj][xi][65]= 7.00e-11*nn[2];              // O2+ - He
    nust[zk][yj][xi][66]= 2.59e-11*nn[3]*Tr12Oi*pow(1.0-0.073*logTrOi, 2.0); // O2+ - O2
    nust[zk][yj][xi][67]= 4.13e-10*nn[4];              // O2+ - N2
    nust[zk][yj][xi][68]= 2.69e-11*nn[5];              // O2+ - NO
    nust[zk][yj][xi][69]= 2.64e-10*nn[6];              // O2+ - N

/* ---- N2+ collision frequencies --------------------------------------------------------*/
    /* N2+ Coulomb collision frequencies */
    nust[zk][yj][xi][70]=ne*ame/(ni[4]*ams[4])*nust[zk][yj][xi][4];   // N2+ - e
    nust[zk][yj][xi][71]=1.50e-1*ni[0]/Ti32;      // N2+ - O+
    nust[zk][yj][xi][72]=4.50e-2*ni[1]/Ti32;      // N2+ - H+
    nust[zk][yj][xi][73]=8.50e-2*ni[2]/Ti32;      // N2+ - He+
    nust[zk][yj][xi][74]=1.80e-1*ni[3]/Ti32;      // N2+ - O2+
    nust[zk][yj][xi][75]=1.70e-1*ni[5]/Ti32;      // N2+ - NO+
    nust[zk][yj][xi][76]=1.40e-1*ni[6]/Ti32;      // N2+ - N+

    /* N2+ - neutral collision frequencies */
    nust[zk][yj][xi][77]= 2.58e-10*nn[0];                 // N2+ - O
    nust[zk][yj][xi][78]= 7.40e-11*nn[1];                 // N2+ - H
    nust[zk][yj][xi][79]= 7.90e-11*nn[2];                 // N2+ - He
    nust[zk][yj][xi][80]= 4.49e-10*nn[3];                 // N2+ - O2
    nust[zk][yj][xi][81]= 5.14e-11*nn[4]*Tr12Oi*pow(1.0-0.069*logTrOi, 2.0);    // N2+ - N2
    nust[zk][yj][xi][82]= 2.69e-11*nn[5];                // N2+ - NO
    nust[zk][yj][xi][83]= 2.95e-10*nn[6];                // N2+ - N

/* ---- NO+ collision frequencies -------------------------------------------------------*/
    /* NO+ Coulomb collision frequencies */
    nust[zk][yj][xi][84]=ne*ame/(ni[5]*ams[5])*nust[zk][yj][xi][5]; // NO+ - e
    nust[zk][yj][xi][85]=1.40e-1*ni[0]/Ti32;    // NO+ - O+
    nust[zk][yj][xi][86]=4.20e-2*ni[1]/Ti32;    // NO+ - H+
    nust[zk][yj][xi][87]=8.00e-2*ni[2]/Ti32;    // NO+ - H+e
    nust[zk][yj][xi][88]=1.70e-1*ni[3]/Ti32;    // NO+ - O2+
    nust[zk][yj][xi][89]=1.60e-1*ni[4]/Ti32;    // NO+ - N2+
    nust[zk][yj][xi][90]=1.30e-1*ni[6]/Ti32;    // NO+ - N+

    /* NO+ - neutral collision frequencies nu_{NO+,n} */
    nust[zk][yj][xi][91]= 2.44e-10*nn[0];                   // NO+ - O
    nust[zk][yj][xi][92]= 6.90e-11*nn[1];                   // NO+ - H
    nust[zk][yj][xi][93]= 7.40e-11*nn[2];                   // NO+ - He
    nust[zk][yj][xi][94]= 4.27e-10*nn[3];                   // NO+ - O2
    nust[zk][yj][xi][95]= 4.34e-10*nn[4];                   // NO+ - N2
    nust[zk][yj][xi][96]= 2.69e-11*nn[5];                   // NO+ - NO
    nust[zk][yj][xi][97]= 2.79e-10*nn[6];                   // NO+ - N

/* ---- N+ collision frequencies -------------------------------------------------------*/
    /* N+ Coulomb collision frequencies */
    nust[zk][yj][xi][98] =ne*ame/(ni[6]*ams[6])*nust[zk][yj][xi][6]; // NO+ - e
    nust[zk][yj][xi][99] =0.25*ni[0]/Ti32;    // N+ - O+
    nust[zk][yj][xi][100]=0.088*ni[1]/Ti32;   // N+ - H+
    nust[zk][yj][xi][101]=0.16*ni[2]/Ti32;    // N+ - H+e
    nust[zk][yj][xi][102]=0.28*ni[3]/Ti32;    // N+ - O2+
    nust[zk][yj][xi][103]=0.28*ni[4]/Ti32;    // N+ - N2+
    nust[zk][yj][xi][104]=0.28*ni[6]/Ti32;    // N+ - NO+

    /* N+ - neutral collision frequencies */
    nust[zk][yj][xi][105]= 4.42e-10*nn[0];                  // N+ - O
    nust[zk][yj][xi][106]= 1.45e-10*nn[1];                  // N+ - H
    nust[zk][yj][xi][107]= 1.49e-10*nn[2];                  // N+ - He
    nust[zk][yj][xi][108]= 7.25e-10*nn[3];                  // N+ - O2
    nust[zk][yj][xi][109]= 7.47e-10*nn[4];                  // N+ - N2
    nust[zk][yj][xi][110]= 2.69e-11*nn[5];                  // N+ - NO
    nust[zk][yj][xi][111]= 3.83e-11*nn[6]*Tr12Oi*pow(1.0-0.063*logTrOi, 2.0);  // N+ - N
}