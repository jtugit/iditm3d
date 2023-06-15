#include <cmath>

#include "param.h"

inline void collision_freq(Field ***xx, double nis[], double nns[], double ne, int i, int j, int k, 
    int xi, int yj, int zk, double nqd[])
{
    int    s, t, tt;
    double ni[7], nn[7], Te, Te12, nu_st[14], amt, Dst;
    double Ti, Tn, Ti12, Ti32, logTi, Tr, logTr, Tr12, logTrOi, Tr12Oi;
    const double n00=n0*1.0e-6;

    for (s = 0; s < sl; s++) {
        ni[s] = nis[s]*n00; nn[s] = nns[s]*n00; //density in cm^{-3}
    }
    ne=ne*n00;

/*--------------------------------------------------*/
/*---------------- collision frequencies -----------*/
/*--------------------------------------------------*/
    //electron Coulomb collision frequencies
    Te=xx[k][j][i].fx[19]*T0;
    Te12=sqrt(Te);

    //electron - ions collisions
    for (s = 0; s < sl; s++) nu_st[s]=coe[0]*ni[s]/(Te*Te12);
    nust[zk][yj][xi][0] = nu_st[0];    //e - O+ collision frequency
    nust[zk][yj][xi][1] = nu_st[1];    //e - H+ collision frequency
    nust[zk][yj][xi][2] = nu_st[2];    //e - He+ collision frequency

    //nu_{e,O+} + nu_{e,O2+} + nu_{e,N2+} + nu_{e,NO+} + nu_{e,N+}
    nust[zk][yj][xi][3] = nu_st[0]+nu_st[3]+nu_st[4]+nu_st[5]+nu_st[6];

    //nu_{e,O+}/m_O+ + nu_{e,O2+}/m_O2+ + nu_{e,N2+}/m_N2+ + nu_{e,NO+}/m_NO+ + nu_{e,N+}/m_N+
    nust[zk][yj][xi][4]= nu_st[0]/ams[0]+nu_st[3]/ams[3]+nu_st[4]/ams[4]
                        +nu_st[5]/ams[5]+nu_st[6]/ams[6];

    /* electron - neutral collision frequencies */
    nu_st[7] =coe[1]*nn[0]*(1.0+5.7e-4*Te)*Te12;    //e - O
    nu_st[8] =coe[2]*nn[1]*(1.0-1.35e-4*Te)*Te12;    //e - H
    nu_st[9] =coe[3]*nn[2]*Te12;                    //e - He
    nu_st[10]=coe[4]*nn[3]*(1.0+3.6e-2*Te12)*Te12; //e - O2 
    nu_st[11]=coe[5]*nn[4]*(1.0-1.21e-4*Te)*Te;    //e - N2
    nu_st[12]=coe[6]*nn[5];                         //e - NO
    nu_st[13]=coe[7]*nn[6];                         //e - N

    //total e - neutral collision frequency
    nust[zk][yj][xi][5]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{e,neutrals}/ms[s]
    nust[zk][yj][xi][6]= nu_st[7]/ams[0]+nu_st[8]/ams[1]+nu_st[9]/ams[2]+nu_st[10]/ams[3]
                        +nu_st[11]/ams[4]+nu_st[12]/ams[5]+nu_st[13]/ams[6];

/* ---- O+ collision frequencies ---------------------------------------------------------*/
    Ti=xx[k][j][i].fx[16]*T0; Tn=xx[k][j][i].fx[30]*T0;
    Ti12=sqrt(Ti); Ti32=Ti*Ti12; logTi=log10(Ti);
    Tr=0.5*(Ti+Tn); logTrOi=log10(Tr); Tr12Oi=sqrt(Tr);

/* O+ Coulomb collision frequencies */
    nu_st[1]=coiO[0]*ni[1]/Ti32;     // O+ - H+
    nu_st[2]=coiO[1]*ni[2]/Ti32;     // O+ - He+
    nu_st[3]=coiO[2]*ni[3]/Ti32;     // O+ - O2+
    nu_st[4]=coiO[3]*ni[4]/Ti32;     // O+ - N2+
    nu_st[5]=coiO[4]*ni[5]/Ti32;     // O+ - NO+
    nu_st[6]=coiO[5]*ni[6]/Ti32;     // O+ - N+

    nust[zk][yj][xi][7] = nu_st[1];  //nu_{O+,H+}
    nust[zk][yj][xi][8] = nu_st[2];  //nu_{O+,He+}

    nust[zk][yj][xi][9] = nu_st[1]/(ams[0]+ams[1]);  //nu_{O+,H+}/(m_O+ + m_H+)
    nust[zk][yj][xi][10]= nu_st[2]/(ams[0]+ams[2]);  //nu_{O+,He+}/(m_O+ + m_He+)

    /* O+ - neutral collision frequencies */
    nu_st[7]=coiO[6]*nn[0]*Tr12Oi*pow(1.0-0.064*logTrOi, 2.0);  // O+ - O
    nu_st[8]=coiO[7]*nn[1]*Ti12*pow(1.0-0.047*logTi, 2.0);  // O+ - H
    nu_st[9]=coiO[8]*nn[2];        // O+ - He
    nu_st[10]=coiO[9]*nn[3];        // O+ - O2
    nu_st[11]=coiO[10]*nn[4];        // O+ - N2
    nu_st[12]=coiO[11]*nn[5];        // O+ - NO
    nu_st[13]=coiO[12]*nn[6];        // O+ - N

    //sum of nu_{O+, neutral species} 
    nust[zk][yj][xi][11]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{O+,q} / (m_O+ + m_q)
    nust[zk][yj][xi][12]= nu_st[7]/(ams[0]+ams[0])+nu_st[8]/(ams[0]+ams[1])
                         +nu_st[9]/(ams[0]+ams[2])+nu_st[10]/(ams[0]+ams[3])
                         +nu_st[11]/(ams[0]+ams[4])+nu_st[12]/(ams[0]+ams[5])
                         +nu_st[13]/(ams[0]+ams[6]);

    //sum of m_q nu_{O+,q} / (m_O+ + m_q)
    nust[zk][yj][xi][13]= ams[0]*nu_st[7]/(ams[0]+ams[0])+ams[1]*nu_st[8]/(ams[0]+ams[1])
                         +ams[2]*nu_st[9]/(ams[0]+ams[2])+ams[3]*nu_st[10]/(ams[0]+ams[3])
                         +ams[4]*nu_st[11]/(ams[0]+ams[4])+ams[5]*nu_st[12]/(ams[0]+ams[5])
                         +ams[6]*nu_st[13]/(ams[0]+ams[6]);

    //quantity for calculating O+ thermo-conductivity
    nqd[1]=0.0;
    for (t = 1; t < sl+sm; t++) {
        if (t < sl) {
            Dst=(3.0*ams[0]*ams[0]-0.2*ams[t]*ams[t]+0.1*ams[0]*ams[t])
                /((ams[0]+ams[t])*(ams[0]+ams[t])); 
            nqd[1] += nu_st[t]*(Dst+1.5*ams[t]/(ams[0]+ams[t]));
        }
        else {
            // ion - neutral interaction Schunk & Nagy eq (4.147b)
            amt=ams[t-7];
            Dst=(3.0*ams[0]*ams[0]+amt*amt+1.6*ams[0]*amt)/((ams[0]+amt)*(ams[0]+amt));
            nqd[1] += nu_st[t]*(Dst+1.5*amt/(ams[0]+amt));
        }
    }
    
/* ---- H+ collision frequencies ---------------------------------------------------------*/
    Ti=xx[k][j][i].fx[17]*T0;
    Ti12=sqrt(Ti); Ti32=Ti*Ti12; logTi=log10(Ti);
    logTr=log10(0.5*(Ti+Tn)); Tr12=sqrt(Tr);

    /* H+ Coulomb collision frequencies */
    nu_st[1]=coiH[0]*ni[0]/Ti32;     // H+ - O+
    nu_st[2]=coiH[1]*ni[2]/Ti32;     // H+ - He+
    nu_st[3]=coiH[2]*ni[3]/Ti32;     // H+ - O2+
    nu_st[4]=coiH[3]*ni[4]/Ti32;     // H+ - N2+
    nu_st[5]=coiH[4]*ni[5]/Ti32;     // H+ - NO+
    nu_st[6]=coiH[5]*ni[6]/Ti32;     // H+ - N+

    //nu_{H+,O+} + nu_{H+,O2+} + nu_{H+,N2+} + nu_{H+,NO+} + nu_{H+,N+}
    nust[zk][yj][xi][14]=nu_st[1]+nu_st[3]+nu_st[4]+nu_st[5]+nu_st[6];
    nust[zk][yj][xi][15]=nu_st[2];  //nu_{H+,He+}

    //sum of nu_{H+,t}/(m_H+ + m_t)
    nust[zk][yj][xi][16]= nu_st[1]/(ams[1]+ams[0])+nu_st[3]/(ams[1]+ams[3])
                         +nu_st[4]/(ams[1]+ams[4])+nu_st[5]/(ams[1]+ams[5])
                         +nu_st[6]/(ams[1]+ams[6]);

    nust[zk][yj][xi][17]= nu_st[2]/(ams[1]+ams[2]);

    //sum of m_t nu_{H+,t}/(m_H+ + m_t)
    nust[zk][yj][xi][18]= ams[0]*nu_st[1]/(ams[1]+ams[0])+ams[3]*nu_st[3]/(ams[1]+ams[3])
                         +ams[4]*nu_st[4]/(ams[1]+ams[4])+ams[5]*nu_st[5]/(ams[1]+ams[5])
                         +ams[6]*nu_st[6]/(ams[1]+ams[6]);

/* H+ - neutral collision frequencies */
    nu_st[7]= coiH[6]*nn[0]*Ti12*pow(1.0-0.047*logTi, 2.0);  // H+ - O
    nu_st[8]= coiH[7]*nn[1]*Tr12*pow(1.0-0.083*logTr, 2.0);  // H+ - H
    nu_st[9]= coiH[8]*nn[2];       // H+ - He
    nu_st[10]= coiH[9]*nn[3];       // H+ - O2
    nu_st[11]= coiH[10]*nn[4];       // H+ - N2
    nu_st[12]= coiH[11]*nn[5];      // H+ - NO
    nu_st[13]= coiH[12]*nn[6];       // H+ - N

    //sum of nu_{H+, neutral species} 
    nust[zk][yj][xi][19]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{H+,q} / (m_H+ + m_q)
    nust[zk][yj][xi][20]= nu_st[7]/(ams[1]+ams[0])+nu_st[8]/(ams[1]+ams[1])
                         +nu_st[9]/(ams[1]+ams[2])+nu_st[10]/(ams[1]+ams[3])
                         +nu_st[11]/(ams[1]+ams[4])+nu_st[12]/(ams[1]+ams[5])
                         +nu_st[13]/(ams[1]+ams[6]);

    //sum of m_q nu_{H+,q} / (m_H+ + m_q)
    nust[zk][yj][xi][21]= ams[0]*nu_st[7]/(ams[1]+ams[0])+ams[1]*nu_st[8]/(ams[1]+ams[1])
                         +ams[2]*nu_st[9]/(ams[1]+ams[2])+ams[3]*nu_st[10]/(ams[1]+ams[3])
                         +ams[4]*nu_st[11]/(ams[1]+ams[4])+ams[5]*nu_st[12]/(ams[1]+ams[5])
                         +ams[6]*nu_st[13]/(ams[1]+ams[6]);

    //quantity for calculating H+ thermo-conductivity
    nqd[2]=0.0;
    for (t = 1; t < sl+sm; t++) {
        if (t < sl) {
            if (t <= 1) tt=t-1; else tt=t;
            Dst=(3.0*ams[1]*ams[1]-0.2*ams[tt]*ams[tt]+0.1*ams[1]*ams[tt])
                /((ams[1]+ams[tt])*(ams[1]+ams[tt])); 
            nqd[2] += nu_st[t]*(Dst+1.5*ams[tt]/(ams[1]+ams[tt]));
        }
        else {
            // ion - neutral interaction Schunk & Nagy eq (4.147b)
            amt=ams[t-7];
            Dst=(3.0*ams[1]*ams[1]+amt*amt+1.6*ams[1]*amt)/((ams[1]+amt)*(ams[1]+amt));
            nqd[2] += nu_st[t]*(Dst+1.5*amt/(ams[1]+amt));
        }
    }

/* ---- He+ collision frequencies ---------------------------------------------------------*/
    Ti=xx[k][j][i].fx[18]*T0;
    logTr=log10(0.5*(Ti+Tn)); Tr12=sqrt(Tr);

    /* He+ Coulomb collision frequencies */
    nu_st[1]=coiHe[0]*ni[0]/Ti32;      // He+ - O+
    nu_st[2]=coiHe[1]*ni[1]/Ti32;      // He+ - H+
    nu_st[3]=coiHe[2]*ni[3]/Ti32;      // He+ - O2+
    nu_st[4]=coiHe[3]*ni[4]/Ti32;      // He+ - N2+
    nu_st[5]=coiHe[4]*ni[5]/Ti32;      // He+ - NO+
    nu_st[6]=coiHe[5]*ni[6]/Ti32;      // He+ - N+

    //nu_{He+,O+} + nu_{He+,O2+} + nu_{He+,N2+} + nu_{He+,NO+} + nu_{He+,N+}
    nust[zk][yj][xi][22]=nu_st[1]+nu_st[3]+nu_st[4]+nu_st[5]+nu_st[6];
    nust[zk][yj][xi][23]=nu_st[2];  //nu_{He+,H+}

    //sum of nu_{He+,t}/(m_He+ + m_t)
    nust[zk][yj][xi][24]= nu_st[1]/(ams[2]+ams[0])+nu_st[3]/(ams[2]+ams[3])
                         +nu_st[4]/(ams[2]+ams[4])+nu_st[5]/(ams[2]+ams[5])
                         +nu_st[6]/(ams[2]+ams[6]);

    nust[zk][yj][xi][25]= nu_st[2]/(ams[2]+ams[1]);

    //sum of m_t nu_{He+,t}/(m_He+ + m_t)
    nust[zk][yj][xi][26]= ams[0]*nu_st[1]/(ams[2]+ams[0])+ams[3]*nu_st[3]/(ams[2]+ams[3])
                         +ams[4]*nu_st[4]/(ams[2]+ams[4])+ams[5]*nu_st[5]/(ams[2]+ams[5])
                         +ams[6]*nu_st[6]/(ams[2]+ams[6]);

/* He+ - neutral collision frequencies */
    nu_st[7]= coiHe[6]*nn[0];         // He+ - O
    nu_st[8]= coiHe[7]*nn[1];        // He+ - H
    nu_st[9]= coiHe[8]*nn[2]*Tr12*pow(1.0-0.083*logTr, 2.0); // He+ - He
    nu_st[10]= coiHe[9]*nn[3];         // He+ - O2
    nu_st[11]= coiHe[10]*nn[4];         // He+ - N2
    nu_st[12]= coiHe[11]*nn[5];        // He+ - NO
    nu_st[13]= coiHe[12]*nn[6];         // He+ - N

    //sum of nu_{He+, neutral species} 
    nust[zk][yj][xi][27]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{He+,q} / (m_He+ + m_q)
    nust[zk][yj][xi][28]= nu_st[7]/(ams[2]+ams[0])+nu_st[8]/(ams[2]+ams[1])
                         +nu_st[9]/(ams[2]+ams[2])+nu_st[10]/(ams[2]+ams[3])
                         +nu_st[11]/(ams[2]+ams[4])+nu_st[12]/(ams[2]+ams[5])
                         +nu_st[13]/(ams[2]+ams[6]);

    //sum of m_q nu_{He+,q} / (m_He+ + m_q)
    nust[zk][yj][xi][29]= ams[0]*nu_st[7]/(ams[2]+ams[0])+ams[1]*nu_st[8]/(ams[2]+ams[1])
                         +ams[2]*nu_st[9]/(ams[2]+ams[2])+ams[3]*nu_st[10]/(ams[2]+ams[3])
                         +ams[4]*nu_st[11]/(ams[2]+ams[4])+ams[5]*nu_st[12]/(ams[2]+ams[5])
                         +ams[6]*nu_st[13]/(ams[2]+ams[6]);

    //quantity for calculating H+ thermo-conductivity
    nqd[3]=0.0;
    for (t = 1; t < sl+sm; t++) {
        if (t < sl) {
            if (t <= 2) tt=t-1; else tt=t;
            Dst=(3.0*ams[2]*ams[2]-0.2*ams[tt]*ams[tt]+0.1*ams[2]*ams[tt])
                /((ams[2]+ams[tt])*(ams[2]+ams[tt])); 
            nqd[3] += nu_st[t]*(Dst+1.5*ams[tt]/(ams[2]+ams[tt]));
        }
        else {
            // ion - neutral interaction Schunk & Nagy eq (4.147b)
            amt=ams[t-7];
            Dst=(3.0*ams[2]*ams[2]+amt*amt+1.6*ams[2]*amt)/((ams[2]+amt)*(ams[2]+amt));
            nqd[3] += nu_st[t]*(Dst+1.5*amt/(ams[2]+amt));
        }
    }

/* ---- O2+ collision frequencies --------------------------------------------------------*/
    /* O2+ - neutral collision frequencies */
    nu_st[7]= 2.31e-10*nn[0];              // O2+ - O
    nu_st[8]= 6.50e-11*nn[1];              // O2+ - H
    nu_st[9]= 7.00e-11*nn[2];              // O2+ - He
    nu_st[10]= 2.59e-11*nn[3]*Tr12Oi*pow(1.0-0.073*logTrOi, 2.0); // O2+ - O2
    nu_st[11]= 4.13e-10*nn[4];              // O2+ - N2
    nu_st[12]= 2.69e-11*nn[5];              // O2+ - NO
    nu_st[13]= 2.64e-10*nn[6];              // O2+ - N

    //sum of nu_{O2+, neutral species} 
    nust[zk][yj][xi][30]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{O2+,q} / (m_O2+ + m_q)
    nust[zk][yj][xi][31]= nu_st[7]/(ams[3]+ams[0])+nu_st[8]/(ams[3]+ams[1])
                         +nu_st[9]/(ams[3]+ams[2])+nu_st[10]/(ams[3]+ams[3])
                         +nu_st[11]/(ams[3]+ams[4])+nu_st[12]/(ams[3]+ams[5])
                         +nu_st[13]/(ams[3]+ams[6]);

/* ---- N2+ collision frequencies --------------------------------------------------------*/
    /* N2+ - neutral collision frequencies */
    nu_st[7]= 2.58e-10*nn[0];                 // N2+ - O
    nu_st[8]= 7.40e-11*nn[1];                 // N2+ - H
    nu_st[9]= 7.90e-11*nn[2];                 // N2+ - He
    nu_st[10]= 4.49e-10*nn[3];                 // N2+ - O2
    nu_st[11]= 5.14e-11*nn[4]*Tr12Oi*pow(1.0-0.069*logTrOi, 2.0);    // N2+ - N2
    nu_st[12]= 2.69e-11*nn[5];                // N2+ - NO
    nu_st[13]= 2.95e-10*nn[6];                // N2+ - N

    //sum of nu_{N2+, neutral species} 
    nust[zk][yj][xi][32]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{N2+,q} / (m_N2+ + m_q)
    nust[zk][yj][xi][33]= nu_st[7]/(ams[4]+ams[0])+nu_st[8]/(ams[4]+ams[1])
                         +nu_st[9]/(ams[4]+ams[2])+nu_st[10]/(ams[4]+ams[3])
                         +nu_st[11]/(ams[4]+ams[4])+nu_st[12]/(ams[4]+ams[5])
                         +nu_st[13]/(ams[4]+ams[6]);

/* ---- NO+ collision frequencies -------------------------------------------------------*/
    /* NO+ - neutral collision frequencies nu_{NO+,n} */
    nu_st[7]= 2.44e-10*nn[0];                   // NO+ - O
    nu_st[8]= 6.90e-11*nn[1];                   // NO+ - H
    nu_st[9]= 7.40e-11*nn[2];                   // NO+ - He
    nu_st[10]= 4.27e-10*nn[3];                   // NO+ - O2
    nu_st[11]= 4.34e-10*nn[4];                   // NO+ - N2
    nu_st[12]= 2.69e-11*nn[5];                   // NO+ - NO
    nu_st[13]= 2.79e-10*nn[6];                   // NO+ - N

    //sum of nu_{NO+, neutral species} 
    nust[zk][yj][xi][34]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{NO+,q} / (m_NO+ + m_q)
    nust[zk][yj][xi][35]= nu_st[7]/(ams[5]+ams[0])+nu_st[8]/(ams[5]+ams[1])
                         +nu_st[9]/(ams[5]+ams[2])+nu_st[10]/(ams[5]+ams[3])
                         +nu_st[11]/(ams[5]+ams[4])+nu_st[12]/(ams[5]+ams[5])
                         +nu_st[13]/(ams[5]+ams[6]);

/* ---- N+ collision frequencies -------------------------------------------------------*/
    /* N+ - neutral collision frequencies */
    nu_st[7]= 4.42e-10*nn[0];                  // N+ - O
    nu_st[8]= 1.45e-10*nn[1];                  // N+ - H
    nu_st[9]= 1.49e-10*nn[2];                  // N+ - He
    nu_st[10]= 7.25e-10*nn[3];                  // N+ - O2
    nu_st[11]= 7.47e-10*nn[4];                  // N+ - N2
    nu_st[12]= 2.69e-11*nn[5];                  // N+ - NO
    nu_st[13]= 3.83e-11*nn[6]*Tr12Oi*pow(1.0-0.063*logTrOi, 2.0);  // N+ - N

    //sum of nu_{N+, neutral species} 
    nust[zk][yj][xi][36]=nu_st[7]+nu_st[8]+nu_st[9]+nu_st[10]+nu_st[11]+nu_st[12]+nu_st[13];

    //sum of nu_{N+,q} / (m_N+ + m_q)
    nust[zk][yj][xi][37]= nu_st[7]/(ams[6]+ams[0])+nu_st[8]/(ams[6]+ams[1])
                         +nu_st[9]/(ams[6]+ams[2])+nu_st[10]/(ams[6]+ams[3])
                         +nu_st[11]/(ams[6]+ams[4])+nu_st[12]/(ams[6]+ams[5])
                         +nu_st[13]/(ams[6]+ams[6]);
}