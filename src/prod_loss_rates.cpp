/*************************************************************************
  Routine prod_loss_rates
     Evaluate normalized production and loss rates for ion and neutral
 *   species considered. Also evaluate ion and neutral heating rates due
 *   to exothermic reactions

 * Jiannan Tu
 * 5/21/2022
*************************************************************************/
#include <cmath>
using namespace std;

#include "param.h"
#include "funcdef.h"

double pesimp(double ALT,double SZA,double F107,double TE,double TN,double XNE,
        double nO,double nO2,double nN2);

void prod_loss_rates(Field ***xx, Field ***uu, int i, int j, int k, int zk, int yj, int xi,
    double &Qeuv, double &Qphoto, double Ps[], double Ls[])
{
    int    m, s, isza;
    double Xp, euvf, Hn[5], Ch[5], Nnn[5], tao, CCi;
    double coschi, tcos, y, pih, Tn, qi[8];

    double T1, T2, T300, T3002, Tr, Ti, Te;
    double niO, niO2, niN2, niH, niHe, niNO, niN;
    double nO, nO2, nN2, nH, nHe, nNO, nN;
    double kk[29], react[29], photoIon[8];

    double R, R2, epsl, ne=0;
    const double cc[5]={12.75, 6.941, 1.166, 8.034e-2, 1.996e-3}, n00=n0*1.0e-6, t0divp0=t0/p0;
    //const double norm=1.6022e-13*t0/p0;

    CCi=kb/(mp*g0);
    pih=0.5*pi;

    coschi=cos(zenith[zk][yj]);
    tcos=abs(coschi);

/*--------------------------------------------------*/
/*-----------------photo-ionization ----------------*/
/*--------------------------------------------------*/
    for (m = 0; m < 8; m++) qi[m]=0.0;

    // ionization rates normalized by n0/t0
    Qeuv=0.0;

    Tn=xx[k][j][i].fx[30]*T0;
    ne =uu[k][j][i].fx[6]*n00;  //electron density in cm^{-3}

    /* day time photo-ionization */
    if (zenith[zk][yj] < pih) {
        Nnn[0]=exp(xx[k][j][i].fx[20])*n00;    // O density  in m^-3
        Nnn[1]=exp(xx[k][j][i].fx[23])*n00;    // O2 density
        Nnn[2]=exp(xx[k][j][i].fx[24])*n00;    // N2 density
        Nnn[3]=exp(xx[k][j][i].fx[26])*n00;    // N density
        Nnn[4]=exp(xx[k][j][i].fx[22])*n00;    // He density

/* optical depth due to absorption by O, O2, N2, N, He */
        for (m = 0; m < 5; m++) {
            // scale height of neutrals O, O2, N2, N, He in meter
            Hn[m]=Tn/(ams[m]*gr[xi])*CCi;

            /* calculate Chapman integration along path from Sun to the observation point */
            Xp=rr[i]*r0/Hn[m];
            y=sqrt(0.5*Xp)*tcos;
            if (y > 13.0) y=13.0;

            /* Chapman function in cm^-2 */
            //if (zenith[zk][yj] < pih) {
            Ch[m]=Nnn[m]*(100.0*Hn[m])*sqrt(pih*Xp)*exp(y*y)*cerfc(y);
            /*}
            else {
                dc=Xp*(1.0-tsin);
                if(dc > 100.0) dc=100.0;
                Ch[m]=Nnn[m]*Hn[m]*sqrt(pi2*Xp)
                            *(sqrt(tsin)*exp(dc)-0.5*exp(y*y)*cerfc(y));
            }*/
            //Ch[m]=Nnn[m]*n0*Hn[m]*chapman(Xp, zenith[zk][yj]);
        }

        /* then, determine optical depth and attenuated solar UV 
         * and EUV flux at each of 37 wavelengths */
        for (s = 0; s < 37; s++) {
            /* EuV flux at the location in photons cm^-2 s^-1 */
            tao=( Ch[0]*segabs[s][0]+Ch[1]*segabs[s][1]+Ch[2]*segabs[s][2]
                 +Ch[3]*segabs[s][3]+Ch[4]*segabs[s][4]);
            if (tao > 100.0) tao=100.0;
            euvf=euvflux[s]*exp(-tao);

            /* accumulate ionization coefficients for all wavelengths in s^-1 */
            qi[0] += segion[s][0]*euvf; /* O  + hv --> O+ + e; */
            qi[1] += segion[s][1]*euvf; /* O2 + hv --> O+ + O +e */
            qi[2] += segion[s][2]*euvf; /* O2 + hv --> O2+ + e */
            qi[3] += segion[s][3]*euvf; /* N2 + hv --> N2+ + e */
            qi[4] += segion[s][4]*euvf; /* N2 + hv --> N+ + N + e */
            qi[5] += segion[s][5]*euvf; /* N  + hv --> N+ + e */
            qi[6] += segion[s][6]*euvf; /* He + hv --> He+ + e */

            /* EUV heating of neutrals due to O, O2, N2 
             * absorption of EUV flux in J m^{-3} s^{-1} (normalized by t0/p0), assuming heating efficiency of 0.45 */
            Qeuv += 0.45e6*(segabs[s][0]*Nnn[0]+segabs[s][1]*Nnn[1]+segabs[s][3]*Nnn[2])*euvf*pene[s]*t0divp0;
        }
    }
    else  {
        /* night time photo-ionization rates for O+, O2+, N2+, and NO+ in s^-1 */
        isza = (int)(zenith[zk][yj]*deg) - 90;
        isza = (int)(max((double)isza, 1.));
        for (s = 0; s < 4; s++) {
            qi[0] += segionn[s][0]*fluxn[xi][isza][s];  //O+
            qi[2] += segionn[s][1]*fluxn[xi][isza][s];  //O2+
            qi[3] += segionn[s][2]*fluxn[xi][isza][s];  //N2+
            qi[7] += segionn[s][3]*fluxn[xi][isza][s];  //NO+
        }
    }

    kk[5] =4.5e-10;        // O2+ + NO --> NO+ + O2
    kk[6] =1.2e-10;        // O2+ + N  --> NO+ + O
    kk[10]=3.3e-10;        // N2+ + NO --> NO+ + N2
    kk[12]=1.9e-9;         // H+  + NO --> NO+ + H
    kk[13]=5.2e-10;        // He+ + N2 --> N2+ + He
    kk[14]=9.7e-10;        // He+ + O2 --> O+ + O + He
    kk[15]=1.25e-9;        // He+ + NO --> N+ + He + O
    kk[16]=9.0e-10;        // N+  + NO --> NO+ + N
    kk[17]=2.6e-10;        // N+  + O2 --> NO+ + O
    kk[18]=3.1e-10;        // N+  + O2 --> O2+ + N
    kk[19]=3.7e-11;        // N+  + O2 --> O+ + NO
    kk[20]=2.0e-11;        // N+  + NO --> NO+ + N

    Ti=xx[k][j][i].fx[16]*T0;

    // O+ + N2 --> NO+ + N (note kk[0] not used)
    T1=(0.6363*Ti+0.3637*Tn);
    T300=T1/300.0;
    if (T1 <= 1700.0) kk[1]=1.533e-12 - 5.920e-13*T300 + 8.600e-14*T300*T300;
    else kk[1]=2.730e-12 - 1.155e-12*T300 + 1.483e-13*T300*T300;

    // O+ + O2--> O2+ + O
    T2=(0.667*Ti+0.333*Tn);
    T300=T2/300.0;
    T3002=T300*T300;
    kk[2]= 2.820e-11 - 7.740e-12*T300 + 1.073e-12*T3002-5.170e-14*T3002*T300 + 9.650e-16*T3002*T3002;

    // O+ + NO --> NO+ + O
    T300=Ti/300.0;
    if (Ti <= 1500.0) kk[3]=8.36e-13 - 2.02e-13*T300 + 6.95e-14*T300*T300;
    else {
        T3002=T300*T300;
        kk[3]=5.33e-13 - 1.64e-14*T300 + 4.72e-14*T3002 -7.05e-16*T3002*T300;
    }

    // O+ + H --> H+ + O
    kk[4]=6.0e-10; //2.5e-17*sqrt(Tn)*n0t0;

    Tr=0.5*(Ti+Tn);
    T300=Tr/300.0;
    if (Tr <= 1500.0) {
        // N2+ + O --> O+ + N2
        kk[7]=1.0e-11*pow(300.0/Ti, 0.23);

        // N2+ + O --> NO+ + N
        kk[8]=1.4e-10*pow(T300, -0.44);
    }
    else {
        kk[7]=3.62e-12*pow(Ti/300.0, 0.41);
        kk[8]=5.2e-11*pow(T300, 0.2);
    }

    // N2+ + O2 --> O2+ + N2
    kk[9]=5.0e-11/T300;

    // H+ + O --> O+ + H
    kk[11]=2.2e-11*sqrt(xx[k][j][i].fx[17]);

    // O+ + e --> O + hv
    Te=xx[k][j][i].fx[19]*T0;
    double Te300 = pow(Te/300.0, 0.7);
    kk[21]=4.0e-12*Te300;

    // O2+ + e --> O + O
    if (Te <= 1200.0) kk[22]=1.95e-7*Te300;
    else kk[22]=7.38e-6*pow(1200.0/Te, 0.56);
    //kk[22]=1.6236e-14*pow(1200.0/Te, 0.56)*n0t0;

    // N2+ + e --> N + N
    Te300=300.0/Te;
    kk[23]=1.8e-7*pow(Te300, 0.39);

    // H+ + e --> H
    kk[24]=4.43e-12*pow(Te300, 0.64);

    // He + e --> He
    kk[25]=4.43e-12/pow(Te, 0.7);

    // NO+ + e --> N + O
    kk[26]=4.2e-7*pow(300.0/Te, 0.85);
    
    //NO+ + e --> NO
    kk[27]=6.0e-12;

    //N+ + e --> N
    kk[28]=4.43e-12*pow(Te, 0.7);

    //ion and neutral densities in cm^{-3}
    niO =xx[k][j][i].fx[0]*n00;
    niH =xx[k][j][i].fx[1]*n00;
    niHe=xx[k][j][i].fx[2]*n00;
    niO2=xx[k][j][i].fx[3]*n00;
    niN2=xx[k][j][i].fx[4]*n00;
    niNO=xx[k][j][i].fx[5]*n00;
    niN =xx[k][j][i].fx[6]*n00;

    nO =exp(xx[k][j][i].fx[20])*n00;
    nH =exp(xx[k][j][i].fx[21])*n00;
    nHe=exp(xx[k][j][i].fx[22])*n00;
    nO2=exp(xx[k][j][i].fx[23])*n00;
    nN2=exp(xx[k][j][i].fx[24])*n00;
    nNO=exp(xx[k][j][i].fx[25])*n00;
    nN =exp(xx[k][j][i].fx[26])*n00;

    /* O+ reactions */
    react[1] =kk[1] *niO*nN2;
    react[2] =kk[2] *niO*nO2;
    react[3] =kk[3] *niO*nNO;
    react[4] =kk[4] *niO*nH;

    /* O2+ reactions */
    react[5] =kk[5] *niO2*nNO;
    react[6] =kk[6] *niO2*nN;

    /* N2+ reactions */
    react[7] =kk[7] *niN2*nO;
    react[8] =kk[8] *niN2*nO;
    react[9] =kk[9] *niN2*nO2;
    react[10]=kk[10]*niN2*nNO;

    /* H+ reactions */
    react[11]=kk[11]*niH*nO;
    react[12]=kk[12]*niH*nNO;

    /* He+ reactions */
    react[13]=kk[13]*niHe*nN2;
    react[14]=kk[14]*niHe*nO2;
    react[15]=kk[15]*niHe*nO2;

    /* N+ reactions */
    react[16]=kk[16]*niN*nNO;
    react[17]=kk[17]*niN*nO2;
    react[18]=kk[18]*niN*nO2;
    react[19]=kk[19]*niN*nO2;
    react[20]=kk[20]*niN*nNO;

    /* O+ radiative recombination */
    react[21]=kk[21]*niO*ne;

    /* O2+ combination */
    react[22]=kk[22]*niO2*ne;

    /* N2+ recombination */
    react[23]=kk[23]*niN2*ne;

    /* H+ radiative recombination */
    react[24]=kk[24]*niH*ne;

    /* He radiative recombination */
    react[25]=kk[25]*niHe*ne;

    /* NO+ recombination */
    react[26]=kk[26]*niNO*ne;
    react[27]=kk[27]*niNO*ne;

    /* N+ recombination */
    react[28]=kk[28]*niN*ne;

    /* photo ionization rate in cm^-3 s^-1 */
    photoIon[0]=qi[0]*nO;
    photoIon[1]=qi[1]*nO2;
    photoIon[2]=qi[2]*nO2;
    photoIon[3]=qi[3]*nN2;
    photoIon[4]=qi[4]*nN2;
    photoIon[5]=qi[5]*nN;
    photoIon[6]=qi[6]*nHe;
    photoIon[7]=qi[7]*nNO;

    const double t0divn0=t0/n0;

    /* production rates cm^{-3} s^{-1} & loss coefficients s^{-1} */
    //O+
    Ps[0]=(photoIon[0] + photoIon[1]+ react[7] +react[11] + react[14] + react[19])*t0divn0;
    Ls[0]=kk[1]*nN2 + kk[2]*nO2 + kk[3]*nNO + kk[4]*nH + kk[21]*ne*t0;
    if (Ps[0] > 1.20*Ls[0]*niO || Ls[0]*niO > Ps[0]*1.20) {
        if (Ps[0] > 1.20*Ls[0]*niO) Ps[0]=1.20*Ls[0]*niO; 
        if (Ls[0]*niO > Ps[0]*1.20) Ls[0]=1.20*Ps[0]/niO;
    }

    //H+
    Ps[1]=react[4]*t0divn0;
    Ls[1]=(kk[11]*nO + kk[12]*nNO + kk[24]*ne)*t0;
    if (Ps[1] > 1.20*Ls[1]*niH || Ls[1]*niH > Ps[1]*1.20) {
        if (Ps[1] > 1.20*Ls[1]*niH) Ps[1]=1.20*Ls[1]*niH; 
        if (Ls[1]*niH > Ps[1]*1.20) Ls[1]=1.20*Ps[1]/niH;
    }

    //He+
    Ps[2]=photoIon[6]*t0divn0;
    Ls[2]=(kk[13]*nN2 + kk[14]*nO2 + kk[15]*nNO + kk[25]*ne)*t0;
    if (Ps[2] > 1.20*Ls[2]*niHe || Ls[2]*niHe > Ps[2]*1.20) {
        if (Ps[2] > 1.20*Ls[2]*niHe) Ps[2]=1.20*Ls[2]*niHe; 
        if (Ls[2]*niHe > Ps[2]*1.20) Ls[2]=1.20*Ps[2]/niHe;
    }

    //O2+
    Ps[3]=(photoIon[2] + react[2] + react[9] + react[18])*t0divn0;
    Ls[3]=(kk[5]*nNO + kk[6]*nN + kk[22]*ne)*t0;
    if (Ps[3] > 1.20*Ls[3]*niO2 || Ls[3]*niO2 > Ps[3]*1.20) {
        if (Ps[3] > 1.20*Ls[3]*niO2) Ps[3]=1.20*Ls[3]*niO2; 
        if (Ls[3]*niO2 > Ps[3]*1.20) Ls[3]=1.20*Ps[3]/niO2;
    }

    //N2+
    Ps[4]=(photoIon[3] + react[13])*t0divn0;
    Ls[4]=((kk[7] + kk[8])*nO + kk[9]*nO2 + kk[10]*nNO + kk[23]*ne)*t0;
    if (Ps[4] > 1.20*Ls[4]*niN2 || Ls[4]*niN2 > Ps[4]*1.20) {
        if (Ps[4] > 1.20*Ls[4]*niN2) Ps[4]=1.20*Ls[4]*niN2; 
        if (Ls[4]*niN2 > Ps[4]*1.20) Ls[4]=1.20*Ps[4]/niN2;
    }

    //NO+
    Ps[5]=( react[1] + react[3] + react[5] + react[6] + react[8] + react[10] 
                       +react[12] +react[16] + react[20])*t0divn0;
    Ls[5]=(kk[26] + kk[26])*ne*t0;
    if (Ps[5] > 1.20*Ls[5]*niNO || Ls[5]*niNO > Ps[5]*1.20) {
        if (Ps[5] > 1.20*Ls[5]*niNO) Ps[5]=1.20*Ls[5]*niNO; 
        if (Ls[5]*niNO > Ps[5]*1.20) Ls[5]=1.20*Ps[5]/niNO;
    }

    Ps[6]= react[15]*t0divn0;
    Ls[6]= (kk[16]*nNO +(kk[17]+kk[18]+kk[19])*nO2 + kk[20]*nNO + kk[28]*ne)*t0;
    if (Ps[6] > 1.20*Ls[6]*niN || Ls[6]*niN > Ps[6]*1.20) {
        if (Ps[6] > 1.20*Ls[6]*niN) Ps[6]=1.20*Ls[6]*niN; 
        if (Ls[6]*niN > Ps[6]*1.20) Ls[6]=1.20*Ps[6]/niN;
    }

    //O
    Ps[7]=( photoIon[1] + react[2] + react[3] + react[4]+react[6] + react[14]
                       +react[15] + react[17] + react[21] + react[22] + react[26])*t0divn0;
    Ls[7]=(qi[0] + (kk[7] + kk[8])*niN2 + kk[11]*niH)*t0;
    if (Ps[7] > 1.20*Ls[7]*nO || Ls[7]*nO > Ps[7]*1.20) {
        if (Ps[7] > 1.20*Ls[7]*nO) Ps[7]=1.20*Ls[7]*nO; 
        if (Ls[7]*nO > Ps[7]*1.20) Ls[7]=1.20*Ps[7]/nO;
    }

    //H
    Ps[8]=(react[11] + react[12] + react[24])*t0divn0;
    Ls[8]=kk[4]*niO*t0;
    if (Ps[8] > 1.20*Ls[8]*nH || Ls[8]*nH > Ps[8]*1.20) {
        if (Ps[8] > 1.20*Ls[8]*nH) Ps[8]=1.20*Ls[8]*nH; 
        if (Ls[8]*nH > Ps[8]*1.20) Ls[8]=1.20*Ps[8]/nH;
    }

    //He
    Ps[9]=(react[13] + react[14] + react[15] + react[25])*t0divn0;
    Ls[9]=qi[6]*t0;
    if (Ps[9] > 1.20*Ls[9]*nHe || Ls[9]*nHe > Ps[9]*1.20) {
        if (Ps[9] > 1.20*Ls[9]*nHe) Ps[9]=1.20*Ls[9]*nHe; 
        if (Ls[9]*nHe > Ps[9]*1.20) Ls[9]=1.20*Ps[9]/nHe;
    }

    //O2
    Ps[10]=react[5]*t0divn0;
    Ls[10]= ((qi[1] + qi[2]) + kk[2]*niO + kk[9]*niN2 +kk[14]*niHe +(kk[17]+kk[18]+kk[19])*niN)*t0;
    if (Ps[10] > 1.20*Ls[10]*nO2 || Ls[10]*nO2 > Ps[10]*1.20) {
        if (Ps[10] > 1.20*Ls[10]*nO2) Ps[10]=1.20*Ls[10]*nO2; 
        if (Ls[10]*nO2 > Ps[10]*1.20) Ls[10]=1.20*Ps[10]/nO2;
    }

    //N2
    Ps[11]=(react[7] + react[9] + react[10])*t0divn0;
    Ls[11]=((qi[3]+qi[4]) + kk[1]*niO + kk[13]*niHe)*t0;
    if (Ps[11] > 1.20*Ls[11]*nN2 || Ls[11]*nN2 > Ps[11]*1.20) {
        if (Ps[11] > 1.20*Ls[11]*nN2) Ps[11]=1.20*Ls[11]*nN2; 
        if (Ls[11]*nN2 > Ps[11]*1.20) Ls[11]=1.20*Ps[11]/nN2;
    }

    //NO
    Ps[12]=(react[19] + react[27])*t0divn0;
    Ls[12]=( qi[6] + kk[3]*niO + kk[5]*niO2 + kk[10]*niN2 + kk[12]*niH + kk[15]*niHe
                        +(kk[16]+kk[20])*niN)*t0;
    if (Ps[12] > 1.20*Ls[12]*nNO || Ls[12]*nNO > Ps[12]*1.20) {
        if (Ps[12] > 1.20*Ls[12]*nNO) Ps[12]=1.20*Ls[12]*nNO; 
        if (Ls[12]*nNO > Ps[12]*1.20) Ls[12]=1.20*Ps[12]/nNO;
    }

    //N
    Ps[13]=( photoIon[4] + react[8] + 2.0*react[23] + react[26] + react[16] + react[20]
                        +react[18] + react[28])*t0divn0;
    Ls[13]=(qi[5] + kk[6]*niO2)*t0;
    if (Ps[13] > 1.20*Ls[13]*nN || Ls[13]*nN > Ps[13]*1.20) {
        if (Ps[13] > 1.20*Ls[13]*nN) Ps[13]=1.20*Ls[13]*nN; 
        if (Ls[13]*nN > Ps[13]*1.20) Ls[13]=1.20*Ps[13]/nN;
    }

/*-----------------------------------------------------------------------------*/
/*--------------- electron heating rate ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
    /* local photoelectron heating Swartz & Nisbet, JGR, 1972 in J m^-3 s^-1 normalzed by t0/p0 */
    if (zh[xi] < 300.0) {                   
        R=log(ne/(0.1*nO + nO2 + nN2));
        R2=R*R;
        epsl=exp(-(cc[0]+cc[1]*R+cc[2]*R2+cc[3]*R2*R+cc[4]*R2*R2)); //man energy per photon (eV)

        Qphoto=q*epsl*( photoIon[0]+photoIon[1]+photoIon[2]+photoIon[3]
                       +photoIon[4]+photoIon[5]+photoIon[6]+photoIon[7])*t0divp0;
    }
    //else {
    /* P. Richards et al., JGR, [1984] */
    //    Qee=0.35*pesimp(zh[xi],zenith[zk][yj],user->f107,Te,Tn,ne,nO,nO2,nN2)*norm;
    //}
}
