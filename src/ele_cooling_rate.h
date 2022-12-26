/*************************************************************************
 * ele_cooling_rate.cpp  
 * calculate electron cooling rate by excitation of neutral states through
 * inelastic electron - neutral collisions. 
 * See Schunk and Nagy, Geophy. Review, 113(8), A08307, 2008
 * 
   Input:  parms  - system parameters
              xx  - solution vector

    Output: non
 *  
 * Jiannan Tu
 * 5/21/2022
*************************************************************************/
#include "param.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

inline double ele_cooling_rate(Field ***xx, double Te, double Tn, double ne, int i, int j, int k)
{
    int    m;
    double Td, fc, gc, he, Z, Dxi[3], Exi[3], dd;
    double nO, nO2, nN2, nH, Le[10], Cee=0.0;
    const double n00=n0*1.0e-6, n0t0_p0=n0*t0/p0;

    const double ei[3]={0.02, 0.028, 0.008};
    const double Ei[3]={228.0, 326.0, 98.0};
    const double Ai[3]={8.58e-6, 7.201e-6, 2.463e-7};
    const double Bi[3]={1.008, 0.9617, 1.1448};
    const double Ci[3]={1.009, 0.9444, 1.466};

    //electron density in m^-3, electron and neutral temperature in K
    Te=Te*T0; Tn=Tn*T0;

    Td=sqrt(Te);

    // density in cm^-3
    nO =exp(xx[k][j][i].fx[20])*n00;
    nO2=exp(xx[k][j][i].fx[23])*n00;
    nN2=exp(xx[k][j][i].fx[24])*n00;
    nH =exp(xx[k][j][i].fx[21])*n00;

    /* rate (divided by ne) of cooling due to impact excitation of
    * N2 rotation in eV s^-1 (not multiplied by ne yet) */
    Le[0]=2.9e-14*nN2*(Te-Tn)/Td;

    /* due to N2 vibration */
    fc=1.06e4+7.51e3*tanh(1.10e-3*(Te-1800.0));
    gc=3300.0+(1.233-2.056e-4*(Te-4000.0))*(Te-1000.0);

    //original expression regards cooling as negative heating
    Le[1]=2.99e-12*nN2*exp(fc*((Te-2000.0)/(2000.0*Te)))*(1.0-exp(gc*(Tn-Te)/(Te*Tn)));

    //cooling due to elastic collisions with N2
    Le[2]=1.77e-19*nN2*(1.0-1.21e-4*Te)*Te*(Te-Tn);

    /* due to O2 rotation */
    Le[3]=6.9e-14*nO2*(Te-Tn)/Td;

    /* due to O2 vibration */
    he=3300.0-839.0*sin(1.91e-4*(Te-2700.0));

    //original expression regards cooling as negative heating
    Le[4]=5.196e-13*nO2*exp(he*(Te-700.0)/(700.0*Te))*(1.0-exp(2770.0*(Tn-Te)/(Te*Tn)));

    //cooling due to elastic collisions with O2
    Le[5]=1.21e-18*nO2*(1.0+3.6e-2*Td)*Td*(Te-Tn);

    /* due to O fine structure */
    Z=5.0+3.0*exp(-228.0/Tn)+exp(-326.0/Tn);
    Dxi[0]=exp(-228.0/Tn);
    Dxi[1]=exp(-326.0/Tn);
    Dxi[2]=exp(-326.0/Tn);

    Exi[0]=exp(-228.0/Te);
    Exi[1]=exp(-326.0/Te);
    Exi[2]=exp(-(98.0/Te+228.0/Tn));

    Le[6]=0.0;
    for (m = 0; m < 3; m++) {
        //original expression regards cooling as negative heating
        Le[6] += Ai[m]*Ci[m]*pow(Te, (Bi[m]-0.5))
                      *(ei[m]*(Exi[m]-Dxi[m])+5.91e-9*(Te-Tn)*((1.0+Bi[m])*Dxi[m]+(Ei[m]/Te+1.0+Bi[m])*Exi[m]));
    }
    Le[6]=Le[6]*8.629e-6*nO/Z;

    /* due to O(1D) excitation */
    dd=2.4e4+(0.3-1.947e-5*(Te-4000.0))*(Te-1500.0);

    //original expression regards cooling as negative heating
    Le[7]=1.57e-12*nO*exp(dd*(Te-3000.0)/(3000.0*Te))*(1.0-exp(22713.0*(Tn-Te)/(Te*Tn)));

    //cooling due to elastic collisions with O
    Le[8]=7.9e-19*nO*(1.0+5.7e-4*Te)*Td*(Te-Tn);

    //cooling due to elastic collisions with H
    Le[9]=9.63e-16*nH*(1.0-1.35e-4*Te)*Td*(Te-Tn);

    //normalized electron cooling rate
    Cee=(Le[0]+Le[1]+Le[2]+Le[3]+Le[4]+Le[5]+Le[6]+Le[7]+Le[8]+Le[9])*ne*q*n0t0_p0;

    return Cee;
}
