/*************************************************************************
 * neu_cooling_rate.cpp  
 *    calculate neutral cooling rates due to infrared radiative loss by
 * O at wavelength of 63 micro-meter and 147 micro-meter, NO at 5.3 micro-meter  
   Input:  parms  - system parameters
              xx  - solution vector

    Output: none
 *  
 * Jiannan Tu
 * 5/21/2022
*************************************************************************/
/*#include <iostream>
#include <fstream>
#include <iomanip>*/

#include <cmath>
#include "funcdef.h"

using namespace std;

inline double neu_cooling_rate(Field ***xx, int i, int j, int k)
{
    double Tn, nO, nNO, Tx1, Tx2, Cn; //chi, , E21, E22;
    const double t0divp0=t0/p0;

    //double *tao = new double[xm];

    /* first loop i is to calculate optical depth on all grids */
    //for (i = xs+xm-2; i >= xs; i--) {
        /* density in cm^-3 & temperature in K (not using arrays for nO,  
        * nNO, Tn because of better efficiency when using openMP */
        //nO =xx[k][j][i].fx[34]*n00;
        //nNO=xx[k][j][i].fx[39]*n00;
        //Tn =xx[k][j][i].fx[48]*T0;

        //Tx1=exp(228.0/Tn);

        //calculate optical depth
        //tao[i-xs]=tao[i-xs+1]+0.5e-14*nO*(Tx1-1.0)*(rr[i+1]-rr[i])*r0
        //                             /(sqrt(Tn)*(0.6+0.2*exp(-98.0/Tn)+Tx1));
    //}

    //O and NO density in m^{-3}
    nO =exp(xx[k][j][i].fx[20])*n0;
    nNO=exp(xx[k][j][i].fx[25])*n0;
    Tn =xx[k][j][i].fx[30]*T0;

    Tx1=exp(-228.0/Tn);
    Tx2=exp(-326.0/Tn);

    //calculate chi for reduction factor
    //E21=expon_integral(tao[i-xs], 60);
    //E22=expon_integral(tao[0]-tao[i-xs],60);
    //chi=0.5*(Tx1-1.0)*((2.0-E21-E22)/(Tx1-1.0)+E22/(exp(228.0/TOb)-1.0));

    /* cooling rate in Joule m^-3 s^-1 (normalized) */
    Cn=nO*( (1.69e-25*Tx1+4.59e-27*Tx2)/(1.0+0.6*Tx1+0.2*Tx2)
           +3.24015e-29*nNO*exp(-2714.57/Tn)/(6.5e-11*nO+13.3))*t0divp0;

    //delete[] tao;

    return Cn;
}

