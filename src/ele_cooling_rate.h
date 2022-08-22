/*************************************************************************
 * ele_cooling_rate.cpp  
 * calculate electron cooling rate by excitation of neutral states through
 * inelastic electron - neutral collisions. 
 * See Schun and Nagy, Geophy. Review, 113(8), A08307, 2008
 * 
   Input:  parms  - system parameters
              xx  - solution vector

    Output: non
 *  
 * Jiannan Tu
 * 5/21/2022
*************************************************************************/
#include <ctime>

#include "reconstruction.h"

inline double ele_cooling_rate(Field ***xx, double ne, double Te, double Tn, int i, int j, int k)
{
    int    m;
    double fc, gc, he, Z, Dxi[3], Exi[3], dd;
    double nO, nO2, nN2, Le[4], temp1; //, temp2;

    const double ei[3]={0.02, 0.028, 0.008};
    //const double Ei[3]={228.0, 326.0, 98.0};
    const double Ai[3]={7.833e-6, 9.466e-6, 1.037e-8};
    const double Bi[3]={1.021, 0.8458, 1.633};
    const double Ci[3]={1.009, 0.9444, 1.466};

    // O, O2, N2, H density in cm^-3
    nO =xx[k][j][i].fx[12];
    nO2=xx[k][j][i].fx[15];
    nN2=xx[k][j][i].fx[16];

    /* due to impact excitation of */
    /* due to N2 vibration */
    fc=1.06e4+7.51e3*tanh(1.10e-3*(Te-1800.0));
    gc=3300.0+(1.233-2.056e-4*(Te-4000.0))*(Te-1000.0);

    Le[0]=2.99e-12*nN2*exp(fc*(Te-2000.0)/(2000.0*Te))*(exp(gc*(Tn-Te)/(Te*Tn))-1.0);

    /* due to O2 vibration */
    he=3300.0-839.0*sin(1.91e-4*(Te-2700.0));
    Le[1]=5.196e-13*nO2*exp(he*(Te-700.0)/(700.0*Te))*(exp(2770.0*(Tn-Te)/(Te*Tn))-1.0);

    /* due to O fine structure: constant part of cooling rate */
    Z=5.0+3.0*exp(-228.0/Tn)+exp(-326.0/Tn);
    Dxi[0]=exp(-228.0/Tn);
    Dxi[1]=exp(-326.0/Tn);
    Dxi[2]=exp(-326.0/Tn);

    Exi[0]=exp(-228.0/Te);
    Exi[1]=exp(-326.0/Te);
    Exi[2]=exp(-(98.0/Te+228.0/Tn));

    Le[2]=0.0;
    for (m = 0; m < 3; m++) {
        //original expression regards cooling as negative heating
        temp1=Ai[m]*Ci[m]*pow(Te, (Bi[m]-0.5));
        //temp2=((1.0+Bi[m])*Dxi[m]+(Ei[m]/Te+1.0+Bi[m])*Exi[m]);

        Le[2] += temp1*ei[m]*(Exi[m]-Dxi[m]); // + 5.91e-9*(Tn-Te)*temp2);
    }
    Le[2]=-fabs(Le[2]*8.629e-6*nO/Z);

    /* due to O(1D) excitation */
    dd=2.4e4+(0.3-1.947e-5*(Te-4000.0))*(Te-1500.0);
    Le[3]=1.57e-12*nO*exp(dd*(Te-3000.0)/(3000.0*Te))*(exp(22713.0*(Tn-Te)/(Te*Tn))-1.0);

    //cooling rate in J cm^{-3} s^{-1}
    return (Le[0]+Le[1]+Le[2]+Le[3])*ne*e;
}
