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

inline void ele_cooling_rate(Field ***xx, double ne, double Te, double Tn, int i, int j, int k, 
    double &Ce, double &dCedTe, double &dCedTn)
{
    int    m;
    double Te12, Te2, fc, gc, he, Z, Dxi[3], Exi[3], dd;
    double nO, nO2, nN2, nH, Tn2, Le[11], temp0, temp, temp1, temp2;
    double ke, ke2, ki2, lnk, dgdTe;

    const double ei[3]={0.02, 0.028, 0.008};
    const double Ei[3]={228.0, 326.0, 98.0};
    const double Ai[3]={8.58e-6, 7.201e-6, 2.463e-7};
    const double Bi[3]={1.008, 0.9617, 1.1448};
    const double Ci[3]={1.009, 0.9444, 1.466};
    const double e2_div_kb=e*(e/kb);
    const double T00=22755.884, T1=48414.304;

    double Ti;

    if (Te > Tn) {
        Te12=sqrt(Te);
        Te2=Te*Te;
        Tn2=Tn*Tn;

        // O, O2, N2, H density in cm^-3
        nO =xx[k][j][i].fx[12];
        nO2=xx[k][j][i].fx[15];
        nN2=xx[k][j][i].fx[16];
        nH =xx[k][j][i].fx[13];

    /* due to impact excitation of
       * N2 rotation in eV s^-1 (not multiplied by ne yet) */
        temp=2.9e-14*nN2/Te12;
        Le[0]=temp*(Tn-Te);

        dCedTe=-0.5*temp*(Tn/Te+1.0);
        dCedTn=temp;

    /* due to N2 vibration */
        fc=1.06e4+7.51e3*tanh(1.10e-3*(Te-1800.0));
        gc=3300.0+(1.233-2.056e-4*(Te-4000.0))*(Te-1000.0);

        temp=2.99e-12*nN2*exp(fc*(Te-2000.0)/(2000.0*Te));
        temp1=exp(gc*(Tn-Te)/(Te*Tn));
        Le[1]=temp*(temp1-1.0);

        temp2=temp*temp1;
        dgdTe=1.233-2.056e-4*(2.0*Te-5000.0);
        dCedTe += Le[1]*(fc/Te2+(0.004105-8.261/Te)*cosh(1.10e-3*(Te-1800.0)))
                 +temp2*(dgdTe-gc/Te2);
        dCedTn += temp2*gc/Tn2;

    //cooling due to elastic collisions with N2
        Le[2]=1.77e-19*nN2*(1.0-1.21e-4*Te)*Te*(Tn-Te); //only treated explicitly

       /* due to O2 rotation */
        temp=6.9e-14*nO2/Te12;
        Le[3]=temp*(Tn-Te);

        dCedTe += -0.5*temp*(Tn/Te+1.0);
        dCedTn += temp;

       /* due to O2 vibration */
        temp1=1.91e-4*(Te-2700.0);
        he=3300.0-839.0*sin(temp1);

        temp=5.196e-13*nO2*exp(he*(Te-700.0)/(700.0*Te));
        temp2=exp(2770.0*(Tn-Te)/(Te*Tn));
        Le[4]=temp*(temp2-1.0);

        temp=temp*temp2;
        dCedTe += Le[4]*(he/Te2-(2.2892714286e-4-0.160249/Te)*cos(temp1))
                 -temp*2770.0/Te2;
        dCedTn += temp*2770/Tn2;

        //cooling due to elastic collisions with O2
        Le[5]=1.21e-18*nO2*(1.0+3.6e-2*Te12)*Te12*(Tn-Te);

        /* due to O fine structure */
        Z=5.0+3.0*exp(-228.0/T1)+exp(-326.0/T00);
        Dxi[0]=exp(-228.0/T1);
        Dxi[1]=exp(-326.0/T00);
        Dxi[2]=exp(-326.0/T00);

        Exi[0]=exp(-228.0/Te);
        Exi[1]=exp(-326.0/Te);
        Exi[2]=exp(-(98.0/Te+228.0/Tn));

        Le[6]=0.0; temp0=0.0, temp=0.0;
        for (m = 0; m < 3; m++) {
            //original expression regards cooling as negative heating
            temp1=Ai[m]*Ci[m]*pow(Te, (Bi[m]-0.5));
            temp2=((1.0+Bi[m])*Dxi[m]+(Ei[m]/Te+1.0+Bi[m])*Exi[m]);

            Le[6] += temp1*(ei[m]*(Exi[m]-Dxi[m]) + 5.91e-9*(Te-Tn)*temp2);

            temp0 += temp1*( (Bi[m]-0.5)*(ei[m]*(Exi[m]-Dxi[m]) + 5.91e-9*(Te-Tn)*temp2)
                            +5.91e-6*( (Tn-Te)*(Ei[m]/Te+Bi[m])*Exi[m]*Ei[m]/Te2
                                      -temp2-ei[m]*Exi[m]*Ei[m]/Te2));
            temp += temp1*temp2;
        }
        Le[6]=Le[6]*8.629e-6*nO/Z;

        dCedTe += 8.629e-6*nO/Z*temp0;
        dCedTn += 5.136972e-14*temp;

    /* due to O(1D) excitation */
        dd=2.4e4+(0.3-1.947e-5*(Te-4000.0))*(Te-1500.0);

        temp=1.57e-12*nO*exp(dd*(Te-3000.0)/(3000.0*Te));
        temp1=exp(22713.0*(Tn-Te)/(Te*Tn));
        Le[7]=temp*(temp1-1.0);

        temp2=temp*temp1;
        dCedTe += Le[7]*(dd/Te2+(1.0/3000.0-1.0/Te)*(0.3-1.94e-5*(2.0*Te-4500.0)))
                 -temp*(22713.0/Te2);
        dCedTn += temp*(22713.0/Tn2);

    //cooling due to elastic collisions with O
        Le[8]=7.9e-19*nO*(1.0+5.7e-4*Te)*Te12*(Tn-Te);

    //cooling due to elastic collisions with H
        Le[9]=9.63e-16*nH*(1.0-1.35e-4*Te)*Te12*(Tn-Te);

    /* cooling due to elastic collisions with positive ions */
        Ti = xx[k][j][i].fx[10]/(ne*kb);
        if (Te > Ti) {
            ki2=4.0*pi*ne*e2_div_kb/Ti;
            ke2=4.0*pi*ne*e2_div_kb/Te;
            ke=sqrt(ke2);
            lnk= log(4.0*Te/(e2_div_kb*ke))-1.154-(ki2+ke2)/ki2*log(sqrt(ki2+ke2)/ke);

            Le[10]=3.2e-8*(Te-Ti)/(Te12*Te)*lnk
                         *( xx[k][j][i].fx[0]+4.0*xx[k][j][i].fx[2]+16.0*xx[k][j][i].fx[1]
                           +0.5*xx[k][j][i].fx[3]+0.53*xx[k][j][i].fx[5]);
        }
        else Le[10]=0.0;

        //electron cooling rate in Joule cm^{-3} s^{-1}
        Ce=(Le[0]+Le[1]+Le[2]+Le[3]+Le[4]+Le[5]+Le[6]+Le[7]+Le[8]+Le[9]+Le[10])*ne*e;

        dCedTe = dCedTe*ne*e;
        dCedTn = dCedTn*ne*e;
    }
    else {
        Ce=0.0; dCedTe=0.0; dCedTn=0.0;
    }

    return;
}
