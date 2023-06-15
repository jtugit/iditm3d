/************************************************************************** 
 * const_normalize
 *   evaluate certain constants and perform normalization
 *
 * Jiannan Tu
 * 5/21/2022
 **************************************************************************/
#include <cmath>
#include "param.h"

using namespace std;

void const_normalize(AppCtx *params)
{
    int   i;

/* evaluate normalization parameters from 4 basic parameters: r0, n0, B0, & mp */
    v0= B0/sqrt(mu0*n0*mp);
    t0= r0/v0;
    E0= v0*B0;
    g0= v0/t0;
    p0= B0*B0/mu0;
    j00= B0/(mu0*r0);
    e0= j00/(n0*v0);
    T0= p0/(kb*n0);
    q0= p0*v0;
    lamda0= q0*r0/ t0;
    beta0= q0/j00;

    //normalized elementary charge
    e=q/e0;

    //normalized time steps
    dt=dt/t0; dt_half=0.5*dt; dt2=dt/params->alpha;

    //normalized gravitational acceleration on the surface of the earth
    gen=ge/g0;

    //normalized earth's rotational frequency
    w0n=w0*t0;

    //normalized charge/ion mass ratio for various ion species
    for (i = 0; i < sl; i++) qms[i]=e/ams[i];

    //normalized (multiplied by t0) electron gyro-frequency using normalization factor B0
    //Omegae=q*B0*t0/me;
    
/* coefficients of all collision frequencies, multiplied by t0 */
    coe[0]=54.5*t0;                  //e - ions
    coe[1]=8.90e-11*t0;              //e - O
    coe[2]=4.50e-9*t0;               //e - H
    coe[3]=4.60e-10*t0;              //e - He
    coe[4]=1.82e-10*t0;              //e - O2
    coe[5]=2.33e-11*t0;              //e - N2
    coe[6]=1.0e-11*t0;               //e - NO
    coe[7]=1.0e-11*t0;               //e - N

    coiO[0]=7.70e-2*t0;        //O+ - H+
    coiO[1]=1.40e-1*t0;        //O+ - He+
    coiO[2]=2.60e-1*t0;        //O+ - O2+
    coiO[3]=2.50e-1*t0;        //O+ - N2+
    coiO[4]=2.60e-1*t0;        //O+ - NO+
    coiO[5]=0.22*t0;           //O+ - N+
    coiO[6]=1.5*3.67e-11*t0;         //O+ - O
    coiO[7]=6.61e-11*t0;             //O+ - H
    coiO[8]=1.32e-10*t0;             //O+ - He
    coiO[9]=6.64e-10*t0;             //O+ - O2
    coiO[10]=6.82e-10*t0;            //O+ - N2
    coiO[11]=2.69e-11*t0;            //O+ - NO
    coiO[12]=4.62e-10*t0;            //O+ - N

    coiH[0]=1.23*t0;        //H+ - O+
    coiH[1]=1.14*t0;        //H+ - He+
    coiH[2]=1.25*t0;        //H+ - O2+
    coiH[3]=1.25*t0;        //H+ - N2+
    coiH[4]=1.25*t0;        //H+ - NO+
    coiH[5]=1.23*t0;        //H+ - N+
    coiH[6]=6.61e-11*t0;          //H+ - O
    coiH[7]=2.65e-10*t0;          //H+ - H
    coiH[8]=10.6e-10*t0;          //H+ - He
    coiH[9]=3.20e-9*t0;           //H+ - O2
    coiH[10]=3.36e-9*t0;          //H+ - N2
    coiH[11]=2.69e-11*t0;         //H+ - NO
    coiH[12]=26.1e-10*t0;         //H+ - N

    coiHe[0]=0.57*t0;        //He+ - O+
    coiHe[1]=0.28*t0;        //He+ - H+
    coiHe[2]=0.60*t0;        //He+ - O2+
    coiHe[3]=0.59*t0;        //He+ - N2+
    coiHe[4]=0.60*t0;        //He+ - NO+
    coiHe[5]=0.56*t0;        //He+ - N+
    coiHe[6]=10.1e-10*t0;       //He+ - O
    coiHe[7]=4.71e-10*t0;       //He+ - H
    coiHe[8]=8.73e-11*t0;       //He+ - He
    coiHe[9]=15.3e-10*t0;       //He+ - O2
    coiHe[10]=16.0e-10*t0;      //He+ - N2
    coiHe[11]=2.69e-11*t0;      //He+ - NO
    coiHe[12]=11.9e-10*t0;      //He+ - N

    con[0]=2.26216e-11*t0;    //O - O2
    con[1]=7.60616e-12*t0;    //O - N2
    con[2]=5.15410e-12*t0;    //O2 - N2
    con[3]=9.05133e-12*t0;    //H - O
    con[4]=5.43880e-12*t0;    //H - O2
    con[5]=6.05369e-12*t0;    //H - N2
    con[6]=2.33451e-11*t0;    //H - He
    con[7]=1.49978e-11*t0;    //He - O
    con[8]=8.63023e-12*t0;    //He - O2
    con[9]=1.00277e-11*t0;    //He - N2

    coiO2[0]=1.30e-1*t0;       //O2+ - O+
    coiO2[1]=3.90e-2*t0;       //O2+ - H+
    coiO2[2]=7.50e-2*t0;       //O2+ - He+
    coiO2[3]=1.50e-1*t0;       //O2+ - N2+
    coiO2[4]=1.60e-1*t0;       //O2+ - NO+
    coiO2[5]=0.12*t0;          //O2+ - N+
    coiO2[6]=2.31e-10*t0;        //O2+ - O
    coiO2[7]=6.50e-11*t0;        //O2+ - H
    coiO2[8]=7.00e-11*t0;        //O2+ - He
    coiO2[9]=2.59e-11*t0;        //O2+ - O2
    coiO2[10]=4.13e-10*t0;       //O2+ - N2
    coiO2[11]=2.69e-11*t0;       //O2+ - NO
    coiO2[12]=2.64e-10*t0;       //O2+ - N
    
    coiN2[0]=1.50e-1*t0;       //N2+ - O+
    coiN2[1]=4.50e-2*t0;       //N2+ - H+
    coiN2[2]=8.50e-2*t0;       //N2+ - He+
    coiN2[3]=1.80e-1*t0;       //N2+ - O2+
    coiN2[4]=1.70e-1*t0;       //N2+ - NO+
    coiN2[5]=0.14*t0;          //N2+ - N+
    coiN2[6]=2.58e-10*t0;         //N2+ - O
    coiN2[7]=7.40e-11*t0;         //N2+ - H
    coiN2[8]=7.90e-11*t0;         //N2+ - He
    coiN2[9]=4.49e-10*t0;         //N2+ - O2
    coiN2[10]=5.14e-11*t0;        //N2+ - N2
    coiN2[11]=2.69e-11*t0;        //N2+ - NO
    coiN2[12]=2.95e-10*t0;        //N2+ - N
 
    coiNO[0]=0.14*t0;          //NO+ - O+
    coiNO[1]=4.20e-2*t0;       //NO+ - H+
    coiNO[2]=8.00e-2*t0;       //NO+ - He+
    coiNO[3]=1.70e-1*t0;       //NO+ - O2+
    coiNO[4]=1.60e-1*t0;       //NO+ - N2+
    coiNO[5]=0.13*t0;          //NO+ - N+
    coiNO[6]=2.44e-10*t0;        //NO+ - O
    coiNO[7]=6.90e-11*t0;        //NO+ - H
    coiNO[8]=7.40e-11*t0;        //NO+ - He
    coiNO[9]=4.27e-10*t0;        //NO+ - O2
    coiNO[10]=4.34e-10*t0;       //NO+ - N2
    coiNO[11]=2.69e-11*t0;       //NO+ - NO
    coiNO[12]=2.79e-10*t0;       //NO+ - N
 
    coiN[0]=0.25*t0;          //N+ - O+
    coiN[1]=0.088*t0;         //N+ - H+
    coiN[2]=0.16*t0;          //N+ - He+
    coiN[3]=0.28*t0;          //N+ - O2+
    coiN[4]=0.28*t0;          //N+ - N2+
    coiN[5]=0.28*t0;          //N+ - NO+
    coiN[6]=4.42e-10*t0;        //N+ - O
    coiN[7]=1.45e-10*t0;        //N+ - H
    coiN[8]=1.49e-10*t0;        //N+ - He
    coiN[9]=7.25e-10*t0;        //N+ - O2
    coiN[10]=7.47e-10*t0;       //N+ - N2
    coiN[11]=2.69e-11*t0;       //N+ - NO
    coiN[12]=3.83e-11*t0;       //N+ - N

    return;
}
