/*****************************************************************************
 *  input_prev
 *   input parameters & simulation results from a previous simulation. 
 *
 *   Jiannan Tu
 *   12/11/2013, 1/29/2014
 ****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cmath>

using namespace std;

#include "param.h"
#include "funcdef.h"

void smooth(Field ***xx, int xs, int xm, int j, int k, int s, int ncomp);
void shapiro(double *x, int n, double cshap);

int input_iri_msis(DM da, Vec X, Field ***xx, AppCtx *params)
{
    Vec      localX;
    Field    ***localxx;
    PetscInt i, j, k;
    int      s;
    double   f20[14], ne, Nn;
    fstream  irifstr, msisfstr;
    PetscInt xs, ys, zs, xm, ym, zm;
    int      i0[7]={0,0,0,0,0,0,0}, im[7]={0,0,0,0,0,0,0};

    //yy=to_string(params->iyr);

    string fname1 = params->workdir + "/inp/" + params->irifln;
    irifstr.open(fname1, fstream::in);
    if(!irifstr) {
        cout << "Can't open file " << fname1 << endl;
        return -1;
    }

    string fname2 = params->workdir + "/inp/" + params->msisfln;
    msisfstr.open(fname2, fstream::in);
    if(!msisfstr) {
        cout << "Can't open file " << fname2 << endl;
        return -1;
    }

    /* skip headers in iri and msis data files */
    for (s = 0; s < 31; s++) irifstr.ignore(200, '\n');
    for (s = 0; s < 7; s++) msisfstr.ignore(200, '\n');

    DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

    //set up local arrays
    DMGetLocalVector(da, &localX);
    DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);
    DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);
    DMDAVecGetArray(da, localX, &localxx);

    double y1, y2;
    int    ims, numHi, numHei;

    for (k = zs; k < zs+zm; k++) {

        for (j = ys; j < ys+ym; j++) {

            for (s = 0; s < sl; s++) {
                i0[s]=0.0; im[s]=0.0;
            }

            numHi=0; numHei=0;
            for (i = xs; i < xs+xm; i++) {
                /* read IRI ion & electron temperature, ion density */
                for (s = 0; s < 14; s++) {
                    irifstr >> f20[s];
                    if (s >=5 && f20[s] < 0.0) {
                        cout<<"Negative data from IRI file at (i, j, k, l) = ";
                        cout<<"("<<i<<", "<<j<<", "<<k<<", "<<s<<")"<<endl;
                    }
                }
                irifstr.ignore(1,'\n');

                // O+, H+, He+, O2+, N2+, NO+, N+ density (cm^{-3}) at (rC, thetaC, phi_k).
                // Density from IRI in m^-3 
                localxx[k][j][i].fx[0]=f20[8]*1.0e-6;             //O+
                localxx[k][j][i].fx[1]=f20[9]*1.0e-6;             //H+
                localxx[k][j][i].fx[2]=f20[10]*1.0e-6;            //He+
                localxx[k][j][i].fx[3]=f20[11]*0.1e-6;            //O2+
                localxx[k][j][i].fx[4]=0.05e-6*(f20[11]+f20[12]); //N2+
                localxx[k][j][i].fx[5]=0.03e-6*f20[12];           //NO+
                localxx[k][j][i].fx[6]=f20[13]*1.0e-6;            //N+

                if (f20[9]>0.0) numHi++;
                if (f20[10]>0.0) numHei++;

                //initialy set (rhoi uir), (rhoi uitheta), (rhoi uiphi) to a small constant
                for (s = 7; s < 10; s++) xx[k][j][i].fx[s] = 1.0e-23;

                /* ion and electron temperatures (K) (to be converted to pressures) */
                xx[k][j][i].fx[10] = f20[6];
                xx[k][j][i].fx[11] = f20[7];

                /* determine indexes between them the IRI density values non-zero */
                for (s = 0; s < sl; s++) {
                    if (i > 0) {
                        if(localxx[k][j][i].fx[s]> 0.0 && localxx[k][j][i-1].fx[s]<=0.0) i0[s]=i;
                        if(localxx[k][j][i].fx[s]<=0.0 && localxx[k][j][i-1].fx[s]> 0.0) im[s]=i-1;
                    }
                }

                /* read neutral density and temperature */
                for (s = 0; s < 11; s++) msisfstr >> f20[s];
                msisfstr.ignore(1,'\n');
                for (s = 3; s < 9; s++) {
                    if (f20[s] <= 0.0) f20[s]=denmin;
                }

                /* neutral density from MSIS00 in cm^-3 */
                xx[k][j][i].fx[12]=f20[4];  // O
                xx[k][j][i].fx[13]=f20[7];  // H
                xx[k][j][i].fx[14]=f20[3];  // He
                xx[k][j][i].fx[15]=f20[6];  // O2
                xx[k][j][i].fx[16]=f20[5];  // N2

                /* NO density from empirical model of Mitra, A. P., A review of 
                 * D-region processes in non-polar latitudes, JATP, 30, 1065-1114, 1968 */
                xx[k][j][i].fx[17]=(0.4*exp(-3700.0/f20[10])*f20[6]+5.0e-7*f20[4]);   // NO
                xx[k][j][i].fx[18]=f20[8];  // N

                /* initialy neutral rhon u_{n,r}, rhon u_{n,theta}, rhon u_{n,phi} set to a small constant */
                for (s = 19; s < 22; s++) xx[k][j][i].fx[s]=1.0e-19;

                if (f20[10] <= 0.0) {
                    cout<<"Negative neutral temperature from MSIS file at ";
                    cout<<"(i, j, k) = ("<<i<<", "<<j<<", "<<k<<") Tn = "<<f20[10]<<endl;
                    exit(-1);
                }

                /* conservative variable (pn = Nn kb Tn) (in unit Newton / cm^{-3})*/
                Nn=0.0;
                for (s = 0; s <sm; s++) Nn += xx[k][j][i].fx[12+s];
                xx[k][j][i].fx[22]=Nn*kb*f20[10];

                /* delta_B */
                xx[k][j][i].fx[23]=0.0;
                xx[k][j][i].fx[24]=0.0;
                xx[k][j][i].fx[25]=0.0;
            }

            /* extrapolation of densities to region where the density from IRI is zero */
            if (numHi < 3) {
                for (i = 0; i < a1; i++) localxx[k][j][i].fx[1] = localxx[k][j-1][i].fx[1];
                i0[1]=0; im[1]=0;
            }
            if (numHei < 3) {
                for (i = 0; i < a1; i++) localxx[k][j][i].fx[2] = localxx[k][j-1][i].fx[2];
                i0[2]=0; im[2]=0;
            }

            for (s = 0; s < sl; s++) {
                if (i0[s] > 0) {
                    for (i = i0[s]-1; i >= 0; i--) {
                        y1 = log(localxx[k][j][i+1].fx[s]);
                        y2 = log(localxx[k][j][i+2].fx[s]);
                        localxx[k][j][i].fx[s] = exp(y1+(rC[i]-rC[i+1])/(rC[i+2]-rC[i+1])*(y2-y1));
                        if (localxx[k][j][i].fx[s] < denmin) localxx[k][j][i].fx[s]=denmin;
                    }
                }
                if (im[s] > 0) {
                    ims=im[s];
                    for (i = im[s]; i > 1; i--) {
                        if (localxx[k][j][i].fx[s] < localxx[k][j][i-1].fx[s]) {
                            ims = i; break;
                        }
                    }
                    for (i = ims+1; i < a1; i++) {
                        y1 = log(localxx[k][j][i-2].fx[s]);
                        y2 = log(localxx[k][j][i-1].fx[s]);
                        localxx[k][j][i].fx[s] = exp(y1+(rC[i]-rC[i-2])/(rC[i-1]-rC[i-2])*(y2-y1));
                        if (localxx[k][j][i].fx[s] < denmin) localxx[k][j][i].fx[s]=denmin;
                        //if (s == 3 && localxx[k][j][i].fx[s] > localxx[k][j][i-1].fx[s])
                        //    localxx[k][j][i].fx[s] = localxx[k][j][i-2].fx[s];
                    }
                }
            }

            //convert ion and electron temperatures to conservative variables (pi = ni kb Ti), (pe = ne kb Te)
            //must be done after zero ion density regions have been filled
            for (i = xs; i < xs+xm; i++) {
                ne=0.0;
                for (s=0; s<sl; s++) {
                    xx[k][j][i].fx[s] = localxx[k][j][i].fx[s];

                    ne += xx[k][j][i].fx[s];
                }

                xx[k][j][i].fx[10] = ne*kb*xx[k][j][i].fx[10];
                xx[k][j][i].fx[11] = ne*kb*xx[k][j][i].fx[11];

                for (s=0; s<nvar; s++) {
                    if (isnan(xx[k][j][i].fx[s]) || isinf(xx[k][j][i].fx[s])) {
                        cout<<"Variable is Nan or inf at ("<<i<<", "<<j<<", "<<k
                            <<", "<<s<<") in input_iri_msis"<<endl;
                        exit(-1);
                    }
                }
            }
        }

        for (j = Nth/3; j < ys+ym; j++) {
            for (i = xs; i < xs+xm; i++) {
                for (s = 3; s < 7; s++) xx[k][j][i].fx[s]=xx[k][j-1][i].fx[s];
            }
        }
    }

    int kc;
    if (ys+ym == a2) {
        for (k = zs; k < zs+zm; k++) {
            kc = (k+a3/2) % a3;

            for (i = xs; i< xs+xm; i++) {
                for (s = 0; s < nvar; s++)
                    if (s != 24) xx[k][Nth][i].fx[s] = xx[kc][0][i].fx[s];
            }
        }
    }

    //release memory allocated
    DMDAVecRestoreArray(da, localX, &localxx);
    DMRestoreLocalVector(da,&localX);

    return 0;
}
