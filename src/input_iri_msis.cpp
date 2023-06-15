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
    double   f20[14];
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

                // O+, H+, He+, O2+, N2+, NO+, N+ density (m^{-3}) Density from IRI in m^-3 
                localxx[k][j][i].fx[0]=f20[8];             //O+
                localxx[k][j][i].fx[1]=f20[9];             //H+
                localxx[k][j][i].fx[2]=f20[10];            //He+
                localxx[k][j][i].fx[3]=f20[11];            //O2+
                localxx[k][j][i].fx[4]=0.03*(f20[11]+f20[12]); //N2+
                localxx[k][j][i].fx[5]=0.03*f20[12];       //NO+
                localxx[k][j][i].fx[6]=f20[13];            //N+

                if (f20[9]>0.0) numHi++;
                if (f20[10]>0.0) numHei++;

                //initialy set three components of O+, H+, and He+ velocity to a small constant (m/s)
                for (s = 7; s < 16; s++) xx[k][j][i].fx[s] = 0.5/v0;

                /* neutral, ion O+, H+, and He+, electron temperatures (K) */
                xx[k][j][i].fx[30] =  f20[5]/T0;
                for (s = 16; s < 19; s++) xx[k][j][i].fx[s] = f20[6]/T0;
                xx[k][j][i].fx[19] = f20[7]/T0;

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

                /* neutral density from MSIS00 in m^-3 in log scale */
                xx[k][j][i].fx[20]=log(f20[4]*1.0e6/n0);  // O
                xx[k][j][i].fx[21]=log(f20[7]*1.0e6/n0);  // H
                xx[k][j][i].fx[22]=log(f20[3]*1.0e6/n0/5.0);  // He
                xx[k][j][i].fx[23]=log(f20[6]*1.0e6/n0);  // O2
                xx[k][j][i].fx[24]=log(f20[5]*1.0e6/n0);  // N2

                /* NO density from empirical model of Mitra, A. P., A review of 
                 * D-region processes in non-polar latitudes, JATP, 30, 1065-1114, 1968 */
                xx[k][j][i].fx[25]=log((0.4*exp(-3700.0/f20[10])*f20[6]+5.0e-7*f20[4])*1.0e6/n0);   // NO
                xx[k][j][i].fx[26]=log(f20[8]*1.0e6/n0);  // N

                /* initialy neutral u_{n,r}, u_{n,theta}, u_{n,phi} set to a small constant (m/s) */
                for (s = 27; s < 30; s++) xx[k][j][i].fx[s]=0.1/v0;

                if (f20[10] <= 0.0) {
                    cout<<"Negative neutral temperature from MSIS file at ";
                    cout<<"(i, j, k) = ("<<i<<", "<<j<<", "<<k<<") Tn = "<<f20[10]<<endl;
                    exit(-1);
                }
                xx[k][j][i].fx[30]=f20[10]/T0;

                /* delta_B */
                xx[k][j][i].fx[31]=1.0e-9; xx[k][j][i].fx[32]=1.0e-9; xx[k][j][i].fx[33]=1.0e-9;

                xx[k][j][i].fx[34]=1.0e-4; xx[k][j][i].fx[35]=1.0e-4; xx[k][j][i].fx[36]=1.0e-4;
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
                        localxx[k][j][i].fx[s] = exp(y1+(rr[i]-rr[i+1])/(rr[i+2]-rr[i+1])*(y2-y1));
                        if (localxx[k][j][i].fx[s] < denmin/n0) localxx[k][j][i].fx[s]=denmin/n0;
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
                        localxx[k][j][i].fx[s] = exp(y1+(rr[i]-rr[i-2])/(rr[i-1]-rr[i-2])*(y2-y1));
                        if (localxx[k][j][i].fx[s] < denmin/n0) localxx[k][j][i].fx[s]=denmin/n0;
                    }
                }
            }

            for (i = xs; i < xs+xm; i++) {
                for (s=0; s<sl; s++) {
                    if (localxx[k][j][i].fx[s] <= 0.0) xx[k][j][i].fx[s] = denmin/n0;
                    else xx[k][j][i].fx[s] = localxx[k][j][i].fx[s]/n0;
                }

                for (s=0; s<a4; s++) {
                    if (isnan(xx[k][j][i].fx[s]) || isinf(xx[k][j][i].fx[s])) {
                        cout<<"Variable is Nan or inf at ("<<i<<", "<<j<<", "<<k<<", "<<s<<") in input_iri_msis"<<endl;
                        exit(-1);
                    }
                }
            }
        }
    }

    //release memory allocated
    DMDAVecRestoreArray(da, localX, &localxx);
    DMRestoreLocalVector(da,&localX);

    return 0;
}
