/*******************************************************************************
 *  array_alloc
 *  allocate memory for multi-dimension arrays
 * 
 *  Jiannan Tu
 *  5/21/2022
 ******************************************************************************/
#include <cmath>
using namespace std;

#include "param.h"
#include "funcdef.h"

void array_allocate(PetscInt xm, PetscInt ym, PetscInt zm)
{
    int  i, j, k;

//mean spherical coordinates, every process has a copy
    rr = new double[a1];
    rh = new double[a1+1];
    rC = new double[a1+1];
    rfavg =  new double[a1];

    theta =new double[a2];
    thetah =new double[a2];
    thetaC =new double[a2];

    phi = new double[a3];
    phih =new double[a3+1];

    rh_d3 = new double[a1];
    rh_d2 = new double[a1];

    sinth_h=new double[a2];
    sinth = new double[a2];

    zh    =new double[a1];
    cotth = new double[a2];

    rh2 = new double[a1];
    rfavg_dth=new double[a1];

// gravity
    gr  = new double[xm];
    
// ---- solar zenith angles
    zenith = new double*[zm];

    vt = new double*[zm];
    vp = new double*[zm];

// ---- Earth's rotation vector in geomagnetic coordinates
    rotat_r =new double*[zm];
    rotat_t=new double*[zm];
    rotat_p=new double*[zm];
    cenf_r =new double**[zm];
    cenf_t=new double**[zm];
    cenf_p=new double**[zm];

    for (k = 0; k < zm; k++) {
        zenith[k]=new double[ym];

        vt[k] = new double[ym];
        vp[k] = new double[ym];

        rotat_r[k] =new double[ym];
        rotat_t[k]=new double[ym];
        rotat_p[k]=new double[ym];

        cenf_r[k] =new double*[ym];
        cenf_t[k]=new double*[ym];
        cenf_p[k]=new double*[ym];

        for (j = 0; j < ym; j++) {
            cenf_r[k][j] =new double[xm];
            cenf_t[k][j]=new double[xm];
            cenf_p[k][j]=new double[xm];
        }
    }

    nust=new double***[zm];
    Omegae=new double**[zm];

    for (k = 0; k< zm; k++) {
        nust[k]=new double**[ym];
        Omegae[k]=new double*[ym];

        for (j = 0; j < ym; j++) {
            nust[k][j]=new double*[xm];
            Omegae[k][j]=new double[xm];

            for (i = 0 ; i < xm; i++) {
                nust[k][j][i]=new double[112];
            }
        }
    }

    nuin_omegae = new double**[a3];
    for (k = 0; k < a3; k++) {
        nuin_omegae[k] = new double*[a2];
        for (j = 0; j < a2; j++) nuin_omegae[k][j]=new double[a1];
    }
 
    for (k=0; k<zm; k++) {
        for (j=0; j<ym; j++) {
            vt[k][j]=-1.0e6;
            vp[k][j]=-1.0e6;
        }
    }

    Ftheta_Rface = new double*[xm];
    for (i = 0; i < xm; i++) Ftheta_Rface[i]=new double[nvar-3];

    Fphi_Rface = new double**[ym];
    for (j = 0; j < ym; j++) {
        Fphi_Rface[j]=new double*[xm];

        for (i = 0; i < xm; i++) Fphi_Rface[j][i]=new double[nvar-3];
    }

    rCsinC=new double*[a2];
    rfavg_costh = new double*[a2];
    rfavg_costh_dth_dph=new double*[a2];
    rfavg_sinth_dph=new double*[a2];
    rh_costh = new double*[a2];
    rh_costh_dth_dph = new double*[a2];
    dAtheta_dV = new double*[a2];

    for (j = 0; j < a2; j++) {
        rCsinC[j] = new double[a1];
        rfavg_costh[j] =new double[a1];
        rfavg_costh_dth_dph[j]=new double[a1];
        rfavg_sinth_dph[j]=new double[a1];
        rh_costh[j] = new double[a1];
        rh_costh_dth_dph[j] = new double[a1];
        dAtheta_dV[j] = new double[a1];
    }

    fluxn=new double**[xm];
    for (i =0 ; i < xm; i++) {
        fluxn[i]=new double*[92];
        for (j=0; j<92; j++) fluxn[i][j]=new double[4];
    }

    return;
}

void array_deallocate(PetscInt xm, PetscInt ym, PetscInt zm)
{
    int  i, j, k;

    delete[] rr;
    delete[] rh;
    delete[] rC;
    delete[] rfavg;

    delete[] theta;
    delete[] thetah;
    delete[] thetaC;

    delete[] phi;
    delete[] phih;

    delete[] rh_d3;
    delete[] rh_d2;

    delete[] sinth_h;
    delete[] sinth;

    delete[] zh;
    delete[] cotth;

    delete[] rh2;
    delete[] rfavg_dth;

    delete[] gr;

    for (k = 0; k < zm; k++) {
        delete[] zenith[k];

        delete[] vt[k];
        delete[] vp[k];

        for (j = 0; j < ym; j++) {
            delete[] cenf_r[k][j];
            delete[] cenf_t[k][j];
            delete[] cenf_p[k][j];
        }

        delete[] rotat_r[k];
        delete[] rotat_t[k];
        delete[] rotat_p[k];

        delete[] cenf_r[k];
        delete[] cenf_t[k];
        delete[] cenf_p[k];
    }
    delete[] zenith;
    delete[] vt;
    delete[] vp;

    delete[] rotat_r;
    delete[] rotat_t;
    delete[] rotat_p;
    delete[] cenf_r;
    delete[] cenf_t;
    delete[] cenf_p;

    for (k = 0; k < zm; k++) {
        for (j = 0; j < ym; j++) {
            for (i = 0; i < xm; i++) {
                delete[] nust[k][j][i];
            }
            delete[] nust[k][j];
            delete[] Omegae[k][j];
        }
        delete[] nust[k];
        delete[] Omegae[k];
    }
    delete[] nust;
    delete[] Omegae;

    for (k = 0; k < a3; k++) {
        for (j = 0; j < a2; j++) delete[] nuin_omegae[k][j];
        delete[] nuin_omegae[k];
    }
    delete[] nuin_omegae;

    for (j = 0; j < a2; j++) {
        delete[] rCsinC[j];
        delete[] rfavg_costh[j];
        delete[] rfavg_costh_dth_dph[j];
        delete[] rfavg_sinth_dph[j];
        delete[] rh_costh[j];
        delete[] rh_costh_dth_dph[j];
        delete[] dAtheta_dV[j];
    }
    delete[] rCsinC;
    delete[] rfavg_costh;
    delete[] rfavg_costh_dth_dph;
    delete[] rfavg_sinth_dph;
    delete[] rh_costh;
    delete[] rh_costh_dth_dph;
    delete[] dAtheta_dV;

    for (i = 0; i < xm; i++) delete[] Ftheta_Rface[i];
    delete[] Ftheta_Rface;

    for (j = 0; j < ym; j++) {
        for (i = 0; i < xm; i++) delete[] Fphi_Rface[j][i];

        delete[] Fphi_Rface[j];
    }
    delete[] Fphi_Rface;

    for (i=0; i< xm; i++) {
        for (j=0; j<92; j++) delete[] fluxn[i][j];
        delete[] fluxn[i];
    }
    delete[] fluxn;
}
