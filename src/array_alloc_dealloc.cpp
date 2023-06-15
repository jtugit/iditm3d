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
    PetscInt  i, j, k;

//mean spherical coordinates, every process has a copy
    rr = new double[a1];
    rh = new double[a1+1];

    theta =new double[a2];
    thetah =new double[a2];

    phi = new double[a3];
    phih =new double[a3+1];

    zh    =new double[a1];

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

        rotat_r[k]=new double[ym];
        rotat_t[k]=new double[ym];
        rotat_p[k]=new double[ym];

        cenf_r[k]=new double*[ym];
        cenf_t[k]=new double*[ym];
        cenf_p[k]=new double*[ym];

        for (j = 0; j < ym; j++) {
            cenf_r[k][j]=new double[xm];
            cenf_t[k][j]=new double[xm];
            cenf_p[k][j]=new double[xm];
        }
    }

    r2sintheta = new double*[ym];
    cot_div_r = new double*[ym];
    rsin = new double*[ym];
    for (j = 0; j < ym; j++) {
        r2sintheta[j] =  new double[xm];
        cot_div_r[j] =  new double[xm];
        rsin[j] =  new double[xm];
    }

    grad_pe = new vector3D**[a3];
    for (k = 0; k < a3; k++) {
        grad_pe[k] = new vector3D*[a2];

        for (j = 0; j < a2; j++) grad_pe[k][j] = new vector3D[a2];
    }

    nust=new double***[zm];

    for (k = 0; k< zm; k++) {
        nust[k]=new double**[ym];

        for (j = 0; j < ym; j++) {
            nust[k][j]=new double*[xm];

            for (i = 0 ; i < xm; i++) {
                nust[k][j][i]=new double[38];
            }
        }
    }
 
    for (k=0; k<zm; k++) {
        for (j=0; j<ym; j++) {
            vt[k][j]=-1.0e6;
            vp[k][j]=-1.0e6;
        }
    }

    fluxn=new double**[xm];
    for (i =0 ; i < xm; i++) {
        fluxn[i]=new double*[92];
        for (j=0; j<92; j++) fluxn[i][j]=new double[4];
    }

    J11 = new double*[zm]; J12 = new double*[zm];

    Jiv11 = new double*[zm]; Jiv21 = new double*[zm];

    K11 = new double*[zm]; K12 = new double*[zm];
    K21 = new double*[zm]; K22 = new double*[zm];

    for (k = 0; k < zm; k++) {
        J11[k] = new double[ym]; J12[k] = new double[ym];

        Jiv11[k] = new double[ym]; Jiv21[k] = new double[ym];

        K11[k] = new double[ym]; K12[k] = new double[ym];
        K21[k] = new double[ym]; K22[k] = new double[ym];
    }

    J13 = new double[ym]; Jiv31 = new double[ym];
    K31 = new double[ym]; K32 = new double[ym];

    J21 = new double**[zm]; J22 = new double**[zm];
    J31 = new double**[zm]; J32 = new double**[zm];

    Jiv12 = new double**[zm]; Jiv13 = new double**[zm];
    Jiv22 = new double**[zm]; Jiv23 = new double**[zm];

    for (k = 0; k < zm; k++) {
        J21[k] = new double*[ym]; J22[k] = new double*[ym];
        J31[k] = new double*[ym]; J32[k] = new double*[ym];

        Jiv12[k] = new double*[ym]; Jiv13[k] = new double*[ym];
        Jiv22[k] = new double*[ym]; Jiv23[k] = new double*[ym];

        for (j = 0; j < ym; j++) {
            J21[k][j] = new double[xm]; J22[k][j] = new double[xm];
            J31[k][j] = new double[xm]; J32[k][j] = new double[xm];

            Jiv12[k][j] = new double[xm]; Jiv13[k][j] = new double[xm];
            Jiv22[k][j] = new double[xm]; Jiv23[k][j] = new double[xm];
        }
    }

    J23 = new double*[ym]; Jiv32 = new double*[ym];
    for (j = 0; j < ym; j++) {
        J23[j] = new double[xm]; Jiv32[j] = new double[xm];
    }

    K13 = new double[zm]; K23 = new double[zm];

    return;
}

void array_deallocate(PetscInt xm, PetscInt ym, PetscInt zm)
{
    PetscInt i, j, k;

    delete[] rr;
    delete[] rh;

    delete[] theta;
    delete[] thetah;

    delete[] phi;
    delete[] phih;

    delete[] zh;

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

    for (j = 0; j < ym; j++) {
        delete[] r2sintheta[j]; delete[] cot_div_r[j]; delete[] rsin[j];
    }
    delete[] r2sintheta; delete[] cot_div_r; delete[] rsin;

    for (k = 0; k < a3; k++) {
        for (j = 0; j < a2; j++) delete[] grad_pe[k][j];
        delete[] grad_pe[k];
    }
    delete[] grad_pe;

    for (k = 0; k < zm; k++) {
        for (j = 0; j < ym; j++) {
            for (i = 0; i < xm; i++) delete[] nust[k][j][i];
            delete[] nust[k][j];
        }
        delete[] nust[k];
    }
    delete[] nust;

    for (i=0; i< xm; i++) {
        for (j=0; j<92; j++) delete[] fluxn[i][j];
        delete[] fluxn[i];
    }
    delete[] fluxn;

    for (k = 0; k < zm; k++) {
        delete[] J11[k]; delete[] J12[k];
        delete[] Jiv11[k]; delete[] Jiv21[k];
        delete[] K11[k]; delete[] K12[k]; delete[] K21[k]; delete[] K22[k];
    }
    delete[] J11; delete[] J12;
    delete[] Jiv11; delete[] Jiv21;
    delete[] K11; delete[] K12; delete[] K21; delete[] K22;

    delete[] J13; delete[] Jiv31;
    delete[] K31; delete[] K32;

    for (k = 0; k < zm; k++) {
        for (j = 0; j < ym; j++) {
            delete[] J21[k][j]; delete[] J22[k][j];
            delete[] J31[k][j]; delete[] J32[k][j];

            delete[] Jiv12[k][j]; delete[] Jiv13[k][j];
            delete[] Jiv22[k][j]; delete[] Jiv23[k][j];
        }
        delete[] J21[k]; delete[] J22[k];
        delete[] J31[k]; delete[] J32[k];

        delete[] Jiv12[k]; delete[] Jiv13[k];
        delete[] Jiv22[k]; delete[] Jiv23[k];
    }
    delete[] J21; delete[] J22; delete[] J31; delete[] J32;

    delete[] Jiv12; delete[] Jiv13; delete[] Jiv22; delete[] Jiv23;

    for (j = 0; j < ym; j++) {
        delete[] J23[j]; delete[] Jiv32[j];
    }
    delete[] J23; delete[] Jiv32;

    delete[] K13; delete[] K23;

    return;
}
