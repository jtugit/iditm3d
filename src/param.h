/************************************************************************
 header file: param.h
 Purpose: declare variables and arrays for global use (as externals)

 By Jiannan Tu 5/21/2022
************************************************************************/
#ifndef INC_PARAM_H
#define INC_PARAM_H

#include "petscdm.h"
#include "petscdmda.h"
#include "petscts.h"
#include <string>
#include <vector>
using namespace std;

#define kbmp 8.25481286619634103e3
#define data_dim 4
#define a4      37
#define nion     7
#define nonu     7
#define denmin   10.0e4

/* variables of interest to be solved at each grid point */
typedef struct {
   double fx[a4];
} Field;

/* User defined variables. Used in solving for variables of interest */
typedef struct {
    int         npre, nout, ntot, smod;
    std::string workdir;
    std::string outpdir;
    std::string prefln;
    std::string irifln, msisfln;
    time_t      start_t;

    int         iyr, mon, idate;
    double      sec, f107, f107a, Ap;

    //lower and upper radial distance
    double      rb, ru;

    int         ndt;
    int         diag_step;
    int         lognum;     //time step interval to output to std

    double      Bximf;
    double      Byimf;
    double      Bzimf;
    double      SWDen;
    double      vgsex;
    double      vgsey;
    double      vgsez;
    int         UseAL;
    double      ALindex;

    double      rurb3;

    char dset_name[11]="dBPeNiViPi";
    char dset_diag[11]="diagnostic";
    char dset_pl[10]="prod_loss";

    Vec  U;

    double ftime;
} AppCtx;

const double pi=3.141592653589793238;
const double pi2=6.283185307179586476;
const double deg=57.29577951308232;
const double rad=0.017453292519943295;

//physical constants
const double q=1.6022e-19, ge=9.80665, Re=6371.2e3;
const double kb=1.3807e-23, eps0=8.85e-12, mu0=1.25663706e-6;
const double gm=3.984997527e14, cspd=3.0e8;
const double w0=7.27220517e-5; //Earth's rotation rate
const double cspeed2=9.0e16*1.0e-6;  //speed of light reduced by a factor of 0.002

const double mp = 1.6605390666e-27;

const double me=9.1093837015e-31;
//electron mass normalized by H mass, i.e., mass in AMU
const double ame=0.0005485799090624057;

const double ms[7]={2.656762874216004e-26, 1.673773562960802e-27, 6.64647366797316e-27, 
                    5.313525748432007e-26, 4.651734508829244e-26, 4.9826301286306257e-26,
                    2.325867254414622e-26};
/* mass of O+, H+, He+, O2+, N2+, NO+, N+ or O, H, He, O2, N2, NO, N normalized 
 * by proton mass, i.e. mass in atomic mass unit (AMU) */
const double ams[7]={15.9994, 1.00797, 4.0026, 31.9988, 28.0134, 30.0061, 14.0067};
//const double e_ms[3]={6.03064735904794e6, 9.57238204411478e7, 2.4106015912173025e7};
//const double ms_kb[3]={1.9242144377605592e-3, 1.21226447668632e-4, 4.813843462001695e-4};
//const double kb_ms[3]={519.6925978602576, 8249.025020789712, 2077.341651440076};

//average photoelectron energies (eps, in eV) in the parameterization of 
//photoelectron heating
const double epi[16]={169.4, 99.07, 69.43, 54.50, 48.37, 43.63, 45.93, 40.88,
                      40.81, 37.38, 33.68, 34.22, 29.62, 29.65, 25.38, 23.69};

//coefficients used in the parameterization of photoelectron heating
const double cc1[7]={1.468, 9.229e-1, 4.956e-2, -1.897e-2, -3.934e-3, -2.634e-4, 
                     -5.980e-6};
const double cc2[7]={1.020, 1.540e-2, -6.858e-3, -8.528e-3, -2.052e-3, -1.634e-4,
                     -4.314e-6};

//for use to interpolate night time euv flux
const double zaltnt[4][2]={{90.0, 220.0},{170.0, 320.0},{138.0, 280.0},{65.0, 220.0}};
const double thetant[4][5]={{0.0, 90.0, 120.0, 150.0, 180.0},
                            {0.0, 90.0, 105.0, 120.0, 135.0},
                            {0.0, 102.0, 120.0, 130.0, 130.0},
                            {0.0, 90.0, 120.0, 150.0, 180.0}};

//number of grids along r, theta, phi
extern int  a1, a2, a3;
extern int Nr, Nth, Np, Nrm, Nthm, Npm;

//number of ion and major neutral species
extern int sl, sm;

//time step
extern double dt, dt_half, dt2;

//change-mass ratios, speed of light
extern double qms[7];

//rotation acceleration of Earth, Centrifugal force
extern double **rotat_r, **rotat_t, **rotat_p;
extern double ***cenf_r, ***cenf_t, ***cenf_p;

//r, theta, phi and geometrical parameters in geographic coordinates
//every process has a copy
extern double *rr, *rh, *theta, *thetah, *phi, *phih;
extern double dr, dth, dph, *zh;

//gravitational acceleration 
extern double *gr;

//EUV emission flux, photo-absorption and photo-ionization cross sections & photon energies
//in given wavelength bins
extern double euvflux[37], segabs[37][5], segion[37][7], pene[37];

//for nighttime photo-absorption & ionization
extern double euvfluxn[4], segabsn[4][4], segionn[4][4];

//solar zenith
extern double **zenith;

//coefficients for collision frequencies, multiplied by n0*t0
extern double coe[15], coiO2[13], coiN2[13], coiO[13], coiH[13], coiHe[13], coiNO[13], coiN[13];
extern double con[10];

//ion and neutral collision frequencies */
extern double ****nust;

extern double ***fluxn;

extern double **vt, **vp;

extern double ****Ps, ****Ls, ***Qee, ***Qeuv;

typedef struct{
    vector<double> J11, J12;     //J11[zk][yj], J12[zk][yj] - use 1-D indexing kj=k*ym+j
    vector<double> J13;          //J13[yj]
    vector<double> J21, J22;     //J21[zk][yj][xi], J22[zk][yj][xi] - use 1-D indexing kji=k*ym*xm + j*xm+i
    vector<double> J23;          //J23[yj][xi] - use 1-D indexing ji = j*xm+i
    vector<double> J31, J32;     //J31[zk][yj][xi], J32[zk][yj][xi] - use 1-D indexing kji=k*ym*xm + j*xm+i
} Jmatrix;

typedef struct{
    vector<double> Jiv11;               //Jiv11[zk][yj] - 1-D indexing kj=k*ym+j
    vector<double> Jiv12, Jiv13;        //Jiv12[zk][yj][xi], Jiv13[zk][yj][xi] - 1-D indexing kji=k*ym*xm + j*xm+i
    vector<double> Jiv21;               //Jiv21[zk][yj] - 1-D indexing kj=k*ym+j
    vector<double> Jiv22, Jiv23;        //Jiv22[zk][yj][xi], Jiv23[zk][yj][xi] - 1-D indexing kji=k*ym*xm + j*xm+i
    vector<double> Jiv31;               //Jiv31[yj]
    vector<double> Jiv32;               //J32[yj][xi] - use 1-D indexing ji = j*xm+i
} Jinvmatrix;

typedef struct{
    vector<double> K11, K12; //K11[zk][yj], K12[zk][yj] - 1-D indexing kj=k*ym+j
    vector<double> K13;      //K13[zk]
    vector<double> K21, K22; //K21[zk][yj], K22[zk][yj] - use 1-D indexing ji = j*a1+i
    vector<double> K23;      //K23[zk]
    vector<double> K31, K32; //K31[yj], K32[yj]
} Kmatrix;

extern Jmatrix Jmat;
extern Jinvmatrix Jinv;
extern Kmatrix Kmat;
extern double **r2sintheta, **cot_div_r, **rsin;

typedef struct {
    double r, t, p;
} vector3D;

extern vector3D ***grad_pe;

/* evaluate normalization parameters from 4 basic parameters: r0, n0, B0, & mp */
extern double n0, B0, r0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0, gen, e, Omegae;
//normalized earth's rotational frequency
extern double w0n;

const PetscInt dfill[676] = {
    //  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 1
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 2
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 3
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 4
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 5
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 6
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  // 7
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  // 8
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  // 9
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,  // 10
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,  // 11
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 12
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 13
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 14
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 15
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  // 17
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  // 18
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  // 19
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  // 20
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  // 21
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,  // 22
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  // 23
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  // 24
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  // 25
    };
const PetscInt ofill[676] = {
    //  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 16 18 29 21 22 23 24 25
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 1
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 2
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 3
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 4
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 5
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 6
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 7
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 8
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 9
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 10
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 11
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 12
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 13
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 14
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 15
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  // 17
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  // 18
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  // 19
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  // 20
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  // 21
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  // 22
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  // 23
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  // 24
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  // 25
    };

#endif  /* INC PARAM */
