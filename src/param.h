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

#define kbmp 8.25481286619634103e3
#define data_dim 4
#define nvar     26
#define nion     7
#define nonu     7
#define denmin   10.0

/* variables of interest to be solved at each grid point */
typedef struct {
   double fx[nvar];
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
    double      RAgamma;

    double      rurb3;

    char dset_name[11]="dBPeNiViPi";
    char dset_diag[11]="diagnostic";
    char dset_pl[10]="prod_loss";

    Vec  U, V, W, Z;
} AppCtx;

const double pi=3.141592653589793238;
const double pi2=6.283185307179586476;
const double deg=57.29577951308232;
const double rad=0.017453292519943295;

//physical constants
const double e=1.6022e-19, ge=9.80665, Re=6371.2e3;
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

const double inv_gama=2.0/3.0; //adiabatic constant in neutral temperature eqs 

//extern double r0, n0, B0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0;

//number of grids along r, theta, phi
extern int  a1, a2, a3, a4;
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
extern double *rr, *rh, *rC, *rfavg, *theta, *thetah, *thetaC, *phi, *phih;
extern double dr, *rh_d3, *rh_d2, *sinth_h, *cotth, *sinth;
extern double **rCsinC, **rh_costh, **rh_costh_dth_dph;
extern double *zh, *rh2, **rfavg_costh, **rfavg_costh_dth_dph, **rfavg_sinth_dph, *rfavg_dth;
extern double **rsinth_dph;
extern double dth, dph, rCm1, rfavgm1;
extern double **dAtheta_dV;

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
extern double ****nust, ***Omegae;

extern double ***fluxn;

extern double **vt, **vp;

extern double **Ftheta_Rface, ***Fphi_Rface;

#endif  /* INC PARAM */
