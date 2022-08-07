/************************************************************************
 header file: param_def.h
 Purpose: declare variables and arrays for global use (as externals)

 By Jiannan Tu 5/21/2022
************************************************************************/
#ifndef INC_PARAM_DEF_H
#define INC_PARAM_DEF_H

//double r0, n0, B0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0;

//number of grids along r, theta, phi
int a1, a2, a3, a4;
int Nr, Nth, Np, Nrm, Nthm, Npm;

//number of ion and major neutral species
int sl, sm;

//time step
double dt, dt_half, dt2;

//change-mass ratios, normalized speed of light
double qms[7];

//normalized elementary charge
//normalized rotation frequency of Earth, normalized acceleration of gravity
double **rotat_r, **rotat_t, **rotat_p;
double ***cenf_r, ***cenf_t, ***cenf_p;

//r, theta, phi and geometrical parameters in geographic coordinates
//every process has a copy
double *rr, *rh, *rC, *rfavg, *theta, *thetah, *thetaC, *phi, *phih;
double dr, *rh_d3, *rh_d2;
double *sinth_h, **rCsinC, **rh_costh, **rh_costh_dth_dph;
double *zh, *rh2, **rfavg_costh, **rfavg_costh_dth_dph, **rfavg_sinth_dph, *rfavg_dth;
double *sinth;
double **rsinth_dph, *cotth;
double dth, dph, rCm1, rfavgm1;
double **dAtheta_dV;

//fast mode speed
extern double ****Vfs;

//gravitational acceleration 
double *gr;

//EUV emission flux, photo-absorption and photo-ionization cross sections & photon energies
//in given wavelength bins
double euvflux[37], segabs[37][5], segion[37][7], pene[37];

//for nighttime photo-absorption & ionization
double euvfluxn[4], segabsn[4][4], segionn[4][4];

//solar zenith
double **zenith;

//photo ionization rates
//double ***Cn;

//coefficients for collision frequencies, multiplied by n0*t0
double coe[15], coiO2[13], coiN2[13], coiO[13], coiH[13], coiHe[13], coiNO[13], coiN[13];
double con[10];

//ion and neutral collision frequencies */
double ****nust, ***nuin_omegae, ***Omegae;

double ***fluxn;

double **vt, **vp;

double **Ftheta_Rface, ***Fphi_Rface;

#endif  /* INC PARAM_DEF */
