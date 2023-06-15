/************************************************************************
 header file: param_def.h
 Purpose: declare variables and arrays for global use (as externals)

 By Jiannan Tu 5/21/2022
************************************************************************/
#ifndef INC_PARAM_DEF_H
#define INC_PARAM_DEF_H

//double r0, n0, B0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0;

//number of grids along r, theta, phi
int a1, a2, a3;
int Nr, Nth, Np, Nrm, Nthm, Npm;

//number of ion and major neutral species
int sl, sm;

//time step
double dt, dt_half, dt2, alpha;

//change-mass ratios, normalized speed of light
double qms[7];

//normalized elementary charge
//normalized rotation frequency of Earth, normalized acceleration of gravity
double **rotat_r, **rotat_t, **rotat_p;
double ***cenf_r, ***cenf_t, ***cenf_p;

//r, theta, phi and geometrical parameters in geographic coordinates
//every process has a copy
double *rr, *rh, *theta, *thetah, *phi, *phih;
double dr, dth, dph, *zh;

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
double ****nust;

double ***fluxn;

double **vt, **vp;

//double ****Ps, ****Ls;

double **J11, **J12;     //J11[zk][yj], J12[zk][yj]
double *J13;             //J13[yj]
double ***J21, ***J22;   //J21[zk][yj][xi], J22[zk][yj][xi]
double **J23;            //J23[yj][xi]
double ***J31, ***J32;   //J31[zk][yj][xi], J32[zk][yj][xi]

double **Jiv11;                   //Jiv11[zk][yj]
double ***Jiv12, ***Jiv13;        //Jiv12[zk][yj][xi], Jiv13[zk][yj][xi]
double **Jiv21;                   //Jiv21[zk][yj]
double ***Jiv22, ***Jiv23;        //Jiv22[zk][yj][xi], Jiv23[zk][yj][xi]
double *Jiv31;                    //Jiv31[yj]
double **Jiv32;                   //J32[yj][xi]

double **K11, **K12;      //K11[zk][yj], K12[zk][yj]
double *K13;              //K13[zk]
double **K21, **K22;      //K21[zk][yj], K22[zk][yj]
double *K23;              //K23[zk]
double *K31, *K32;        //K31[yj], K32[yj]

double **r2sintheta, **cot_div_r, **rsin;

vector3D ***grad_pe;

/* evaluate normalization parameters from 4 basic parameters: r0, n0, B0, & mp */
double n0, B0, r0, v0, t0, E0, g0, p0, j00, e0, T0, q0, lamda0, beta0, gen, e; //, Omegae;
//normalized earth's rotational frequency
double w0n;

#endif  /* INC PARAM_DEF */
