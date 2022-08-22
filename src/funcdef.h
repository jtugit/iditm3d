//Global scope function definition
//#ifndef INC_FUNCDEF_H
//#define INC_FUNCDEF_H

#include <string>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include "hdf5.h"
#include "param.h"

int input_param(AppCtx*);
void output_param(char **, int, int, int, AppCtx *);
int initialize(DM, Vec, AppCtx*);
void top_bc_vel(AppCtx *, PetscInt, PetscInt, PetscInt, PetscInt);
int imex_leap_frog(DM, Vec, Vec, Vec, AppCtx *);
int forward_scheme(DM da, Vec X, Vec Xn, AppCtx *params);

void array_allocate(PetscInt, PetscInt, PetscInt);
void array_deallocate(PetscInt, PetscInt, PetscInt);

int parameters(DM, Field***, AppCtx *);

int grids(DM, AppCtx*);

int euvflux_seg(DM, AppCtx*);

void dipole_magnetic(DM, AppCtx*);

void hyperslab_set(PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, 
        int, hsize_t[], hsize_t[], hsize_t[], hsize_t[]);
int hdf5parallelwrite(MPI_Comm, hsize_t*, hsize_t*, hsize_t*, hsize_t*, char*, char *, double*);
int hdf5parallelread(MPI_Comm, hsize_t*, hsize_t*, hsize_t*,  char*, char *, double*);

int input_psolutions(DM, Field ***, AppCtx*);
int output_solution(DM,Field***,AppCtx*);

double neu_cooling_rate(Field ***, Field***, int, int, int);

void prod_loss_rates(Field***, Field ***, int, int, int, int, int, int, double[], double[], double &, double &);

void solar_zenith(AppCtx*, PetscInt, PetscInt, PetscInt, PetscInt);

int SUN_08(int,int,int,int,int,double&,double&,double &,double&);
int dayno(int, int, int);
void update_timedate(AppCtx *);

double heat_fluxes_top(int, int, int, AppCtx *);

double cerfc(double);
double lgrg(double[], double[], PetscInt, double);
double expon_PetscIntegral(double, PetscInt);
void shapiro(DM da, Vec X, int cshap);
int data_attribute(char *, AppCtx *, int);
int print_auxiliary_hdf5(DM, Vec, Field ***, AppCtx *);

double expon_integral(double z, int n);

int rhsFunctions(TS ts, double ftime, Vec X, Vec F, void* ctx);

int boundary(DM da, Vec X, AppCtx *, int);
void boundary_bc(Field ***xx, int xs, int xm, int ys, int ym, int zs, int zm, AppCtx *);
void efd_boundary_bc(Field ***uu, int xs, int xm, int ys, int ym, int zs, int zm);
void boundary_V_T(Field ***xx, int xs, int xm, int ys, int ym, int zs, int zm);

//#endif  /* INC FUNCDEF */
