/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <cmath>
using namespace std;

#include "funcdef.h"

inline void lower_boundary_bc(Field ***xx, Field***gg, PetscInt j, PetscInt k)
{
    int     s;
    double  y1, y2;

    double rrb=(rr[0]-rr[1])/(rr[2]-rr[1]);

/*---------------------------------------------------------------*/
/*------- lower boundary conditions -----------------------------*/
/*---------------------------------------------------------------*/
    //density of ion and neutral species
    for (s=0; s< sl; s++) {
        y1=xx[k][j][1].fx[s];
        y2=xx[k][j][2].fx[s];
        gg[k][j][0].fx[s] = (y1 + rrb*(y2-y1));

        y1=xx[k][j][1].fx[20+s];
        y2=xx[k][j][2].fx[20+s];
        gg[k][j][0].fx[20+s] = (y1 + rrb*(y2-y1));
    }

    for (s=7; s <= 15; s++) {
        //us_r, us_theta, us_phi
        y1= xx[k][j][1].fx[s];
        y2= xx[k][j][2].fx[s];
        gg[k][j][0].fx[s] = (y1 + rrb*(y2-y1));
    }

    for (s=27; s <= 29; s++) {
        //un_r, un_theta, un_phi
        y1= xx[k][j][1].fx[s];
        y2= xx[k][j][2].fx[s];
        gg[k][j][0].fx[s] = (y1 + rrb*(y2-y1));
    }

    // TO+, TH+, THe+, Te
    for (s=16; s <= 19; s++) {
        y1= xx[k][j][1].fx[s];
        y2= xx[k][j][2].fx[s];
        gg[k][j][0].fx[s] = (y1 + rrb*(y2-y1));
    }

    // Tn
    y1= xx[k][j][1].fx[30];
    y2= xx[k][j][2].fx[30];
    gg[k][j][0].fx[30] = (y1 + rrb*(y2-y1));

    //Br,  Btheta and Bphi
    for (s = 31; s < 33; s++) {
        y1=xx[k][j][1].fx[s];
        y2=xx[k][j][2].fx[s];
        gg[k][j][0].fx[s] = (y1+rrb*(y2-y1));
    }
}

inline void upper_boundary_bc(Field ***xx, Field***gg, PetscInt j, PetscInt k, PetscInt yj, PetscInt zk)
{
    int     s, Nrm2=Nr-2;
    double  y1, y2, sgn, rho, rrt;

    rrt=(rr[Nr]-rr[Nrm])/(rr[Nrm]-rr[Nrm2]);

/*---------------------------------------------------------------*/
/*------- upper boundary conditions -----------------------------*/
/*---------------------------------------------------------------*/
    //density of ion and neutral species
    for (s=0; s<sl; s++) {
        y1= (xx[k][j][Nrm].fx[s]);
        y2= (xx[k][j][Nrm2].fx[s]);
        gg[k][j][Nr].fx[s] = (y1+rrt*(y1-y2));

        y1=(xx[k][j][Nrm].fx[20+s]);
        y2=(xx[k][j][Nrm2].fx[20+s]);
        gg[k][j][Nr].fx[20+s] = (y1+rrt*(y1-y2));
    }

    y1=xx[k][j][Nrm].fx[7];
    y2=xx[k][j][Nrm2].fx[7];
    gg[k][j][Nr].fx[7] = (y1+rrt*(y1-y2));

    y1=xx[k][j][Nrm].fx[10];
    y2=xx[k][j][Nrm2].fx[10];
    gg[k][j][Nr].fx[10] = (y1+rrt*(y1-y2));
    
    y1=xx[k][j][Nrm].fx[13];
    y2=xx[k][j][Nrm2].fx[13];
    gg[k][j][Nr].fx[13] = (y1+rrt*(y1-y2));

    //Br
    y1=xx[k][j][Nrm].fx[31];
    y2=xx[k][j][Nrm2].fx[31];
    gg[k][j][Nr].fx[31] = (y1+rrt*(y1-y2));

    if (vt[zk][yj] > -1.e5) {
        rho=0.0;
        for (s=0; s<sl; s++) rho += xx[k][j][Nr].fx[s]*ms[s];

        // us_theta                      
        gg[k][j][Nr].fx[8]  = rho*vt[zk][yj];
        gg[k][j][Nr].fx[11] = rho*vt[zk][yj];
        gg[k][j][Nr].fx[14] = rho*vt[zk][yj];

        // Btheta
        if (j<Nth/2) sgn=-1.0; else sgn=1.0;
        gg[k][j][Nr].fx[32] = sgn*vt[zk][yj]*sqrt(rho*mu0);
    }
    else {
        y1=xx[k][j][Nrm].fx[8];
        y2=xx[k][j][Nrm2].fx[8];
        gg[k][j][Nr].fx[8] = (y1+rrt*(y1-y2));

        y1=xx[k][j][Nrm].fx[11];
        y2=xx[k][j][Nrm2].fx[11];
        gg[k][j][Nr].fx[11] = (y1+rrt*(y1-y2));

        y1=xx[k][j][Nrm].fx[14];
        y2=xx[k][j][Nrm2].fx[14];
        gg[k][j][Nr].fx[14] = (y1+rrt*(y1-y2));

        y1=xx[k][j][Nrm].fx[32];
        y2=xx[k][j][Nrm2].fx[32];
        gg[k][j][Nr].fx[32] = (y1+rrt*(y1-y2));
    }

    if (vp[zk][yj] > -1.e5) {
        // us_phi
        rho=0.0;
        for (s=0; s<sl; s++) rho += xx[k][j][Nr].fx[s]*ms[s];

        gg[k][j][Nr].fx[9]  = rho*vp[zk][yj];
        gg[k][j][Nr].fx[12] = rho*vp[zk][yj];
        gg[k][j][Nr].fx[15] = rho*vp[zk][yj];

        // Bphi
        if (j<Nth/2) sgn=-1.0; else sgn=1.0;
        gg[k][j][Nr].fx[33] = sgn*vp[zk][yj]*sqrt(rho*mu0);
    }
    else {
        y1=xx[k][j][Nrm].fx[9];
        y2=xx[k][j][Nrm2].fx[9];
        gg[k][j][Nr].fx[9] = (y1+rrt*(y1-y2));

        y1=xx[k][j][Nrm].fx[12];
        y2=xx[k][j][Nrm2].fx[12];
        gg[k][j][Nr].fx[12] = (y1+rrt*(y1-y2));

        y1=xx[k][j][Nrm].fx[15];
        y2=xx[k][j][Nrm2].fx[15];
        gg[k][j][Nr].fx[15] = (y1+rrt*(y1-y2));

        y1=xx[k][j][Nrm].fx[33];
        y2=xx[k][j][Nrm2].fx[33];
        gg[k][j][Nr].fx[33] = (y1+rrt*(y1-y2));
    }

    // TO+, TH+, THe+ and Te
    for (s = 16; s <= 19; s++) {
        y1= log(xx[k][j][Nrm].fx[s]);
        y2= log(xx[k][j][Nrm2].fx[s]);
        gg[k][j][Nr].fx[s] = exp(y1 + rrt*(y1-y2));
    }

    //un_r, un_theta, un_phi
    for (s = 27; s <= 29; s++) {
        gg[k][j][Nr].fx[s] = xx[k][j][Nrm].fx[s];
        /*y1= xx[k][j][Nrm].fx[s];
        y2= xx[k][j][Nrm2].fx[s];
        gg[k][j][Nr].fx[s]= y1 + rrt*(y1-y2);*/
    }

    y1= log(xx[k][j][Nrm].fx[30]);
    y2= log(xx[k][j][Nrm2].fx[30]);
    gg[k][j][Nr].fx[30] = exp(y1 + rrt*(y1-y2));
}
