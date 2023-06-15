/******************************************************************************
 *  input_param
 *   Input various parameters
 *
 * Jiannan Tu 5/21/2022
 ******************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

#include "param.h"

#define tag1 1

PetscErrorCode input_param(AppCtx *params)
{
    string fname;
    PetscMPIInt  rank, nprocs, file_free=0;
    MPI_Status status;

    fname = params->workdir + "/iditm3d.dat";

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

    if (rank==0) file_free=1;
    else {
        MPI_Recv(&file_free,1, MPI_INT, rank-1, tag1, MPI_COMM_WORLD, &status);
    }

    if (file_free==1) {
        fstream infstr(fname, fstream::in);
        if(!infstr) {
            cout << "Can't open file " << fname << endl;
            return -1;
        }
        infstr >> params->smod;
        infstr.ignore(200, '\n');
        infstr >> params->lognum;
        infstr.ignore(200, '\n');
        infstr >> params->iyr;
        infstr.ignore(200, '\n');
        infstr >> params->mon;
        infstr.ignore(200, '\n');
        infstr >> params->idate;
        infstr.ignore(200, '\n');
        infstr >> params->sec;
        infstr.ignore(200, '\n');
        infstr >> params->f107;
        infstr.ignore(200, '\n');
        infstr >> params->f107a;
        infstr.ignore(200, '\n');
        infstr >> params->Ap;
        infstr.ignore(200, '\n');
        infstr >> params->Bximf;
        infstr.ignore(200, '\n');
        infstr >> params->Byimf;
        infstr.ignore(200, '\n');
        infstr >> params->Bzimf;
        infstr.ignore(200, '\n');
        infstr >> params->SWDen;
        infstr.ignore(200, '\n');
        infstr >> params->vgsex;
        infstr.ignore(200, '\n');
        infstr >> params->vgsey;
        infstr.ignore(200, '\n');
        infstr >> params->vgsez;
        infstr.ignore(200, '\n');
        infstr >> params->UseAL;
        infstr.ignore(200, '\n');
        infstr >> params->ALindex;
        infstr.ignore(200, '\n');
        infstr >> r0;
        infstr.ignore(200, '\n');
        infstr >> n0;
        infstr.ignore(200, '\n');
        infstr >> B0;
        infstr.ignore(200, '\n');
        infstr >> params->ntot;
        infstr.ignore(200, '\n');
        infstr >> params->nout;
        infstr.ignore(200, '\n');
        infstr >> params->npre;
        infstr.ignore(200, '\n');
        infstr >> params->diag_step;
        infstr.ignore(200, '\n');
        infstr >> dt;
        infstr.ignore(200, '\n');
        infstr >> params->rb;
        infstr.ignore(200, '\n');
        infstr >> params->ru;
        infstr.ignore(200, '\n');
        infstr >> a1;
        infstr.ignore(200, '\n');
        infstr >> a2;
        infstr.ignore(200, '\n');
        infstr >> a3;
        infstr.ignore(200, '\n');
        infstr >> sl;
        infstr.ignore(200, '\n');
        infstr >> sm;
        infstr.ignore(200, '\n');
        infstr >> params->alpha;
        infstr.ignore(200, '\n');
        infstr >> params->outpdir;
        infstr.ignore(200, '\n');
        infstr >> params->prefln;
        infstr.ignore(200, '\n');
        infstr >> params->irifln;
        infstr.ignore(200, '\n');
        infstr >> params->msisfln;
        infstr.ignore(200, '\n');

        infstr.close();

        if(rank == 0 && params->lognum > 1) cout<<"lognum > 1. Not output to std device every time step"<<endl;
    }

    if(rank != nprocs-1) MPI_Send(&file_free,1, MPI_INT, rank+1, tag1, MPI_COMM_WORLD);

    if (params->smod == 0) {
        if (nprocs > 1) {
            cout << "When reading data from IRI and MSIS, must be run with one process"<<endl;
            return -1;
        }
    }
    //if (a3 % 2 !=0) {
    //    cout << "Number of grids, a3, along longitude must be even number"<< endl;
    //    return -1;
    //}

    if (!rank) {
        if(a1 % 5 != 0 || a2 % 5 !=0 || a3 % 5 !=0) {
            cout<<"Number of grids along r, theta and phi must be multiples of 5"<<endl;
            return -1;
        }
    }

    Nr =a1-1; Nrm=Nr-1;
    Nth=a2-1; Nthm=Nth-1;
    Np =a3-1; Npm=Np-1;

    return 0;
}

/* output_param
 *   Output various system constants & parameters, including background magnetic
 *   field. All vectors are in curvilinear coordinates. 
 *
 *   Jiannan Tu
 *   12/9/2013, 1/29/2014, 3/26/2020
 */
void output_param(char **argv, int nprocx, int nprocy, int nprocz, AppCtx *params)
{
    ofstream logfstr;
    int nprocs, ierr;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    time_t start_t = time(NULL);
    struct tm *now = localtime(&start_t );

    //open the log file to write execution records
    string fname = params->outpdir + "/iditm3d.log";
    logfstr.open(fname,fstream::out);
    if(!logfstr) {
        cout << "Can't open file " <<fname << endl;
        ierr=-1;
        MPI_Abort(MPI_COMM_WORLD,ierr);
    }

    cout <<"Program started: "
         <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday
         <<" " <<now->tm_hour<<":"<<now->tm_min<< ":"<<now->tm_sec<<endl;

    logfstr <<"Program started: "
            <<now->tm_year+1900<<"-"<<now->tm_mon+1<<"-"<<now->tm_mday
            <<" " <<now->tm_hour<<":"<<now->tm_min<< ":"<<now->tm_sec<<endl;

    logfstr << "Run command: " <<argv[0]<<endl << endl;
    logfstr<<"Number of processes: "<<nprocs<<endl;
    logfstr<<"Number of grids on r, theta and phi coordinates: "
           <<a1<<" "<<a2<<" "<<a3<<endl;
    logfstr << "Number of processes on r, theta, and phi coordinate: " 
            <<nprocx<<" "<<nprocy<<" "<<nprocz<<endl<<endl;

    cout<<"Number of processes: "<<nprocs<<endl;
    cout<<"Number of grids on r, theta and phi coordinates: "<<a1<<" "<<a2<<" "<<a3<<endl;
    cout << "Number of processes on r, theta and phi coordinate: " 
         <<nprocx<<" "<<nprocy<<" "<<nprocz<<endl;

/* basic parameters*/
    logfstr <<a1<<" "<<a2<<" "<<a3<<"   a1, a2 & a3: # of grids" <<endl;
    logfstr <<params->iyr <<" "<<params->mon<<" "<< params->idate<<" "
            <<params->sec<<"   iyr, mon, idate & UT sec"<<endl;
    logfstr <<params->Ap<<" "<<params->f107<<" "<< params->f107a
            <<"   Ap, f107 & f107a"<<endl;
    logfstr <<dt<<" "<<(params->rb-Re)*1.0e-3<<" "<<(params->ru-Re)*1.0e-3
            <<"   dt (s), zmin & zmax (km)"<<endl;
    logfstr <<params->ntot<<" "<< params->nout
            <<"   number of time step to run & number of time step between outputs"<< endl;
    logfstr <<sl<<" "<<sm<<"   sl & sm # of ion and neutral species" << endl;

    if (params->smod == 1) logfstr <<params->npre<<" time step at which continuous run starts"
                           <<endl;

    logfstr.close();

    return;
}
