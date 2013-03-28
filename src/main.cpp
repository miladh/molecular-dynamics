#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/system/system.h>
#include <src/atom/atom.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"


int main(int argc, char* argv[])
{

    /************************************************************
    Initializing
    */
    int procID, nProc;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    Config cfg;
    for(int i=0; i<nProc; i++){
        if(procID==i){
            if(argc <= 1){
                cfg.readFile("../md/src/config.cfg");
            }else{
                cfg.readFile(argv[1]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    /***************************************************************/



    Atom **ArAtoms= new Atom*[NEMAX];
    for(int i=0; i<NEMAX; i++) {
        ArAtoms[i] = new Atom(&cfg);
    }


    System Ar(procID, nProc, ArAtoms);
    Ar.loadConfiguration(&cfg);
    Ar.MDrun();


    MPI_Finalize();
}
