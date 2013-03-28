#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/system/system.h>
#include <src/atom/atom.h>
#include <src/generator/generator.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"


int main(int argc, char* argv[])
{

    /********************************************************************
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

    int nLocalResAtoms;
    string velocityDist= cfg.lookup("initialVelocitySetting.velocityDist");
    /***********************************************************************/


    Atom **ArAtoms= new Atom*[NEMAX];
    for(int i=0; i<NEMAX; i++) {
        ArAtoms[i] = new Atom(&cfg);
    }

    Generator ArGenerator(procID, nProc);
    ArGenerator.loadConfiguration(&cfg);
    ArGenerator.fccLatticeGenerator(ArAtoms);
    ArGenerator.setVelocity(ArAtoms,velocityDist);
    nLocalResAtoms = ArGenerator.getNLocalResAtoms();

    System Ar(procID, nProc,  nLocalResAtoms ,ArAtoms);
    Ar.loadConfiguration(&cfg);
    Ar.MDrun();

    MPI_Finalize();
}
