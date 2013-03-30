#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/mdApp/mdapp.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

int main(int argc, char* argv[])
{
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

    if(procID==0){
        cout << "Starting molecular-dynamics" << endl;
    }

    MDApp mdApp(procID, nProc);
    mdApp.loadConfiguration(&cfg);
    mdApp.runMDApp();
    if(procID==0){
        cout << "Molecular-dynamics done!" << endl;
    }
    MPI_Finalize();
}
