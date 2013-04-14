#include <iostream>
#include <libconfig.h++>
#include <src/mdApp/mdapp.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace std;
using namespace libconfig;

int main(int argc, char* argv[])
{
    Config cfg;
    if(argc <= 1){
        cfg.readFile("../md/src/config.cfg");
    }else{
        cfg.readFile(argv[1]);
    }

    int procID, nProc;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    if(procID==0){
        cout << "Starting molecular-dynamics" << endl;
        cout << "Number of processors: " << nProc <<endl;
    }


    MDApp mdApp(procID, nProc);
    mdApp.loadConfiguration(&cfg);
    mdApp.runMDApp();


    if(procID==0){
        cout << "Molecular-dynamics done!" << endl;
    }
    MPI_Finalize();

    return 0;
}

