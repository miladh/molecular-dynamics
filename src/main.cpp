#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include<src/system/system.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

int main(int argc, char* argv[])
{
    Config cfg;
    if(argc <= 1){
          cfg.readFile("../md/src/config.cfg");
    }
    else{
        cfg.readFile(argv[1]);
    }

    System Ar;
    Ar.loadConfiguration(&cfg);
    Ar.MDrun();

}
