#ifndef GENERATOR_H
#define GENERATOR_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>
#include "src/atom/atom.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace arma;
using namespace std;
using namespace libconfig;

class Generator
{
public:
    Generator(const int &procID, const int &nProc);

    void loadConfiguration(Config* cfg);
    void fccLatticeGenerator(Atom **atoms);
    void setVelocity(Atom **atoms, string distribution);
    int getNLocalResAtoms();

private:
    int procID,nProc;
    int Nc, nAtoms, nLocalResAtoms;
    double sigma,Temperator,T_0,latticeConstant, density;
};

#endif // GENERATOR_H
