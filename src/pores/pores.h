#ifndef PORES_H
#define PORES_H

#include <iostream>
#include <armadillo>
#include <libconfig.h++>
#include <src/atom/atom.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace arma;
using namespace std;
using namespace libconfig;

class Pores
{
public:
    Pores(const int &procID, const int &nProc, const int &nLocalResAtoms);
    virtual void makePores(Atom** atoms)=0;

protected:
    int procID, nProc;
    double sigma, latticeConstant;
    int nLocalResAtoms, Nc;
    int nLocalFrozenAtoms,nFrozenAtoms;
};

#endif // PORES_H
