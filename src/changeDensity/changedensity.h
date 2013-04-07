#ifndef CHANGEDENSITY_H
#define CHANGEDENSITY_H

#include <iostream>
#include <armadillo>
#include <libconfig.h++>

#include <src/atom/atom.h>
#include <src/includes/defines.h>
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace arma;
using namespace std;
using namespace libconfig;

class ChangeDensity
{
public:
    ChangeDensity(const int &procID, const int &nProc, const int &nLocalResAtoms);

    void reduceDensity(Atom** atoms);
    void loadConfiguration(Config* cfg);
    int getnLocalResAtoms();

private:
    double reductionRatio;
    int procID, nProc;
    int nLocalResAtoms, nLocalRemovedAtoms,nRemovedAtoms;

    Config* cfg;
    Atom** newAtoms;
};

#endif // CHANGEDENSITY_H
