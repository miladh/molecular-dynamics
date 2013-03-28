#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/atom/atom.h>
#include <src/includes/defines.h>
#include <src/fileManager/filemanager.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"


using namespace arma;
using namespace std;
using namespace libconfig;

class System
{
public:
    System(const int &procID, const int &nProc, const int &nLocalResAtoms, Atom **atoms);

    void MDrun();
    void initilizeParameters();
    void loadConfiguration(Config* cfg);
    void setTopology();
    void atomCopy();
    void atomMove();
    void singleStep();
    void halfKick();
    void computeAccel() ;
    int atomDidMove(rowvec r, int neighborID);
    int atomIsBoundary(rowvec r, int neighborID);

    int procID, nProc, Nc;
    int nAtoms,nLocalResAtoms;
    int nBounAtoms;
    double rCut;
    double Uc, dUc;    /* Potential cut-off parameters */
    double cpu,comt;
    double kinEnergy,potEnergy,totEnergy,temperature;

    int cellList[NCLMAX], atomList[NEMAX];
    double dbuf[NDBUF],dbufr[NDBUF];

    ivec procVec, IDVec, neigID;
    ivec nLocalCells, myParity;
    vec cellSize, subsysSize, origo;

    mat aPosition, aVelocity, aAcceleration;
    mat shiftVec;
    imat boundaryAtomList;

    Config* cfg;
    MPI_Status status;

    Atom** atoms;

    ivec systemSize;   /* Number of unit cells per processor */
    double density;     /* Number density of atoms (in reduced unit) */
    double dt;      /* Size of a time step (in reduced unit) */
    int stepLimit;      /* Number of time steps to be simulated */
    int stepAvg;        /* Reporting interval for statistical data */

    double T_0;
    string path;

    int nX,nY,nZ;

};

#endif // SYSTEM_H
