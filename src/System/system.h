#ifndef SYSTEM_H
#define SYSTEM_H

class Modifier;

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/atom/atom.h>
#include <src/includes/defines.h>
#include <src/fileManager/filemanager.h>
#include <src/force/twobodyforce.h>
#include "src/modifier/modifier.h"

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
    System(const int &procID, const int &nProc, const int &nLocalResAtoms, const int &nAtoms, Atom **atoms);

    void simulateSystem();
    void initilizeParameters();
    void loadConfiguration(Config* cfg);
    void setTopology();
    void atomCopy();
    void atomMove();
    void singleStep();
    void halfKick();
    void computeAccel() ;
    void restForce();
    int atomDidMove(rowvec r, int neighborID);
    int atomIsBoundary(rowvec r, int neighborID);
    void evaluateSystemProperties();
    double getTemperature();
    void addModifiers(Modifier* modifier);
    void applyModifier();


    int procID, nProc, Nc;
    int nLocalResAtoms,nAtoms;
    int nBounAtoms;
    double rCut;
    double cpu,comt;
    vec time,kinEnergy,potEnergy,totEnergy,temperature,pressure, displacement;

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
    double density,latticeConstant,sigma;     /* Number density of atoms (in reduced unit) */
    double dt;      /* Size of a time step (in reduced unit) */
    int stepLimit;      /* Number of time steps to be simulated */
    int stepAvg;        /* Reporting interval for statistical data */

    double T_0;
    string path;
    double volume;
    int nX,nY,nZ;

    int state;
    int loadState;
    TwoBodyForce* force;
     vector <Modifier*> modifiers;
     rowvec vdt;

     double localKinEnergy, localPotEnergy, localPressure, localDisplacement;

};

#endif // SYSTEM_H
