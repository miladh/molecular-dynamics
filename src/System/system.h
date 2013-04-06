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
#include <src/force/force.h>
#include <src/force/onebodyforce.h>
#include <src/force/twobodyforce.h>
#include <src/modifier/modifier.h>

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
    void addModifiers(Modifier* modifier);
    void addForces(Force* force);
    void evaluateSystemProperties();
    void setTopology();
    void atomCopy();
    void atomMove();
    void singleStep();
    void halfKick();
    void computeAccel() ;
    void restForce();
    void applyModifier();
    void applyForces(int atomI, int atomJ, int atomIsResident, int pairIsNotEvaluated);
    double getTemperature();
    int atomDidMove(rowvec r, int neighborID);
    int atomIsBoundary(rowvec r, int neighborID);

    Config* cfg;
    Atom** atoms;
    vector <Force*> forces;
    vector <Modifier*> modifiers;

    int procID, nProc, Nc, nX,nY,nZ;
    int nLocalResAtoms,nAtoms,nBounAtoms;
    double rCut;
    double cpu,comt;

    int cellList[NCLMAX], atomList[NEMAX];
    double dbuf[NDBUF],dbufr[NDBUF];

    ivec procVec, IDVec, neigID;
    ivec systemSize, nLocalCells, myParity;
    vec cellSize, subsysSize, origo;
    rowvec vdt;

    mat aPosition, aVelocity, aAcceleration;
    mat shiftVec;
    imat boundaryAtomList;

    MPI_Status status;
    string path;

    vec time,kinEnergy,potEnergy,totEnergy,temperature,pressure, displacement;

    double localKinEnergy, localPotEnergy, localPressure, localDisplacement;
    double density,latticeConstant, sigma, dt, T_0, volume;
    int stepLimit,stepAvg;
    int state;
    int loadState;

};

#endif // SYSTEM_H
