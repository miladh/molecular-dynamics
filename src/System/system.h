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
    void add1BForces(Force* force);
    void add2BForces(Force* force);
    void evaluateSystemProperties();
    void setTopology();
    void atomCopy();
    void atomMove();
    void singleStep();
    void halfKick();
    void computeAccel() ;
    void restForce();
    void applyModifier();
    void apply1BForces(int atomI);
    void apply2BForces(int atomI, int atomJ, int atomIsResident, int pairIsNotEvaluated);
    double getTemperature();
    int atomDidMove(rowvec r, int neighborID);
    int atomIsBoundary(rowvec r, int neighborID);

    Config* cfg;
    Atom** atoms;
    vector <Force*> forces1B;
    vector <Force*> forces2B;
    vector <Modifier*> modifiers;

    int procID, nProc, Nc, nX,nY,nZ;
    int nLocalResAtoms,nAtoms,nBounAtoms;
    int stepLimit,stepAvg,actingStep,state;
    int loadState;
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

    vec time,kinEnergy,potEnergy,totEnergy,temperature,pressure, displacement,meanVelocity;

    double localKinEnergy, localPotEnergy, localPressure, localDisplacement,localMeanVelocity;
    double density,latticeConstant, sigma, dt, T_0, volume;


};

#endif // SYSTEM_H
