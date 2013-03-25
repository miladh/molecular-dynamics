#ifndef MD_H
#define MD_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>
#include <src/includes/defines.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"


using namespace arma;
using namespace std;
using namespace libconfig;

Config cfg;

int procID, nProc, N;
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
vec cellSize, subsysSize;

mat aPosition, aVelocity, aAcceleration;
mat shiftVec;
imat boundaryAtomList;


MPI_Status status;








/* Input data-----------------------------------------------------------
----------------------------------------------------------------------*/
ivec systemSize;   /* Number of unit cells per processor */
double density;     /* Number density of atoms (in reduced unit) */
double initTemp;    /* Starting temperature (in reduced unit) */
double dt;      /* Size of a time step (in reduced unit) */
int stepLimit;      /* Number of time steps to be simulated */
int stepAvg;        /* Reporting interval for statistical data */



/*****************************
 * Functions
 *****************************/
void initilizeParameters();
void loadConfiguration();
void setTopology();
void initializeConfiguration();
void atomCopy();
void atomMove();
void singleStep();
void halfKick();
void computeAccel() ;

int atomDidMove(rowvec r, int neighborID);
int atomIsBoundary(rowvec r, int neighborID);


#endif // MD_H
