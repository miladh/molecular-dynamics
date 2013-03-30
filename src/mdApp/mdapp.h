#ifndef MDAPP_H
#define MDAPP_H

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/includes/defines.h>
#include <src/system/system.h>
#include <src/atom/atom.h>
#include <src/generator/generator.h>
#include <src/force/twobodyforce.h>
#include <src/modifier/modifier.h>


class MDApp
{
public:
    MDApp(const int& procID, const int& nProc);
    void runMDApp();
    void loadConfiguration(Config* cfg);

private:
    Config* cfg;
    string stateDir;
    stringstream outName;
    int loadState;
    int procID, nProc;
    int nLocalResAtoms, nAtoms;
    int modifierType, forceType;
    double tau,targetTemperature, T_0;

    void setModifierType(System *system);
    TwoBodyForce* setForceType();
};

#endif // MDAPP_H
