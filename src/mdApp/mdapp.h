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

#include <src/analysis/analysis.h>
#include <src/includes/defines.h>
#include <src/System/system.h>
#include <src/atom/atom.h>
#include <src/generator/generator.h>
#include <src/force/force.h>
#include <src/force/twobodyforce.h>
#include <src/force/onebodyforce.h>
#include <src/modifier/modifier.h>
#include <src/pores/pores.h>
#include <src/changeDensity/changedensity.h>


class MDApp
{
public:
    MDApp(const int& procID, const int& nProc);
    void runMDApp();
    void analyzeData();
    void loadConfiguration(Config* cfg);

private:
    Config* cfg;
    string stateDir;
    stringstream outName;
    int loadState, makePores,changeDensity;
    int procID, nProc;
    int nLocalResAtoms, nAtoms;
    int modifierType, forceType, poresShape;
    double tau,targetTemperature, T_0;

    void setModifierType(System *system);
    void setForceType(System *system);
    Pores* setPoresShape();
};

#endif // MDAPP_H
