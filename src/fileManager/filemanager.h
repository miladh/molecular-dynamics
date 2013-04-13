#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <libconfig.h++>
#include <armadillo>
#include <iostream>

#include <src/atom/atom.h>
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace arma;
using namespace std;
using namespace libconfig;

class FileManager
{
public:
    FileManager(const int &procID, const int &nProc);
    void loadConfiguration(Config* cfg);
    void writeAtomProperties(const int &state, const int &nLocalResAtoms, const vec &origo, Atom **atoms);
    void writeSystemProperties(int numStates, const vec &t , const vec &Ek, const vec &Ep, vec const &Etot, const vec &T, const vec &P, const vec &D, const vec &meanV);
    void readDataFromFile(Atom** atoms);
    void writeAnalysisData(const vec &radius, const vec &numDensity, const vec &velo);

    int nLocalResAtoms;
private:
    int procID, nProc, stepLimit, nBins;
    double t_0,T_0;
    double sigma, epsilon, pressureFactor;

    string statisticsDir, statesDir,rawDataDir;
    stringstream outName;
    ofstream myfile;



};

#endif // FILEMANAGER_H
