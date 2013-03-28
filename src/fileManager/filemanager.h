#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <libconfig.h++>
#include <armadillo>
#include <iostream>
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
    FileManager(const int &procID, const int &nLocalResAtoms, const int &nProc);
    void loadConfiguration(Config* cfg);
    void writeAtomProperties(const int &state, const mat &aPosition, const mat &aVelocity);
    void readDataFromFile();
private:
    int procID, nLocalResAtoms, nProc, stepLimit;
    double t_0,T_0;
    double sigma, epsilon, pressureFactor;

    string statisticsDir, statesDir;
    stringstream outName;
    ofstream myfile;



};

#endif // FILEMANAGER_H