#ifndef CIRCULARPORES_H
#define CIRCULARPORES_H

#include<src/pores/pores.h>

class Circularpores : public Pores
{
public:
    Circularpores(const int &procID, const int &nProc, const int &nLocalResAtoms);
    void makePores(Atom**atoms);
    void loadConfiguration(Config*cfg);

private:
    int nPores;
    double rMin, rMax;
    vec radius;
    vec3 distance, center, systemSize;



};

#endif // CIRCULARPORES_H
