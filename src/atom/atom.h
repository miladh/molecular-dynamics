#ifndef ATOM_H
#define ATOM_H

#include <armadillo>
#include <iostream>
#include <libconfig.h++>

using namespace arma;
using namespace std;
using namespace libconfig;

class Atom
{
public:
    Atom(Config* cfg);
    int nDimension;
    string aType;
    rowvec aPosition,aVelocity, aAcceleration,aDisplacement;
};

#endif // ATOM_H
