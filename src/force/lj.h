#ifndef LJ_H
#define LJ_H


#include <libconfig.h++>
#include <armadillo>

#include <src/force/twobodyforce.h>

using namespace arma;
using namespace std;
using namespace libconfig;



class LJ : public TwoBodyForce
{
public:
    LJ();
};

#endif // LJ_H
