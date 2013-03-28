#ifndef TWOBODYFORCE_H
#define TWOBODYFORCE_H

#include <libconfig.h++>
#include <armadillo>

#include <src/atom/atom.h>


using namespace arma;
using namespace std;
using namespace libconfig;


class TwoBodyForce
{
public:
    TwoBodyForce();

    virtual void calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated) = 0;
    virtual void setParameters(Config* cfg)=0;

protected:
    double rCut;
};

#endif // TWOBODYFORCE_H
