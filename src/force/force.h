#ifndef FORCE_H
#define FORCE_H

#include <libconfig.h++>
#include <armadillo>

#include <src/atom/atom.h>


using namespace arma;
using namespace std;
using namespace libconfig;

class Force
{
public:
    Force();
    virtual void calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated) = 0;
    virtual void setParameters(Config* cfg)=0;
    virtual void restPotentialEnergy()=0;
    virtual double getPotentialEnergy()=0;
    virtual void restPressure()=0;
    virtual double getPressure()=0;

protected:
    double potEnergy, pressure;
};

#endif // FORCE_H
