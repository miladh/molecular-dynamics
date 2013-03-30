#ifndef NOFORCE_H
#define NOFORCE_H

#include <src/force/twobodyforce.h>

class NoForce : public TwoBodyForce
{
public:
    NoForce();
    void calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated){}
    void setParameters(Config* cfg){}
    void restPressure(){pressure = 0.0; }
    void restPotentialEnergy(){potEnergy = 0.0;}
    double getPotentialEnergy(){return potEnergy;}
    double getPressure(){return pressure;}
};

#endif // NOFORCE_H
