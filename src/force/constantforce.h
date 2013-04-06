#ifndef SIMPLEFORCE_H
#define SIMPLEFORCE_H

#include <src/force/onebodyforce.h>

class ConstantForce : public OneBodyForce
{
public:
    ConstantForce();

    void calculateAndApplyForce(Atom* atom, Atom*, int, int);
    void setParameters(Config* cfg);
    void restPotentialEnergy(){potEnergy = 0.0;}
    void restPressure(){pressure = 0.0; }
    double getPressure(){return pressure;}
    double getPotentialEnergy(){return potEnergy;}

private:
    double forceMagnitude;
    int dir1,dir2;
};

#endif // SIMPLEFORCE_H
