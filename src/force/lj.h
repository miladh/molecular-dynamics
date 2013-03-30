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
    void calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated);
    void setParameters(Config* cfg);
    void restPotentialEnergy(){potEnergy = 0.0;}
    void restPressure(){pressure = 0.0; }
    double getPressure(){return pressure;}
    double getPotentialEnergy(){return potEnergy;}


private:
    double dr2,dri2,dri6,dr1,fcVal,f,vVal;
    double rCut2,rCuti2,rCuti6,r1;
    double Uc, dUc;
    double dRx,dRy,dRz;
    vec dr;
};

#endif // LJ_H
