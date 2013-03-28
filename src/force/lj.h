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

private:
    int bintra;
    double dr2,dri2,dri6,dr1,fcVal,f,vVal;
    double rCut2,rCuti2,rCuti6,r1;
    double Uc, dUc;
    vec dr;


    double dRx,dRy,dRz;
};

#endif // LJ_H
