#ifndef NOFORCE_H
#define NOFORCE_H

#include <src/force/twobodyforce.h>

class NoForce : public TwoBodyForce
{
public:
    NoForce();
    void calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated);
    void setParameters(Config* cfg);
};

#endif // NOFORCE_H
