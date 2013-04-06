#ifndef TWOBODYFORCE_H
#define TWOBODYFORCE_H

#include <src/force/force.h>

class TwoBodyForce: public Force
{
public:
    TwoBodyForce();

protected:
    double rCut;
};

#endif // TWOBODYFORCE_H
