#include "atom.h"

Atom::Atom(Config *cfg):
    nDimension(cfg->lookup("systemSettings.dim")),
    aPosition(zeros<rowvec>(nDimension)),
    aVelocity(zeros<rowvec>(nDimension)),
    aAcceleration(zeros<rowvec>(nDimension)),
    aDisplacement(zeros<rowvec>(nDimension)),
    localPressure(0)
{
}
