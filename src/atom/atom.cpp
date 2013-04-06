#include "atom.h"

Atom::Atom(Config*):
    nDimension(3),
    frozen(0),
    aPosition(zeros<rowvec>(nDimension)),
    aVelocity(zeros<rowvec>(nDimension)),
    aAcceleration(zeros<rowvec>(nDimension)),
    aDisplacement(zeros<rowvec>(nDimension))
{
}
