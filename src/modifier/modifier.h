#ifndef MODIFIER_H
#define MODIFIER_H


class System;
#include <src/system/system.h>

class Modifier
{
public:
    Modifier(System* sys);
    virtual void apply() = 0;

protected:
    System* moleculeSystem;
};

#endif // MODIFIER_H
