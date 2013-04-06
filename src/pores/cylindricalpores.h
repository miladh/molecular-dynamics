#ifndef CYLINDRICALPORES_H
#define CYLINDRICALPORES_H

#include<src/pores/pores.h>

class Cylindricalpores : public Pores
{
public:
    Cylindricalpores(const int &procID, const int &nProc, const int &nLocalResAtoms);
    void makePores(Atom** atoms);
    void loadConfiguration(Config* cfg);

private:
    double radius;
};

#endif // CYLINDRICALPORES_H
