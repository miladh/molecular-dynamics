#include "pores.h"

Pores::Pores(const int &procID, const int &nProc,const int &nLocalResAtoms):
    procID(procID),
    nProc(nProc),
    nLocalResAtoms(nLocalResAtoms)
{
}
