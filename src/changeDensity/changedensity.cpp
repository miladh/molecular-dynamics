#include <src/changeDensity/changedensity.h>

ChangeDensity::ChangeDensity(const int &procID, const int &nProc,const int &nLocalResAtoms):
    procID(procID),
    nProc(nProc),
    nLocalResAtoms(nLocalResAtoms)
{
}


/***************************************************************
Name:            loadConfiguration
Description:     Load system variables
*/
void ChangeDensity::reduceDensity(Atom** atoms)
{

    nLocalRemovedAtoms = 0;
    newAtoms = new Atom*[NEMAX];
    for(int i=0; i<NEMAX; i++) {
        newAtoms[i] = new Atom(cfg);
    }


    int i=0;
    for(int atom=0; atom < nLocalResAtoms; atom++ ){
        if(!atoms[atom]->frozen && randu() < reductionRatio){
            nLocalRemovedAtoms++;
        }
        else{
           newAtoms[i] = atoms[atom];
           i++;
        }
    }


    for(int i=0; i<NEMAX; i++) {
        atoms[i] = newAtoms[i];
    }

    nLocalResAtoms-=nLocalRemovedAtoms;
    MPI_Allreduce(&nLocalRemovedAtoms, &nRemovedAtoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(procID==0){
        cout << "---------------Reducing denisty-----------------" << endl;
        cout << "Number of removed atoms: "<< nRemovedAtoms <<endl;
    }


}


/***************************************************************
Name:            loadConfiguration
Description:     Load system variables
*/
int ChangeDensity::getnLocalResAtoms(){
    return nLocalResAtoms;
}

/***************************************************************
Name:            loadConfiguration
Description:     Load system variables
*/
void ChangeDensity::loadConfiguration(Config* cfg){
    this->cfg=cfg;
    reductionRatio = cfg->lookup("densitySettings.reductionRatio");
}


