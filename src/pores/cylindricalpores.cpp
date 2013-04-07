#include <src/pores/cylindricalpores.h>

Cylindricalpores::Cylindricalpores(const int &procID, const int &nProc,const int &nLocalResAtoms):
    Pores(procID, nProc,nLocalResAtoms)
{
}

/***************************************************************
 * Name:
 * Description:
 */
void Cylindricalpores::makePores(Atom **atoms){

    vec2 distance;
    vec2 center;

    center(0) = Nc*latticeConstant/2;
    center(1) = Nc*latticeConstant/2;

    nLocalFrozenAtoms=0;
    for(int atom=0; atom < nLocalResAtoms; atom++ ){
        distance(0)=atoms[atom]->aPosition(0) - center(0);
        distance(1)=atoms[atom]->aPosition(1) - center(1);

        if (dot(distance, distance) > radius*radius){
            nLocalFrozenAtoms++;
            atoms[atom]->frozen = 1;
            atoms[atom]->aVelocity(0) = 0;
            atoms[atom]->aVelocity(1) = 0;
            atoms[atom]->aVelocity(2) = 0;
        }
    }

    MPI_Allreduce(&nLocalFrozenAtoms, &nFrozenAtoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(procID==0){
        cout << "------------------Making Pores--------------------" << endl;
        cout << "Number of frozen atoms: "<< nFrozenAtoms <<endl;
    }
}


/***************************************************************
 * Name:            loadConfiguration
 * Description:     Load system variables
 */
void Cylindricalpores::loadConfiguration(Config* cfg){
    radius = cfg->lookup("PoresSetting.cylindrical.radius");
    sigma = cfg->lookup("conversionFactors.sigma");
    latticeConstant = cfg->lookup("systemSettings.latticeConstant");
    Nc = cfg->lookup("systemSettings.Nc");
    radius/=sigma;
    latticeConstant/=sigma;
}
