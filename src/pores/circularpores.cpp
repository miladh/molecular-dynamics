#include <src/pores/circularpores.h>

Circularpores::Circularpores(const int &procID, const int &nProc,const int &nLocalResAtoms):
    Pores(procID, nProc, nLocalResAtoms)
{
}


/***************************************************************
 * Name:
 * Description:
 */
void Circularpores::makePores(Atom** atoms){
    nLocalFrozenAtoms = 0;
    radius = zeros<vec>(nPores);

    for (int i = 0; i < nPores; i++){
        radius(i) = rMin + randu()*(rMax - rMin);
        center = randu(3)*systemSize;
        for(int atom=0; atom < nLocalResAtoms; atom++ ){
            distance(0) = atoms[atom]->aPosition(0) - center(0);
            distance(1) = atoms[atom]->aPosition(1) - center(1);
            distance(2) = atoms[atom]->aPosition(2) - center(2);

            if (dot(distance, distance) < radius[i]*radius[i]){
                nLocalFrozenAtoms++;
                atoms[atom]->frozen = 1;
                atoms[atom]->aVelocity(0) = 0;
                atoms[atom]->aVelocity(1) = 0;
                atoms[atom]->aVelocity(2) = 0;
            }
        }
    }

    MPI_Allreduce(&nLocalFrozenAtoms, &nFrozenAtoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(procID==0){
        cout << "------------------Making Pores--------------------" << endl;
        cout << "Number of frozen atoms: "<< nFrozenAtoms <<endl;
        cout << "Radius of pores: "<< radius.t() << endl;
    }

}




/***************************************************************
 * Name:            loadConfiguration
 * Description:     Load system variables
 */
void Circularpores::loadConfiguration(Config* cfg){
    rMin = cfg->lookup("PoresSetting.circular.rMin");
    rMax = cfg->lookup("PoresSetting.circular.rMax");
    sigma = cfg->lookup("conversionFactors.sigma");
    latticeConstant = cfg->lookup("systemSettings.latticeConstant");
    Nc = cfg->lookup("systemSettings.Nc");
    nPores = cfg->lookup("PoresSetting.circular.nPores");

    rMin /=sigma;
    rMax /=sigma;
    latticeConstant/= sigma;
    systemSize << Nc << Nc << Nc;
    systemSize*=latticeConstant;
}
