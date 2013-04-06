#include "generator.h"

Generator::Generator(const int &procID, const int &nProc):
    procID(procID),
    nProc(nProc)
{
}

/************************************************************
Name:           fcc_lattice_creator
Description:    Creates a FCC lattice
*/
void Generator::fccLatticeGenerator(Atom **atoms)
{
    vec dR = zeros<vec>(3,1);
    vec dr = zeros<vec>(3,1);
    mat origAtom= zeros<mat>(4,3);

    ivec systemSize  = zeros<ivec>(3,1);
    vec subsysSize  = zeros<vec>(3,1);


    systemSize << Nc/nX << Nc/nY << Nc/nZ ;
    density = 4/pow(latticeConstant/sigma,3);

    for (int i=0; i < 3; i++){
        subsysSize(i) = (double)systemSize(i) / pow(density/4.0, 1.0/3.0);
    }

    // FCC atoms in the original unit cell
    origAtom << 0.0 << 0.0 << 0.0 << endr
             << 0.0 << 0.5 << 0.5 << endr
             << 0.5 << 0.0 << 0.5 << endr
             << 0.5 << 0.5 << 0.0 << endr;

    // Set up a face-centered cubic (fcc) lattice
    for (int i=0; i<3; i++){
        dr[i] = (double)subsysSize(i)/systemSize(i);
    }
    int aNum = 0;
    for (int nZ=0; nZ < systemSize(2); nZ++) {
        dR[2] = nZ*dr[2];
        for (int nY=0; nY<systemSize(1); nY++) {
            dR[1] = nY*dr[1];
            for (int nX=0; nX<systemSize(0); nX++) {
                dR[0] = nX*dr[0];
                for (int j=0; j < 4; j++) {
                    for (int i=0; i < 3; i++){
                        atoms[aNum]->aPosition(i) = dR[i] + dr[i]*origAtom(j,i);
                    }
                    aNum++;
                }
            }
        }
    }

    nLocalResAtoms = aNum;
    // Total # of atoms summed over processors
    MPI_Allreduce(&nLocalResAtoms,&nAtoms,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);


    //Write
    if(procID==0){
        cout << "------Initilize configuration-------" << endl;
        cout << "Number of Local atoms: "<< nLocalResAtoms << endl;
        cout << "Number of atoms:       "<< nAtoms << endl;
    }
}

/************************************************************
Name:           setVelocity
Description:    Sets initial velocity
*/
void Generator::setVelocity(Atom **atoms)
{
    setInitVelocityDistribution();
    vec localSumVelocities = zeros<vec>(3,1);
    vec sumVelocities = zeros<vec>(3,1);
    idum = idum - procID- time(NULL);
    srand(-idum);

    if(velocityDist =="uniform"){
        double v=2.0;
        for(int i = 0; i < nLocalResAtoms; i++){
            for(int j=0; j < 3 ; j++){
                atoms[i]->aVelocity(j) = -v+2*v*randu();
                localSumVelocities(j) +=atoms[i]->aVelocity(j);
            }
        }
    }else if(velocityDist=="normal"){
        double std =sqrt(Temperator/T_0);
        for(int i = 0; i < nLocalResAtoms; i++){
            for(int j=0; j < 3 ; j++){
                atoms[i]->aVelocity(j) = randn()*std;
                localSumVelocities(j) +=atoms[i]->aVelocity(j);
            }
        }
    }

    MPI_Allreduce(&localSumVelocities[0],&sumVelocities[0],3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    //Removing initial linear momentum
    sumVelocities/=nAtoms;
    for(int i=0; i<nLocalResAtoms; i++){
        for(int j=0; j<atoms[i]->nDimension; j++){
            atoms[i]->aVelocity(j) -= sumVelocities(j);
        }
    }
}

/************************************************************
Name:
Description:
*/
int Generator::getNLocalResAtoms()
{
    return nLocalResAtoms;
}

/************************************************************
Name:
Description:
*/
int Generator::getNAtoms()
{
    return nAtoms;
}



/************************************************************
Name:
Description:
*/
void Generator::setInitVelocityDistribution()
{
    switch (initVelocityDist) {
    case uniform:
        velocityDist = "uniform";
        break;
    case normal:
        velocityDist = "normal";
        break;
    }
}

/************************************************************
Name:          loadConfiguration
Description:    loadConfiguration
*/
void Generator::loadConfiguration(Config* cfg){
    Nc = cfg->lookup("systemSettings.Nc");
    nX = cfg->lookup("systemSettings.nX");
    nY = cfg->lookup("systemSettings.nY");
    nZ = cfg->lookup("systemSettings.nZ");
    latticeConstant = cfg->lookup("systemSettings.latticeConstant");
    sigma=cfg->lookup("conversionFactors.sigma");
    Temperator=cfg->lookup("initialVelocitySetting.initTemp");
    T_0=cfg->lookup("conversionFactors.T_0");
    idum = cfg->lookup("initialVelocitySetting.idum");
    initVelocityDist= cfg->lookup("initialVelocitySetting.initVelocityDist");
}
