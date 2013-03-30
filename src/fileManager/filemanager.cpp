#include "filemanager.h"

FileManager::FileManager(const int &procID, const int &nLocalResAtoms, const int &nProc):
    procID(procID),
    nLocalResAtoms(nLocalResAtoms),
    nProc(nProc),
    stepLimit(stepLimit)
{
}

/************************************************************
Name:           writeToFile
Description:    Writes the states to files
*/
void FileManager::writeAtomProperties(const int &state, const vec &origo, Atom** atoms){

    outName << statesDir << "state" << state << ".xyz";
    myfile.open(outName.str().c_str());

    if(procID==0){
        myfile << nLocalResAtoms*nProc << endl;
        myfile << "Argon atoms" << endl;
    }


    for(int node=0; node<nProc ; node++){
        MPI_Barrier(MPI_COMM_WORLD);
        if(procID==node){
            for(int i=0;  i < nLocalResAtoms; i++){
                atoms[i]->aType="Ar";
                myfile << atoms[i]->aType <<"  "<<join_rows(atoms[i]->aPosition-origo.t(), atoms[i]->aVelocity);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}


/************************************************************
Name:
Description:
*/
void FileManager::writeSystemProperties(int numStates, const vec &t,
                                        const vec &Ek, const vec &Ep, vec const &Etot, const vec &T,const vec &P, const vec &D){

    outName << statisticsDir << "/statistics.txt";
    myfile.open (outName.str().c_str());
    myfile << "Time  "  <<"Kinetic " << "  Potential  "
           <<"  Total Energy  "<<"  Temperature  "<<"  Pressure  "<<"  Displacement  "<<endl;

    for(int state=0; state<numStates; state++){
        myfile << t_0*t[state] <<"      "<< epsilon*Ek[state] <<"  "<< epsilon*Ep[state]
           << "     "<< epsilon*Etot[state]<< "     "<<T_0*T[state]
           << "     "<< pressureFactor*epsilon/pow(sigma,3)*P[state]<<"  "<< pow(sigma,2)*D[state]<< endl;
    }
    myfile.close();
}


/************************************************************
Name:
Description:
*/
void FileManager::readDataFromFile()
{
    mat data;
    mat data_i = zeros(nLocalResAtoms,6);
    // Loading the first file
    stringstream file;

    int nRows= nLocalResAtoms;

    for(int s=1; s < stepLimit; s++){
        outName << statesDir <<"../STATES/state" << s << ".xyz";
        for(int p=0; p < nProc; p++){
            file.str("");
            file << statesDir << "s" << s << "p" << p << ".xyz";
            data_i.load( file.str() );
            data.insert_rows(p*nRows, data_i);
        }
        data.insert_rows(0,nLocalResAtoms*2);
        data.save(outName.str(),arma_ascii);
        outName.str( std::string() );
        outName.clear();
    }
}


/************************************************************
Name:           loadConfiguration
Description:    loadConfiguration
*/
void FileManager::loadConfiguration(Config* cfg){
    t_0=cfg->lookup("conversionFactors.t_0");
    sigma=cfg->lookup("conversionFactors.sigma");
    T_0= cfg->lookup("conversionFactors.T_0");
    epsilon= cfg->lookup("conversionFactors.epsilon");
    pressureFactor=cfg->lookup("conversionFactors.pressureFactor");
    stepLimit = cfg->lookup("statisticsSettings.numStates");

    cfg->lookupValue("fileManagerSettings.statisticsDir",statisticsDir);
    cfg->lookupValue("fileManagerSettings.statesDir",statesDir);

}
