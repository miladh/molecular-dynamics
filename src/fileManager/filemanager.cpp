#include "filemanager.h"

FileManager::FileManager(const int &procID, const int &nProc):
    procID(procID),
    nProc(nProc),
    stepLimit(stepLimit)
{
}


/************************************************************
Name:           writeToFile
Description:    Writes the states to files
*/
void FileManager::writeAtomProperties(const int &state, const int &nLocalResAtoms,const vec &origo, Atom** atoms)
{
    outName << statesDir << "s" << state << "p" << procID << ".bin";
    myfile.open(outName.str().c_str(),ios::binary);

    for(int i=0;  i < nLocalResAtoms; i++){
        myfile <<join_rows(atoms[i]->aPosition+origo.t(), atoms[i]->aVelocity);
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}


/************************************************************
Name:
Description:
*/
void FileManager::readDataFromFile(Atom** atoms)
{
    mat data;
    stringstream file;

    file << rawDataDir << "s99p" << procID << ".bin";
    data.load( file.str() );

    nLocalResAtoms = data.n_rows;

    for(int p=0; p <nProc; p++ ){
        if(procID==p){
            for(int atom=0; atom < nLocalResAtoms; atom++ ){
                for(int i=0; i < 3; i++){
                    atoms[atom]->aPosition(i) = data(atom,i);
                    atoms[atom]->aVelocity(i) = data(atom,3+i);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/************************************************************
Name:
Description:
*/
void FileManager::writeSystemProperties(int numStates, const vec &t, const vec &Ek,
                                        const vec &Ep, vec const &Etot, const vec &T,
                                        const vec &P, const vec &D)
{

    outName << statisticsDir << "/statistics.bin";
    myfile.open (outName.str().c_str(),ios::binary);
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
    cfg->lookupValue("fileManagerSettings.rawDataDir",rawDataDir);

}



/************************************************************
    Name:
    Description:
    */
//void FileManager::readDataFromFile(Atom** atoms)
//{
//    mat data,data_i;
//    stringstream file;

//    for(int i=0;  i < nLocalResAtoms; i++){
//        for(int k=0; k <3 ; k++){
//            if(atoms[i]->frozen && atoms[i]->aVelocity(k)!=0){
//                cerr << "Atom not freezed!!  " << "procID: " << procID << endl;
//            }
//        }
//    }

//    file << rawDataDir << "s99p0" << ".bin";
//    data.load( file.str() );
//    cout << data.n_rows<<endl;
//    for(int p=1; p <nProc; p++){
//        file.str("");
//        file << rawDataDir << "s99p" << p << ".bin";
//        data_i.load( file.str() );
//        data = join_cols(data, data_i);
//    }

//    nLocalResAtoms = data.n_rows/nProc;
//    for(int atom=0; atom < nLocalResAtoms; atom++ ){
//        for(int i=0; i < 3; i++){
//            atoms[atom]->aPosition(i) = data(atom+procID*nLocalResAtoms,i);
//            atoms[atom]->aVelocity(i) = data(atom+procID*nLocalResAtoms,3+i);
//        }
//    }
//}
