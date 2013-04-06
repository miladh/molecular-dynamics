#include "mdapp.h"

#include <src/force/lj.h>
#include <src/force/noforce.h>
#include <src/force/constantforce.h>
#include <src/modifier/berendsenthermostat.h>
#include <src/modifier/andersenthermostat.h>
#include <src/fileManager/filemanager.h>
#include <src/pores/cylindricalpores.h>
#include <src/pores/circularpores.h>

MDApp::MDApp(const int &procID, const int &nProc):
    procID(procID),
    nProc(nProc)
{
}


/************************************************************
Name:
Description:
*/
void MDApp::runMDApp()
{
    Atom **ArAtoms= new Atom*[NEMAX];
    for(int i=0; i<NEMAX; i++) {
        ArAtoms[i] = new Atom(cfg);
    }


    if(loadState){
        FileManager filemanager(procID,nProc);
        filemanager.loadConfiguration(cfg);
        filemanager.readDataFromFile(ArAtoms);

        nLocalResAtoms = filemanager.nLocalResAtoms;
        MPI_Allreduce(&nLocalResAtoms, &nAtoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(procID==0){
        cout << "----------------Initilize configuration-------------------" << endl;
        cout << "Number of atoms:       "<< nAtoms << endl;
        }

        if(makePores){
            Pores* pores = setPoresShape();
            pores->makePores(ArAtoms);
        }

    }
    else{
        Generator ArGenerator(procID, nProc);
        ArGenerator.loadConfiguration(cfg);
        ArGenerator.fccLatticeGenerator(ArAtoms);
        ArGenerator.setVelocity(ArAtoms);

        nLocalResAtoms = ArGenerator.getNLocalResAtoms();
        nAtoms = ArGenerator.getNAtoms();
    }

    System Ar(procID, nProc, nLocalResAtoms, nAtoms ,ArAtoms);

    setForceType(&Ar);
    setModifierType(&Ar);
    Ar.loadConfiguration(cfg);
    Ar.simulateSystem();
}


/************************************************************
Name:
Description:
*/
Pores* MDApp::setPoresShape()
{
    switch (poresShape) {
    case cylindrical:
        Cylindricalpores *Cylpores;
        Cylpores = new Cylindricalpores(procID,nProc,nLocalResAtoms);
        Cylpores->loadConfiguration(cfg);
        return Cylpores;
        break;

    case circular:
        Circularpores *Cirpores;
        Cirpores = new Circularpores(procID,nProc,nLocalResAtoms);
        Cirpores->loadConfiguration(cfg);
        return Cirpores;
        break;

    default:
        if(procID==0){
        cerr << "Unknown pores Shape!" << endl;
        }
        exit(0);

    }
}


/************************************************************
Name:
Description:
*/
void MDApp::setForceType(System *system)
{
    if(procID==0){
        cout << "-----------------Setting Forces------------------" <<endl;
    }

    Force* force;
    switch (forceType) {
    case noInteraction:
        force = new NoForce();
        if(procID==0){
            cout << "No force" <<endl;
        }
        break;

    case constant:
        force = new ConstantForce();
        force->setParameters(cfg);
        system->addForces(force);
        if(procID==0){
            cout << "Constant force" <<endl;
        }
        break;

    case lennardJones:
        force = new LJ();
        force->setParameters(cfg);
        system->addForces(force);
        if(procID==0){
            cout << "Lennard Jones" <<endl;
        }
        break;

    case LJ_Constant:
         force = new LJ();
         force->setParameters(cfg);
         system->addForces(force);
         force = new ConstantForce();
         force->setParameters(cfg);
         system->addForces(force);
         if(procID==0){
             cout << "Lennard Jones+Constant" <<endl;
         }
         break;

    default:
        if(procID==0){
        cerr << "Unknown force type!" << endl;
        }
        exit(0);
    }
}

/************************************************************
Name:
Description:
*/
void MDApp::setModifierType(System *system)
{
    if(procID==0){
        cout << "-----------------Setting Modifieres---------------" <<endl;
    }

    switch (modifierType) {
    case noModifier:
        if(procID==0){
            cout << "No modifiers" <<endl;
        }
        break;

    case Andersen:
        AndersenThermostat* andersenMod;
        andersenMod = new AndersenThermostat(system);
        andersenMod->setTargetTemperature(targetTemperature/T_0);
        andersenMod->setCollisionTime(tau);
        system->addModifiers(andersenMod);
        if(procID==0){
            cout << " Andersen Thermostat" <<endl;
        }
        break;

    case Berendsen:
        BerendsenThermostat* berendsenMod;
        berendsenMod = new BerendsenThermostat(system);
        berendsenMod->setTargetTemperature(targetTemperature/T_0);
        berendsenMod->setRelaxationTime(tau);
        system->addModifiers(berendsenMod);
        if(procID==0){
            cout << "Berendsen Thermostat" <<endl;
        }
        break;
    }
}


/***************************************************************
Name:            loadConfiguration
Description:     Load system variables
*/
void MDApp::loadConfiguration(Config* cfg){
    this->cfg = cfg;
    tau = cfg->lookup("ModifierSettings.tau");
    T_0 = cfg->lookup("conversionFactors.T_0");
    targetTemperature = cfg->lookup("ModifierSettings.targetTemperature");
    modifierType = cfg->lookup("ModifierSettings.modifierType");
    forceType = cfg->lookup("forceSettings.forceType");
    loadState = cfg->lookup("fileManagerSettings.loadState");
    makePores = cfg->lookup("PoresSetting.makePores");
    poresShape = cfg->lookup("PoresSetting.poresShape");
    cfg->lookupValue("fileManagerSettings.statesDir",stateDir);
}




