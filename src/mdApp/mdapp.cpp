#include "mdapp.h"

#include <src/force/lj.h>
#include <src/force/noforce.h>
#include <src/modifier/berendsenthermostat.h>
#include <src/modifier/andersenthermostat.h>

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


    if(!loadState){
        ofstream myfile;
        outName << stateDir << "state10.xyz";
        myfile.open(outName.str().c_str());

    }else{
    Generator ArGenerator(procID, nProc);
    ArGenerator.loadConfiguration(cfg);
    ArGenerator.fccLatticeGenerator(ArAtoms);
    ArGenerator.setVelocity(ArAtoms);

    nLocalResAtoms = ArGenerator.getNLocalResAtoms();
    nAtoms = ArGenerator.getNAtoms();
}
    TwoBodyForce* force = setForceType();
    force->setParameters(cfg);

    System Ar(procID, nProc,  nLocalResAtoms, nAtoms ,ArAtoms);
    Ar.force = force;
    setModifierType(&Ar);
    Ar.loadConfiguration(cfg);
    Ar.simulateSystem();
}


/************************************************************
Name:
Description:
*/
TwoBodyForce* MDApp::setForceType()
{
    TwoBodyForce* force;
    switch (forceType) {
    case noInteraction:
        force = new NoForce();
        break;

    case lennardJones:
        force = new LJ();
        break;
    }
    return force;
}

/************************************************************
Name:
Description:
*/
void MDApp::setModifierType(System *system)
{
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
        break;

    case Berendsen:
        BerendsenThermostat* berendsenMod;
        berendsenMod = new BerendsenThermostat(system);
        berendsenMod->setTargetTemperature(targetTemperature/T_0);
        berendsenMod->setRelaxationTime(tau);
        system->addModifiers(berendsenMod);
        break;
    }
}


/***************************************************************
Name:            loadConfiguration
Description:     Load system variables
*/
void MDApp::loadConfiguration(Config* cfg){
    this->cfg=cfg;
    tau= cfg->lookup("ModifierSettings.tau");
    T_0 = cfg->lookup("conversionFactors.T_0");
    targetTemperature=cfg->lookup("ModifierSettings.targetTemperature");
    modifierType= cfg->lookup("ModifierSettings.modifierType");
    forceType= cfg->lookup("forceSettings.forceType");
    loadState = cfg->lookup("initialStateSettings.loadState");
    cfg->lookupValue("fileManagerSettings.statesDir",stateDir);
}




