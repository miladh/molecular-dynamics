#include <iostream>
#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <libconfig.h++>

#include <src/includes/defines.h>
#include <src/system/system.h>
#include <src/atom/atom.h>
#include <src/generator/generator.h>
#include <src/force/twobodyforce.h>
#include <src/force/lj.h>
#include <src/force/noforce.h>
#include <src/modifier/modifier.h>
#include <src/modifier/berendsenthermostat.h>
#include <src/modifier/andersenthermostat.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"



TwoBodyForce* setForceType(Config *cfg);
void setModifierType(Config *cfg , System* system);


int main(int argc, char* argv[])
{

    /********************************************************************
    Initializing
    */
    int procID, nProc;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    Config cfg;
    for(int i=0; i<nProc; i++){
        if(procID==i){
            if(argc <= 1){
                cfg.readFile("../md/src/config.cfg");
            }else{
                cfg.readFile(argv[1]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int nLocalResAtoms, nAtoms;
    string velocityDist= cfg.lookup("initialVelocitySetting.velocityDist");
    /***********************************************************************/


    Atom **ArAtoms= new Atom*[NEMAX];
    for(int i=0; i<NEMAX; i++) {
        ArAtoms[i] = new Atom(&cfg);
    }

    Generator ArGenerator(procID, nProc);
    ArGenerator.loadConfiguration(&cfg);
    ArGenerator.fccLatticeGenerator(ArAtoms);
    ArGenerator.setVelocity(ArAtoms,velocityDist);
    nLocalResAtoms = ArGenerator.getNLocalResAtoms();
    nAtoms = ArGenerator.getNAtoms();



    TwoBodyForce* force = setForceType(&cfg);
    force->setParameters(&cfg);



    System Ar(procID, nProc,  nLocalResAtoms,nAtoms ,ArAtoms);
    Ar.force = force;
    setModifierType(&cfg, &Ar);
    Ar.loadConfiguration(&cfg);
    Ar.MDrun();



    MPI_Finalize();
}


/************************************************************
Name:               setSolverMethod
Description:
*/
TwoBodyForce* setForceType(Config *cfg)
{
    TwoBodyForce* force;
    int forceType= cfg->lookup("forceSettings.forceType");

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
Name:               setSolverMethod
Description:
*/
void setModifierType(Config *cfg, System *system)
{
    double tau= cfg->lookup("ModifierSettings.tau");
    double targetTemperature=cfg->lookup("ModifierSettings.targetTemperature");
    double T_0 = cfg->lookup("conversionFactors.T_0");
    int modifierType= cfg->lookup("ModifierSettings.modifierType");

    switch (modifierType) {
    case noModifier:
        cout << "No modifiers" <<endl;
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






