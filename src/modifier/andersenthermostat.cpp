#include <src/modifier/andersenthermostat.h>

AndersenThermostat::AndersenThermostat(System *system):
    Modifier(system)
{
}

void AndersenThermostat::apply(){

    for(int i=0; i<moleculeSystem->nLocalResAtoms; i++) {
        Atom** atoms= moleculeSystem->atoms;
        if(!atoms[i]->frozen && randu() < moleculeSystem->dt / collisionTime) {
            atoms[i]->aVelocity=randn<rowvec>(atoms[i]->nDimension)*sqrt(targetTemperature);
        }
    }

}

void AndersenThermostat::setCollisionTime(double collTime){
    collisionTime = collTime;
}


void AndersenThermostat::setTargetTemperature(double targetT) {
    targetTemperature = targetT;
}
