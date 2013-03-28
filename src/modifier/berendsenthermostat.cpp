#include <src/modifier/berendsenthermostat.h>


BerendsenThermostat::BerendsenThermostat(System * system):
    Modifier(system)
{
}

void BerendsenThermostat::apply()
{
    currentTemperature=moleculeSystem->getTemperature();
    gamma = sqrt(1 + moleculeSystem->dt / tau * (targetTemperature / currentTemperature - 1));

    for(int i=0; i<moleculeSystem->nLocalResAtoms; i++) {
        Atom** atoms= moleculeSystem->atoms;
        atoms[i]->aVelocity*=gamma;
    }
}


void BerendsenThermostat::setTargetTemperature(double targetT) {
    targetTemperature = targetT;
}


void BerendsenThermostat::setRelaxationTime(double relaxTime) {
    tau= relaxTime;
}
