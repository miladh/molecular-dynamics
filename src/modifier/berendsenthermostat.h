#ifndef BERENDSENTHERMOSTAT_H
#define BERENDSENTHERMOSTAT_H

#include <src/modifier/modifier.h>

class BerendsenThermostat : public Modifier
{
public:
    BerendsenThermostat(System *system);

    void apply();

    void setTargetTemperature(double targetT);
    void setRelaxationTime(double relaxTime);

private:
    double currentTemperature,targetTemperature;
    double gamma, tau;
};

#endif // BERENDSENTHERMOSTAT_H
