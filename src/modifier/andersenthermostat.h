#ifndef ANDERSENTHERMOSTAT_H
#define ANDERSENTHERMOSTAT_H


#include <src/modifier/modifier.h>

class AndersenThermostat : public Modifier
{
public:
    AndersenThermostat(System *system);
    void apply();

    void setTargetTemperature(double targetT);
    void setCollisionTime(double collTime);


private:
    double targetTemperature;
    double collisionTime;
};

#endif // ANDERSENTHERMOSTAT_H

