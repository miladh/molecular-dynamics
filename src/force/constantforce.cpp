#include <src/force/constantforce.h>

ConstantForce::ConstantForce()
{
}


/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void ConstantForce::calculateAndApplyForce(Atom* atom, Atom *, int, int){

    for(int i = dir1 ; i <= dir2; i++){
        atom->aAcceleration(i) += forceMagnitude;
    }
    potEnergy += 0;
    pressure  += 0;
}

/***************************************************************
 * Name:            loadConfiguration
 * Description:     Load system variables
 ***************************************************************/
void ConstantForce::setParameters(Config* cfg){
    forceMagnitude = cfg->lookup("forceSettings.constant.forceMag");
    dir1 = cfg->lookup("forceSettings.constant.dir1");
    dir2 = cfg->lookup("forceSettings.constant.dir2");
}
