#include "noforce.h"

NoForce::NoForce()
{
}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void NoForce::calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated){
    atomI->aAcceleration = zeros<rowvec>(1,3);
    atomJ->aAcceleration = zeros<rowvec>(1,3);
}

/***************************************************************
 * Name:            loadConfiguration
 * Description:     Load system variables
 ***************************************************************/
void NoForce::setParameters(Config* cfg){
}
