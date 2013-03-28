#include "lj.h"

LJ::LJ():
    dr(zeros<vec>(3,1))
{
}


/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void LJ::calculateAndApplyForce(Atom *atomI, Atom *atomJ, int atomIsResident, int pairIsNotEvaluated){

    dr2=0.0;
    dRx = atomI->aPosition(0) - atomJ->aPosition(0);
    dRy = atomI->aPosition(1) - atomJ->aPosition(1);
    dRz = atomI->aPosition(2) - atomJ->aPosition(2);
    dr2= dRx*dRx + dRy*dRy + dRz*dRz;

    if (pairIsNotEvaluated && dr2 < rCut2) {
        dri2  = 1.0/dr2; dri6 = dri2*dri2*dri2; dr1 = sqrt(dr2);
        fcVal = 48.0*dri2*dri6*(dri6 - 0.5) + dUc / dr1;
        vVal  = 4.0*dri6*(dri6 - 1.0) - Uc - dUc*(dr1 - rCut);

        atomI->aAcceleration(0) += fcVal*dRx;
        atomI->aAcceleration(1) += fcVal*dRy;
        atomI->aAcceleration(2) += fcVal*dRz;

        if (atomIsResident){
            atomI->aPotential += vVal;

            atomJ->aAcceleration(0) -= fcVal*dRx;
            atomJ->aAcceleration(1) -= fcVal*dRy;
            atomJ->aAcceleration(2) -= fcVal*dRz;

        } else{
            atomI->aPotential += 0.5*vVal;
        }

        atomI->localPressure += fcVal*dr2;
        atomJ->localPressure += fcVal*dr2;
    }

}

/***************************************************************
 * Name:            loadConfiguration
 * Description:     Load system variables
 ***************************************************************/
void LJ::setParameters(Config* cfg){

    rCut = cfg->lookup("systemSettings.rCut");
    rCut2 = rCut*rCut;
    rCuti2 = 1.0/rCut2;
    rCuti6 = rCuti2*rCuti2*rCuti2;
    r1=sqrt(rCut2);
    Uc  =  4.0  * rCuti6 *(rCuti6 - 1.0);
    dUc = -48.0 * rCuti6 *(rCuti6 - 0.5) / r1;
}











