#include <src/System/system.h>
#include <src/force/lj.h>

System::System(const int &procID, const int &nProc, const int &nLocalResAtoms, const int &nAtoms,Atom** atoms):
    atoms(atoms),
    procID(procID),
    nProc(nProc),
    nLocalResAtoms(nLocalResAtoms),
    nAtoms(nAtoms)
{
}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::simulateSystem()
{
    initilizeParameters();
    setTopology();
    atomCopy();
    computeAccel();

    FileManager fileManager(procID, nProc);
    fileManager.loadConfiguration(cfg);

    double cpu1;
    cpu1 = MPI_Wtime();


    if(procID==0){
        cout << "-----------------Run integrator-------------------" << endl;
    }

    for (int i=1; i <= stepLimit; i++) {
        state = i-1;
        if(procID==0){
            cout << i << endl;
        }

        evaluateSystemProperties();

        if(state%stepAvg == 1 || state > 1000){
            fileManager.writeAtomProperties(state, nLocalResAtoms, origo,atoms);
        }
        singleStep();

        if(i < actingStep){
        applyModifier();
        }
    }


    cpu = MPI_Wtime() - cpu1;
    if (procID == 0){
        fileManager.writeSystemProperties(stepLimit,time,kinEnergy,potEnergy,totEnergy,
                                          temperature,pressure,displacement,meanVelocity);
        cout << "Elapsed time: "        << cpu  << " s" << endl;
        cout << "Communication time: "  << comt << " s" << endl;
    }
}


/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::setTopology(){

    imat deltaVec = zeros<imat>(6,3);
    shiftVec = zeros<mat>(6,3);
    neigID   = zeros<ivec>(6);
    myParity = zeros<ivec>(3,1);

    deltaVec(0,0) = -1;
    deltaVec(1,0) =  1;
    deltaVec(2,1) = -1;
    deltaVec(3,1) =  1;
    deltaVec(4,2) = -1;
    deltaVec(5,2) =  1;

    ivec tmp = zeros<ivec>(3,1);
    // Set up neighbor tables
    for (int k = 0; k < 6; k++) {
        // Vector index of neighbour k
        for (int a = 0; a < 3; a++){
            tmp(a) = (IDVec(a) + deltaVec(k,a) + procVec(a)) % procVec(a);
        }
        // Scalar neighbor ID
        neigID(k) = tmp(0)*procVec(1)*procVec(2) + tmp(1)*procVec(2) + tmp(2);

        // Shift vector
        for (int a=0; a<3; a++){
            shiftVec(k,a) = subsysSize(a)*deltaVec(k,a);
        }
    }

    // Set up the node parity table, myParity
    for (int i=0; i<3; i++) {
        if (procVec(i) == 1)
            myParity(i) = 2;
        else if (IDVec(i)%2 == 0)
            myParity(i) = 0;
        else
            myParity(i) = 1;
    }

    //Write
    if(procID==0){
        cout << "-----------------Setting topology------------------" << endl;
        cout << "Processor parity: " << myParity.t();
    }

}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::initilizeParameters(){

    procVec     = zeros<ivec>(3,1);
    IDVec       = zeros<ivec>(3,1);
    systemSize  = zeros<ivec>(3,1); //Px, Py, Pz2
    subsysSize  = zeros<vec>(3,1);  //Lx,Ly,Lz
    cellSize    = zeros<vec>(3,1);  //rcx,rcy,rcz
    origo       = zeros<vec>(3,1);
    nLocalCells = zeros<ivec>(3,1);
    boundaryAtomList = zeros<imat>(6,NBMAX);

    time = zeros<vec>(stepLimit,1);
    kinEnergy    = zeros<vec>(stepLimit,1);
    potEnergy    = zeros<vec>(stepLimit,1);
    totEnergy    = zeros<vec>(stepLimit,1);
    temperature  = zeros<vec>(stepLimit,1);
    pressure     = zeros<vec>(stepLimit,1);
    displacement = zeros<vec>(stepLimit,1);
    meanVelocity = zeros<vec>(stepLimit,1);

    vdt          = zeros<rowvec>(1,3);


    // Set system size and number of processors in x,y,z direction
    systemSize << Nc/nX << Nc/nY << Nc/nZ ;
    procVec << nX << nY << nZ;
    if(procVec(0)*procVec(1)*procVec(2)!=nProc){
        cerr << "Number of processors doesn't match!" << endl;
        exit(1);
    }

    // Vector index of this processor
    IDVec(0) = procID / (procVec(1) * procVec(2));
    IDVec(1) = (procID / procVec(2)) % procVec(1);
    IDVec(2) = procID % procVec(2);

    density = 4/pow(latticeConstant/sigma,3);

    // Compute subsystem size (Lx,Ly,Lz)
    for (int i=0; i < 3; i++){  //density = N^3 * 4 / (Lx*Ly*Lz) (not Nc but N!!)
        subsysSize(i) = (double)systemSize(i) / pow(density/4.0, 1.0/3.0);
    }

    // Compute the # of cells (lc, rc)
    for (int i=0; i<3; i++) {
        nLocalCells(i) = int(subsysSize(i) / rCut);
        cellSize(i)    = subsysSize(i) / nLocalCells(i);
        origo(i)       = (double)IDVec(i) * subsysSize(i);
    }

    if(loadState){
        for(int atom=0;  atom < nLocalResAtoms; atom++){
            atoms[atom]->aPosition -= origo.t();
        }
    }

    volume = nAtoms/density;
    //Write
    if(procID==0){
        cout << "-------------Initilizing Parameters---------------" << endl;
        cout << "System size: "     << systemSize.t();
        cout << "Subsystem size: "  << subsysSize.t();
        cout << "Number of cells: " << nLocalCells.t();
        cout << "cellSize:        " << cellSize.t();
    }
}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::computeAccel()
{
    ivec nInteractionCells = zeros<ivec>(3,1);
    ivec cellVec = zeros<ivec>(3,1);
    ivec neigCellVec = zeros<ivec>(3,1);

    int Lyz, Lxyz;
    int atomI, atomJ;
    int neigCell;
    int atomIsResident, pairIsNotEvaluated;

    // Reset the potential & forces
    restForce();

    // Make a linked-cell list, lscl
    for (int i=0; i < 3; i++){
        nInteractionCells(i) = nLocalCells(i) + 2;
    }
    Lyz  = nInteractionCells(1)*nInteractionCells(2);
    Lxyz = nInteractionCells(0)*Lyz;

    // Reset the headers, head
    int cell;
    for (cell=0; cell < Lxyz; cell++){
        cellList[cell] = EMPTY;
    }

    // Scan atoms to construct headers, cellVec, & linked lists, lscl
    for (int i=0; i < nLocalResAtoms + nBounAtoms; i++) {
        for (int j = 0; j < 3; j++){
            cellVec(j) = floor(abs( atoms[i]->aPosition(j) + cellSize(j) ) / cellSize(j));
        }

        // Translate the vector cell index, cellVec, to a scalar cell index
        cell = cellVec(0)*Lyz + cellVec(1)*nInteractionCells(2) + cellVec(2);
        // Link to the previous occupant (or EMPTY if you're the 1st)
        atomList[i] = cellList[cell];
        // The last one goes to the header
        cellList[cell] = i;
    } // Endfor atom i


    // Calculate pair interaction-------------------------------
    // Scan local cells
    for (cellVec(0) = 1; cellVec(0) <= nLocalCells(0); (cellVec(0))++ ){
        for (cellVec(1) = 1; cellVec(1) <= nLocalCells(1); (cellVec(1))++ ){
            for (cellVec(2) = 1; cellVec(2) <= nLocalCells(2); (cellVec(2))++ ){

                // Calculate a scalar cell index
                cell = cellVec(0)*Lyz + cellVec(1)*nInteractionCells(2) + cellVec(2);
                // Skip this cell if empty
                if (cellList[cell] == EMPTY){
                    continue;
                }

                // Scan the neighbor cells (including itself) of cell
                for (neigCellVec(0)=cellVec(0)-1; neigCellVec[0]<=cellVec(0)+1; (neigCellVec(0))++){
                    for (neigCellVec(1)=cellVec(1)-1; neigCellVec(1)<=cellVec(1)+1; (neigCellVec(1))++){
                        for (neigCellVec(2)=cellVec(2)-1; neigCellVec(2)<=cellVec(2)+1; (neigCellVec(2))++) {

                            // Calculate the scalar cell index of the neighbor cell
                            neigCell = neigCellVec(0)*Lyz + neigCellVec(1)*nInteractionCells(2) + neigCellVec(2);
                            // Skip this neighbor cell if empty
                            if (cellList[neigCell] == EMPTY){
                                continue;
                            }

                            // Scan atom i in cell
                            atomI = cellList[cell];
                            while (atomI != EMPTY) {

                                // Scan atom j in neigCell
                                atomJ = cellList[neigCell];
                                while (atomJ != EMPTY) {

                                    // No calculation with itself
                                    if (atomJ != atomI) { //atoms[atomI]->frozen != atoms[atomJ]->frozen
                                        // Check if resident atom
                                        atomIsResident = (atomJ < nLocalResAtoms);
                                        pairIsNotEvaluated = (atomI < atomJ);

                                        apply2BForces(atomI, atomJ,atomIsResident,pairIsNotEvaluated);

                                    } // Endif not self interaction

                                    atomJ = atomList[atomJ];
                                } // Endwhile j not empty

                                atomI = atomList[atomI];
                            } // Endwhile i not empty

                        }
                    }
                } // Endfor neighbor cells, c1

            }
        }
    }// Endfor central cell, c
}


/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::restForce()
{
    for(Force* force: forces2B){
        force->restPotentialEnergy();
        force->restPressure();
    }

    for (int i=0; i < nLocalResAtoms; i++){
        for (int j=0; j<3; j++){
            atoms[i]->aAcceleration(j) = 0.0;
        }
    }
}



/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::halfKick()
{    
    // Calculate onebody forces-------------------------------
    for(int atom = 0; atom <nLocalResAtoms+nBounAtoms; atom++){
        if(!atoms[atom]->frozen){
        apply1BForces(atom);
        }
    }

    rowvec3 vect;
    vect << 0 << 0 << 0.1;
    for (int i=0; i < nLocalResAtoms; i++){
        if(!atoms[i]->frozen){
            atoms[i]->aVelocity+=dt/2*atoms[i]->aAcceleration;
        }
    }
}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::singleStep()
{
    halfKick();
    for (int i=0; i<nLocalResAtoms; i++){
        if(!atoms[i]->frozen){
            vdt = dt*atoms[i]->aVelocity;
            atoms[i]->aPosition += vdt;
            atoms[i]->aDisplacement += vdt;
        }
    }
    atomMove();
    atomCopy();
    computeAccel();
    halfKick();

}


/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::evaluateSystemProperties()
{
    localKinEnergy = 0.0;
    localPotEnergy = 0.0;
    localPressure  = 0.0;
    localDisplacement = 0.0;
    localMeanVelocity =0.0;


    for (int i=0; i<nLocalResAtoms; i++) {
        if(!atoms[i]->frozen){
        localKinEnergy   += dot(atoms[i]->aVelocity,atoms[i]->aVelocity);
        localMeanVelocity     += atoms[i]->aVelocity(2);
        localDisplacement+= dot(atoms[i]->aDisplacement, atoms[i]->aDisplacement);

        }
    }

    localKinEnergy *= 0.5;

    for(Force* force: forces2B){
        localPotEnergy  += force->getPotentialEnergy();
        localPressure   += force->getPressure();
    }

    MPI_Allreduce(&localKinEnergy,&kinEnergy[state],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&localPotEnergy,&potEnergy[state],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&localPressure, &pressure[state],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&localDisplacement, &displacement[state],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&localMeanVelocity, &meanVelocity[state],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    time[state] = state*dt;
    displacement[state] /= nAtoms;
    meanVelocity[state] /=nAtoms;
    totEnergy[state]     = kinEnergy[state] + potEnergy[state];
    temperature[state]   = kinEnergy[state]/nAtoms*2.0/3.0;
    pressure[state]     *= 1/(3*volume);
    pressure[state]     += density*temperature[state];
}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
void System::atomCopy()
{
    int nRecAtoms = 0; // # of "received" boundary atoms
    int neighbor, neighborID;
    int nAtomsToBeSend, nAtomsToBeRecv;
    double com1;


    // Main loop over x, y & z directions starts
    for (int i=0; i<3; i++) {
        // Make a boundary-atom list
        // Reset the # of to-be-copied atoms for lower&higher directions
        for (int j=0; j<2; j++){
            boundaryAtomList(2*i + j, 0) = 0;
        }

        // Scan all the residents & copies to identify boundary atoms
        for (int atom=0; atom < nLocalResAtoms + nRecAtoms; atom++) {
            for (int j=0; j<2; j++) {
                neighbor = 2*i + j; // Neighbor ID
                // Add an atom to the boundary-atom list
                if (atomIsBoundary(atoms[atom]->aPosition,neighbor)){
                    boundaryAtomList(neighbor, ++(boundaryAtomList(neighbor, 0))) = atom;
                }
            }
        }
        // Message passing
        com1 = MPI_Wtime(); // To calculate the communication time

        // Loop over the lower & higher directions
        for (int j=0; j<2; j++) {
            neighbor = 2*i + j;
            neighborID = neigID(neighbor); //Neighbor node ID


            // Send & receive the # of boundary atoms
            nAtomsToBeSend = boundaryAtomList(neighbor, 0); // # of atoms to be sent


            // Even node: send & recv
            if (myParity(i) == 0) {
                MPI_Send(&nAtomsToBeSend,1,MPI_INT,neighborID,10,MPI_COMM_WORLD);
                MPI_Recv(&nAtomsToBeRecv,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
            }
            // Odd node: recv & send
            else if (myParity(i) == 1) {
                MPI_Recv(&nAtomsToBeRecv,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
                MPI_Send(&nAtomsToBeSend,1,MPI_INT,neighborID,10,MPI_COMM_WORLD);
            }
            // Single layer: Exchange information with myself
            else{
                nAtomsToBeRecv = nAtomsToBeSend;
            }

            // Send & receive information on boundary atoms

            // Message buffering
            for (int atom=1; atom <= nAtomsToBeSend; atom++){
                dbuf[4*(atom-1)] = atoms[boundaryAtomList(neighbor, atom)]->frozen;
                for (int k=0; k<3; k++){ // Shift the coordinate origin
                    dbuf[4*(atom-1) + k+1] = atoms[boundaryAtomList(neighbor, atom)]->aPosition(k) - shiftVec(neighbor, k);
                }
            }


            //Even node: send & recv
            if (myParity(i) == 0) {
                MPI_Send(dbuf,4*nAtomsToBeSend,MPI_DOUBLE,neighborID,20,MPI_COMM_WORLD);
                MPI_Recv(dbufr,4*nAtomsToBeRecv,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
            }
            // Odd node: recv & send
            else if (myParity(i) == 1) {
                MPI_Recv(dbufr,4*nAtomsToBeRecv,MPI_DOUBLE,MPI_ANY_SOURCE,20,MPI_COMM_WORLD,&status);
                MPI_Send(dbuf, 4*nAtomsToBeSend,MPI_DOUBLE,neighborID,20,MPI_COMM_WORLD);
            }
            // Single layer: Exchange information with myself
            else{
                for (int atom=0; atom < 4*nAtomsToBeRecv; atom++){
                    dbufr[atom] = dbuf[atom];
                }
            }

            // Message storing
            for (int atom=0; atom < nAtomsToBeRecv; atom++){
                atoms[nLocalResAtoms + nRecAtoms + atom]->frozen =  dbufr[4*atom];
                for (int k=0; k<3; k++){
                    atoms[nLocalResAtoms + nRecAtoms + atom]->aPosition(k) = dbufr[4*atom+k+1];
                }
            }

            // Increment the # of received boundary atoms
            nRecAtoms +=nAtomsToBeRecv;

            // Internode synchronization
            MPI_Barrier(MPI_COMM_WORLD);

        } // Endfor lower & higher directions

        comt += MPI_Wtime()-com1; // Update communication time

    } // Endfor x, y & z directions


    // Update the # of received boundary atoms
    nBounAtoms = nRecAtoms;

}


/***************************************************************
 * Name:
 * Description:
***************************************************************/
void System::atomMove()
{
    imat moveOutAtomList = zeros<imat>(6,NBMAX);
    int nImAtoms = 0; //# of new immigrants
    int lowDirNeig, highDirNeig;
    int neighbor, neighborID;
    int nAtomsToBeSend, nAtomsToBeRecv;
    double com1;

    //Main loop over x, y & z directions starts
    for (int i=0; i < 3; i++) {

        // Make a moved-atom list
        // Scan all the residents & immigrants to list moved-out atoms
        for (int atom =0; atom < nLocalResAtoms + nImAtoms; atom++) {
            lowDirNeig  = 2*i  ; // Neighbor ID
            highDirNeig = 2*i + 1;

            // Register a to-be-copied atom
            if(atoms[atom]->aPosition(0) > MOVED_OUT) {
                // Move to the lower direction
                if(atomDidMove(atoms[atom]->aPosition , lowDirNeig)){
                    moveOutAtomList(lowDirNeig, ++(moveOutAtomList(lowDirNeig,0)))  = atom;
                    if(atoms[atom]->frozen){
                        cerr << "Atom is forzen, but has moved! " << endl;
                    }
                }
                // Move to the higher direction
                else if(atomDidMove(atoms[atom]->aPosition , highDirNeig)){
                    moveOutAtomList(highDirNeig, ++(moveOutAtomList(highDirNeig,0))) = atom;
                    if(atoms[atom]->frozen){
                        cerr << "Atom is forzen, but has moved! " << endl;
                    }
                }
            }
        }

        // Message passing with neighbor nodes
        com1 = MPI_Wtime();

        // Loop over the lower & higher directions
        for (int j=0; j < 2; j++){
            neighbor = 2*i + j;
            neighborID = neigID(neighbor); // Neighbor node ID


            // Send atom-number information
            nAtomsToBeSend = moveOutAtomList(neighbor , 0); // # of atoms to-be-sent

            // Even node: send & recv
            if (myParity(i) == 0) {
                MPI_Send(&nAtomsToBeSend,1,MPI_INT,neighborID,110,MPI_COMM_WORLD);
                MPI_Recv(&nAtomsToBeRecv,1,MPI_INT,MPI_ANY_SOURCE,110,MPI_COMM_WORLD,&status);
            }
            // Odd node: recv & send
            else if (myParity(i) == 1) {
                MPI_Recv(&nAtomsToBeRecv,1,MPI_INT,MPI_ANY_SOURCE,110,MPI_COMM_WORLD,&status);
                MPI_Send(&nAtomsToBeSend,1,MPI_INT,neighborID,110,MPI_COMM_WORLD);
            }
            //Single layer: Exchange information with myself
            else{
                nAtomsToBeRecv = nAtomsToBeSend;
            }

            // Send & receive information on boundary atoms
            //Message buffering
            for (int atom=1; atom <= nAtomsToBeSend; atom++){
                dbuf[7*(atom-1)] = atoms[moveOutAtomList(neighbor ,atom)]->frozen;
                for (int k=0; k<3; k++) {
                    dbuf[7*(atom-1) + k+1] = atoms[moveOutAtomList(neighbor ,atom)]->aPosition(k) - shiftVec(neighbor ,k);
                    dbuf[7*(atom-1) + k+3+1] = atoms[moveOutAtomList(neighbor ,atom)]->aVelocity(k);
                    atoms[moveOutAtomList(neighbor ,atom)]->aPosition(0) = MOVED_OUT;

                    // Mark the moved-out atom
                }
            }

            // Even node: send & recv, if not empty
            if (myParity(i) == 0) {
                MPI_Send(dbuf,7*nAtomsToBeSend,MPI_DOUBLE,neighborID,120,MPI_COMM_WORLD);
                MPI_Recv(dbufr,7*nAtomsToBeRecv,MPI_DOUBLE,MPI_ANY_SOURCE,120,MPI_COMM_WORLD,&status);
            }
            // Odd node: recv & send, if not empty
            else if (myParity(i) == 1) {
                MPI_Recv(dbufr,7*nAtomsToBeRecv,MPI_DOUBLE,MPI_ANY_SOURCE,120, MPI_COMM_WORLD,&status);
                MPI_Send(dbuf,7*nAtomsToBeSend,MPI_DOUBLE,neighborID,120,MPI_COMM_WORLD);
            }
            // Single layer: Exchange information with myself
            else{
                for (int atom=0; atom < 7*nAtomsToBeRecv; atom++){
                    dbufr[atom] = dbuf[atom];
                }
            }

            // Message storing
            for (int atom=0; atom < nAtomsToBeRecv; atom++){
                atoms[nLocalResAtoms + nImAtoms + atom]->frozen = dbufr[7*atom];
                for (int k=0; k<3; k++) {
                    atoms[nLocalResAtoms + nImAtoms + atom]->aPosition(k) = dbufr[7*atom  + k+1];
                    atoms[nLocalResAtoms + nImAtoms + atom]->aVelocity(k) = dbufr[7*atom  + k+1+3];
                }
            }

            // Increment the # of new immigrants
            nImAtoms += nAtomsToBeRecv;

            // Internode synchronization
            MPI_Barrier(MPI_COMM_WORLD);

        } // Endfor lower & higher directions, id

        comt=comt + MPI_Wtime() - com1;
    } // Endfor x, y & z directions, i


    // Compress resident arrays including new immigrants
    int ipt = 0;
    for (int i=0; i < nLocalResAtoms + nImAtoms; i++) {
        if (atoms[i]->aPosition(0) > MOVED_OUT) {
            for (int j=0; j<3; j++) {
                atoms[ipt]->aPosition(j) = atoms[i]->aPosition(j);
                atoms[ipt]->aVelocity(j) = atoms[i]->aVelocity(j);
            }
            atoms[ipt]->frozen = atoms[i]->frozen;
            ipt++;
        }
    }
    /* Update the compressed # of resident atoms */
    nLocalResAtoms = ipt;

}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
int System::atomDidMove(rowvec r, int neighborID) {
    int i,id;
    i  = int(neighborID / 2); // x(0)|y(1)|z(2) direction
    id = neighborID % 2; // Lower(0)|higher(1) direction

    if (id == 0)
        return r(i) < 0.0;
    else
        return subsysSize(i) < r(i);
}

/***************************************************************
 * Name:
 * Description:
 ***************************************************************/
int System::atomIsBoundary(rowvec r, int neighborID)
{
    int i,id;
    i  = int(neighborID / 2);  // x(0)|y(1)|z(2) direction
    id = neighborID % 2;       // Lower(0)|higher(1) direction

    if (id == 0){
        return r(i) < rCut;
    }else{
        return subsysSize(i) - rCut < r(i);
    }
}

/************************************************************
Name:           AddModifiers
Description:
*/
void System::addModifiers(Modifier* modifier){
    modifiers.push_back(modifier);

}

/************************************************************
Name:           computeDynamics
Description:    Computes the dynamics of the system.
*/
void System::applyModifier(){

    for(Modifier* modifier: modifiers){
        modifier->apply();
    }
}


/************************************************************
Name:           AddModifiers
Description:
*/
void System::add1BForces(Force* force){
    forces1B.push_back(force);
}

/************************************************************
Name:           AddModifiers
Description:
*/
void System::add2BForces(Force* force){
    forces2B.push_back(force);
}

/************************************************************
Name:           computeDynamics
Description:    Computes the dynamics of the system.
*/
void System::apply1BForces(int atomI)
{
    for(Force* force: forces1B){
        force->calculateAndApplyForce(atoms[atomI],0,0,0);
    }
}


/************************************************************
Name:           computeDynamics
Description:    Computes the dynamics of the system.
*/
void System::apply2BForces(int atomI, int atomJ, int atomIsResident, int pairIsNotEvaluated){

    for(Force* force: forces2B){
        force->calculateAndApplyForce(atoms[atomI], atoms[atomJ],atomIsResident,pairIsNotEvaluated);
    }
}
/************************************************************
Name:           computeDynamics
Description:    Computes the dynamics of the system.
*/
double System::getTemperature(){
    return temperature[state];
}


/***************************************************************
 * Name:            loadConfiguration
 * Description:     Load system variables
 ***************************************************************/
void System::loadConfiguration(Config* cfg){
    this->cfg=cfg;
    Nc = cfg->lookup("systemSettings.Nc");
    nX = cfg->lookup("systemSettings.nX");
    nY = cfg->lookup("systemSettings.nY");
    nZ = cfg->lookup("systemSettings.nZ");
    latticeConstant = cfg->lookup("systemSettings.latticeConstant");
    sigma=cfg->lookup("conversionFactors.sigma");
    rCut = cfg->lookup("systemSettings.rCut");
    T_0 =cfg->lookup("conversionFactors.T_0");
    dt = cfg->lookup("statisticsSettings.dt");
    stepLimit = cfg->lookup("statisticsSettings.numStates");
    stepAvg = cfg->lookup("statisticsSettings.stepAvg");
    loadState = cfg->lookup("fileManagerSettings.loadState");
    actingStep = cfg->lookup("ModifierSettings.actingStep");
    cfg->lookupValue("fileManagerSettings.statesDir",path);

}


