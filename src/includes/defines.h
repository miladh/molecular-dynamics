#ifndef DEFINES_H
#define DEFINES_H

/* Constants------------------------------------------------------------
NMAX = Maximum # of atoms per processor
NEMAX = Maximum # of augmented (= resident + copied) atoms
NDBUF = Size of a double-precision buffer, dbuf
      > 6*(# of boundary atoms for each neighbor)
NBMAX = Maximum # of copied boundary atoms per neighbor.
NCLMAX = Maximum # of cells per processor.
RCUT = Potential cut-off length
MOVED_OUT: Signifies a moved-out resident atom in function atom_move.
EMPTY: Signifies the end of a linked list.
----------------------------------------------------------------------*/
#define NMAX 100000
#define NEMAX 200000
#define NDBUF 300000
#define NBMAX 100000
#define NCLMAX 100000
#define MOVED_OUT -1.0e10
#define EMPTY -1

#define PI 3.14159265359



enum force {
   noInteraction, constant, lennardJones, LJ_Constant
};

enum modifier {
    noModifier, Andersen, Berendsen
};

enum velocityDist {
    uniform, normal
};

enum poresShapes {
    cylindrical, circular
};

#endif // DEFINES_H
