#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************************************
 * 📋 neighborlist.h
 * -------------------------------------------------------------------------------------------------
 * Neighbor list (Verlet list) for Lennard-Jones interactions between monomers.
 *******************************************************************************************************/

typedef struct {
    int *neighbors;
    int count;
    int capacity;
} NeighborList;

void build_neighbor_list(double **R, NeighborList *neighbor_lists, int N, double RCUT, double SKIN);

#endif // NEIGHBORLIST_H
