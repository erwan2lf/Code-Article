#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <stdio.h>
#include "neighborlist.h"

/*******************************************************************************************************
 * ⚛️  potentials.h
 * -------------------------------------------------------------------------------------------------
 * Various potentials (Lennard-Jones, barrier, bead-bead spring, bending). Historical functions,
 * independent of RNAP, kept as-is.
 *******************************************************************************************************/

double ForceLJ(double sigma, double epsilon, double x, double d);

void potential_barrier(double** R, double* origin_point, double radius, double displacement_amount, double thickness, int N);

void force_bead_bead(double** R1, double** R2, double K_cohesin, double dist,
                       int particle1, int particle2, int cohesin_eq_distance, double dt);

void lennard_jones_forces(double **R, NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          double Delta, int fixed_ends,
                          int period_force_record,
                          FILE *file_force_lj, int t);

void f_bending_forces(double **R, double **t_link, double **bending_forces, double K_bend, int N, int t);

#endif // POTENTIALS_H
