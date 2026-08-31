#ifndef FORCES_H
#define FORCES_H

#include "neighborlist.h"

/*******************************************************************************************************
 * 💪 forces.h
 * -------------------------------------------------------------------------------------------------
 * Force accumulation for the classic polymer chain: harmonic spring, Lennard-Jones,
 * spherical confinement, and Euler-Maruyama update (overdamped Brownian dynamics).
 *******************************************************************************************************/

void accumulate_spring_forces(double **R, double (*F)[3], double K, int N);

void accumulate_lj_forces(double **R, double (*F)[3],
                          NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          int fixed_ends);

void accumulate_conf_forces(double **R, double (*F)[3], int N,
                            const double cdm[3],
                            double r_conf, double sigma_conf,
                            double epsilon_conf);

void euler_maruyama_update(double **R, double (*F)[3], int N,
                           double Delta, int temperature,
                           int fixed_ends, int plan, double gamma_fric);

#endif // FORCES_H
