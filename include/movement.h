#ifndef MOVEMENT_H
#define MOVEMENT_H

#include <stdio.h>

/*******************************************************************************************************
 * 🚶 movement.h
 * -------------------------------------------------------------------------------------------------
 * "All-in-one" Brownian dynamics (spring + bending + thermal noise) used during the equilibration
 * phase. ideal_gas_motion is a free diffusion utility (for testing / debugging).
 *******************************************************************************************************/

void polymer_brownian_motion(double **R, double K, double Delta, int N,
                              double K_bend, double **bending_forces, int fixed_ends, int plan,
                              int t, FILE *test, int bending, int unused_legacy_param, int T,
                              FILE *file_force, int period_force_record,
                              FILE *file_force_thermal, int temperature);

double** ideal_gas_motion(double **R, int N, double **r_new, double Delta, int plan, int fixed_ends);

#endif // MOVEMENT_H
