#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include "config.h"

/*******************************************************************************************************
 * 📂 file.h
 * -------------------------------------------------------------------------------------------------
 * Groups all the output file descriptors of the simulation (positions, forces, end-to-end,
 * neighbors, binary trajectory, ...). All RNAP files have been removed.
 *******************************************************************************************************/

typedef struct {
    FILE *file;              // main LAMMPS trajectory
    FILE *file_equilibrium;    // LAMMPS trajectory of the equilibration phase

    FILE *test;
    FILE *test2;
    FILE *center_of_mass;

    FILE *file_force;
    FILE *file_force_thermal;
    FILE *file_force_lj;

    FILE *file_endtoend_segment;
    FILE *file_endtoend_before;
    FILE *file_endtoend_after;
    FILE *file_endtoend;

    FILE *file_neighbor;
    FILE *file_correl_segment;

    FILE *param;

    FILE *traj_bin;             // binary trajectory (for MSD computation)
} Files;

void open_simulation_files(const Config *cfg, Files *f);
void close_simulation_files(const Config *cfg, Files *f);

#endif // FILE_H
