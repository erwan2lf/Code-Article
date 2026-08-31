#ifndef SIMULATION_H
#define SIMULATION_H

#include "config.h"
#include "file.h"
#include "neighborlist.h"

/*******************************************************************************************************
 * 🧮 simulation.h
 * -------------------------------------------------------------------------------------------------
 * Dynamic state (SimVars) and functions of the simulation engine for a classic polymer chain:
 * allocation, main loop, equilibration, checkpointing, confinement.
 * No RNAP physics here (no rnap.h, no RnapState field).
 *******************************************************************************************************/

typedef struct {

    // --- Chain positions and forces ---
    double **R;             // current positions [N][3]
    double **R_new;         // updated positions [N][3] (legacy / equilibration)
    double **t_link;        // bond vectors [N-1][3] (bending)
    double **bending_forces;// bending forces [N][3]
    double (*F_chrom)[3];   // accumulated forces on the chain [N][3]

    // --- Confinement ---
    double *com_conf;       // center of mass used for confinement
    int num_large_moves;

    // --- Counters ---
    int large_displacement_counter;
    int test_variable_2;
    int test_variable_3;
    int unused_legacy_param;               // legacy parameter passed to polymer_brownian_motion

    // --- Correlations and time ---
    double *correlation_history;
    double *correlation_segment_history;
    double *time;
    double *log10_time;
    double **R_endtoend;
    double **R_segment;

    // --- End-to-end ---
    double **stock;
    double **R_endtoend_segment;
    double **R_endtoend_before;
    double **R_endtoend_after;
    int *neighbor_count;

    // --- Center of mass & gyration radius ---
    double **R_center_of_mass;
    double **com_history;
    double *gyration_radius;

    // --- MSD (legacy, in-memory arrays — see also msd.c / traj_binary.c) ---
    double ***msd_storage;
    double ***monomer_position_arrays;
    int *monomer_list;

    // --- Durations (copied from Config at initialization) ---
    int T_record;
    int T_msd;
    int T_correlation;
    int T_endtoend;
    int T_center_of_mass;
    int T_force;
    int T_neighbor;

} SimVars;

void init_sim_vars(SimVars *sv, Config *cfg);
void cleanup_sim_vars(SimVars *sv, Config *cfg);

void record_simulation_data(SimVars *sv, const Config *cfg, const Files *f, int t);

void run_calculation(SimVars *sv, const Config *cfg, const Files *f,
           NeighborList *neighbor_lists, int t_start);

void f_equilibrate(SimVars *sv, const Config *cfg, const Files *f,
                    NeighborList *neighbor_lists);

void save_checkpoint(SimVars *sv, const Config *cfg, int t);

void confinement_sphere(const Config *cfg, SimVars *sv, int t);

// Main simulation entry point (replaces the old simu_LJ_RNAP_erwan)
void run_simulation(const Config *cfg, SimVars *sv, const Files *f, int t_start);

#endif // SIMULATION_H
