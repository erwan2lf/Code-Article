#ifndef CONFIG_H
#define CONFIG_H

/*******************************************************************************************************
 * ⚙️  config.h
 * -------------------------------------------------------------------------------------------------
 * Physical, numerical parameters and output file names for a classic polymer chain simulation
 * (Brownian dynamics, harmonic springs + Lennard-Jones + spherical confinement).
 *
 * "Classic chain" version: no RNAP physics here (no transcription, no motor force).
 *******************************************************************************************************/

typedef struct {

    // === System size ===
    int N;                  // number of monomers in the chain
    int Nm;                 // number of monomers tracked for the legacy MSD arrays (legacy, see basic_functions.c)
    int total_num_simulations;     // number of simulations (legacy)

    // === Physical parameters ===
    double a;                // size of a monomer
    double Delta;            // time step
    double gamma_fric;       // friction coefficient
    double K;                // harmonic spring stiffness, chromatin-chromatin
    double K_bend;           // bending modulus
    double epsilon;          // Lennard-Jones potential strength, chromatin-chromatin
    double sigma, sigma6, sigma12;  // LJ sigma and its precomputed powers

    // === Spherical confinement ===
    double r_conf;           // radius of the confining sphere (0 = no confinement)
    double sigma_conf;       // range of the confinement potential
    double epsilon_conf;     // strength of the confinement potential

    // === Segment recorded in the binary trajectory (for MSD computation) ===
    int segment_start;
    int segment_end;

    // === Options ===
    int fixed_ends;           // 1 = chain endpoints fixed
    int confinement;          // 1 = spherical confinement active
    int plan;                 // 1 = reflective boundary at z<0
    int bending;               // 1 = bending forces active
    int temperature;          // 1 = thermal noise active
    int equilibrate;          // 1 = run an equilibration phase before measurement
    int resume_from_checkpoint; // 1 = resume from a checkpoint

    // === Random seed ===
    unsigned long seed;

    // === Durations and periods ===
    int T;                       // total number of time steps (measurement)
    int T_eq;                    // number of equilibration time steps
    int period_record;
    int period_msd;
    int period_correlation;
    int period_endtoend;
    int period_lammps;
    int period_center_of_mass;
    int period_neighbor;
    int period_force;

    int T_record;
    int T_msd;
    int T_correlation;
    int T_endtoend;
    int T_center_of_mass;
    int T_neighbor;
    int T_force;

    // === Output file names (in ./Results/) ===
    char file_name[256];
    char file_name_equilibrium[256];

    char file_name_test[256];
    char file_name_test2[256];
    char file_name_center_of_mass[256];

    char file_name_force[256];
    char file_name_force_thermal[256];
    char file_name_force_lj[256];

    char file_name_endtoend_segment[256];
    char file_name_endtoend_before[256];
    char file_name_endtoend_after[256];
    char file_name_endtoend[256];
    char file_name_neighbor[256];
    char file_name_correl_segment[256];

    char file_name_lammps[256];
    char file_name_equilibrium_lammps[256];

} Config;

Config parse_config(int argc, char *argv[]);
void print_header(const char *title);
void print_banner(void);

#endif // CONFIG_H
