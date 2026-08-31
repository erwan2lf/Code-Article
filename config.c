#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"


Config parse_config(int argc, char *argv[])
{
    Config cfg = {0};



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////// Reading arguments (defined in .bash) ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if (argc < 6)
    {
        fprintf(stderr,
            "Usage : %s seed Delta gamma_fric r_conf T\n",
            argv[0]);
        exit(1);
    }

    cfg.seed            = strtoul(argv[1], NULL, 10);
    cfg.Delta           = atof(argv[2]);
    cfg.gamma_fric      = atof(argv[3]);
    cfg.r_conf          = atof(argv[4]);
    long T_arg          = (long)atof(argv[5]);  // atof also handles scientific notation (e.g. 1e7)

    // print_header("Reading arguments (defined in .bash)");

    // printf("Seed  : %lu \n", cfg.seed);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////            Default Parameters           ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cfg.total_num_simulations  = 1;
    cfg.N               = 1000;
    cfg.a               = 1;
    cfg.K               = 10.0;
    cfg.K_bend          = 0.0;
    cfg.epsilon         = 0.0024;

    cfg.sigma           = cfg.a / 1.112; // a / 2^(1/6)
    cfg.sigma6          = pow(cfg.sigma,6);
    cfg.sigma12         = pow(cfg.sigma,12);

    // Chain segment recorded in the binary trajectory (for the MSD)
    cfg.segment_start   = 300;
    cfg.segment_end     = 400;

    cfg.epsilon_conf = 0.0024;
    cfg.sigma_conf = cfg.a;

    print_header("Default Parameters");


    printf("Number of monomers: %d\n", cfg.N);
    printf("Monomer size: %lf\n", cfg.a);
    printf("Chromatin spring stiffness: %lf\n", cfg.K);
    printf("Bending modulus: %lf\n", cfg.K_bend);
    printf("Time step: %lf\n", cfg.Delta);
    printf("Chromatin-chromatin LJ strength: %lf\n",cfg.epsilon);
    printf("LJ sigma (in units of a): %lf\n", cfg.sigma);
    printf("Segment recorded for MSD: [%d, %d) \n", cfg.segment_start, cfg.segment_end);
    printf("Confinement sphere size (0 = no confinement): %lf\n",cfg.r_conf);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////               Fixed Options                ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    cfg.fixed_ends = 0; // fixed_ends
    cfg.confinement = 1; // confinement
    cfg.plan = 0; // plan
    cfg.bending = 0; // bending
    cfg.temperature = 1; // temperature
    cfg.equilibrate = 1; // Equilibrate the system before computation


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////           Durations and Periods           ///////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // print_header("Durations and Periods");

    // T is provided directly on the command line (number of time steps) and used AS-IS
    // (no rounding): the actual simulation duration is therefore exactly the one requested.
    if (T_arg <= 0) {
        fprintf(stderr, "Error: T must be a positive integer (number of time steps)\n");
        exit(1);
    }
    cfg.T = (int)T_arg;
    cfg.T_eq = cfg.T / 10;

    // N_rec = desired number of recordings over the duration T (capped to T for short runs,
    // otherwise period_record would exceed T and nothing would ever be recorded).
    int N_rec = 100000;
    if (N_rec > cfg.T) N_rec = cfg.T;
    if (N_rec < 1) N_rec = 1;

    int N_rec10 = N_rec / 10;
    if (N_rec10 < 1) N_rec10 = 1;

    int k = (cfg.T + N_rec - 1) / N_rec;   // >= 1, depends only on T (not on a fixed multiple)

    printf("T = %d\n", cfg.T);

    cfg.period_record = k;
    printf("record period = %d \n", cfg.period_record);
    cfg.period_msd = k;
    cfg.period_correlation = (cfg.T + N_rec - 1) / N_rec;
    cfg.period_endtoend = (cfg.T + N_rec10 - 1) / N_rec10;
    cfg.period_lammps = cfg.T / 1000;
    if (cfg.period_lammps <= 0) cfg.period_lammps = 1;
    printf("lammps period = %d \n", cfg.period_lammps);
    cfg.period_center_of_mass = cfg.T;
    cfg.period_neighbor = cfg.T;
    cfg.period_force = cfg.T;

    cfg.T_record = cfg.T / cfg.period_record;
    printf("T_record = %d \n", cfg.T_record);
    cfg.T_msd = cfg.T / cfg.period_msd;
    cfg.T_correlation = cfg.T / cfg.period_correlation;
    cfg.T_endtoend = cfg.T / cfg.period_endtoend;
    cfg.T_center_of_mass = cfg.T / cfg.period_center_of_mass;
    cfg.T_neighbor = cfg.T / cfg.period_neighbor;
    cfg.T_force = cfg.T / cfg.period_force;


    // === File names ===
    snprintf(cfg.file_name, sizeof(cfg.file_name), "brownian_LJ.lammpstrj");

    snprintf(cfg.file_name_equilibrium, sizeof(cfg.file_name_equilibrium), "brownian_LJ_equilibrium.lammpstrj");

    snprintf(cfg.file_name_test, sizeof(cfg.file_name_test), "test.txt");
    snprintf(cfg.file_name_test2, sizeof(cfg.file_name_test2), "test2.txt");
    snprintf(cfg.file_name_center_of_mass, sizeof(cfg.file_name_center_of_mass), "center_of_mass.txt");

    snprintf(cfg.file_name_force, sizeof(cfg.file_name_force), "file_force.txt");
    snprintf(cfg.file_name_force_thermal, sizeof(cfg.file_name_force_thermal), "file_force_thermal.txt");
    snprintf(cfg.file_name_force_lj, sizeof(cfg.file_name_force_lj), "file_force_lj.txt");

    snprintf(cfg.file_name_endtoend_segment, sizeof(cfg.file_name_endtoend_segment), "endtoend_segment.txt");
    snprintf(cfg.file_name_endtoend_before, sizeof(cfg.file_name_endtoend_before), "endtoend_before.txt");
    snprintf(cfg.file_name_endtoend_after, sizeof(cfg.file_name_endtoend_after), "endtoend_after.txt");
    snprintf(cfg.file_name_endtoend, sizeof(cfg.file_name_endtoend), "endtoend.txt");
    snprintf(cfg.file_name_neighbor, sizeof(cfg.file_name_neighbor), "neighbor.txt");
    snprintf(cfg.file_name_correl_segment, sizeof(cfg.file_name_correl_segment), "correl_segment.txt");

    // === LAMMPS files ===
    snprintf(cfg.file_name_lammps, sizeof(cfg.file_name_lammps), "brownian_LJ.lammpstrj");
    snprintf(cfg.file_name_equilibrium_lammps, sizeof(cfg.file_name_equilibrium_lammps), "brownian_LJ_equilibrium.lammpstrj");


    ///// Simulation restart

    cfg.resume_from_checkpoint = 0;


    return cfg;

}

// === ANSI color codes ===
#define C_RESET   "\033[0m"
#define C_CYAN    "\033[36m"
#define C_YELLOW  "\033[1;33m"
#define C_BOLD    "\033[1m"

// === Display function ===
void print_header(const char *title)
{
    const int total_width = 180;   // total line width
    const int inner_padding = 4;   // spaces around the title
    const char border_char = '/';  // border character

    int title_len = strlen(title);
    int total_inner = title_len + 2 * inner_padding;

    // If the title is too long to fit within the width
    if (total_inner >= total_width - 4) {
        printf(C_CYAN "%c%c %s %c%c\n" C_RESET,
               border_char, border_char, title, border_char, border_char);
        return;
    }

    int side_width = (total_width - total_inner) / 2;

    // Top line
    printf(C_CYAN);
    for (int i = 0; i < total_width; i++) putchar(border_char);
    printf(C_RESET "\n");

    // Middle line with centered, colored title
    printf(C_CYAN);
    for (int i = 0; i < side_width; i++) putchar(border_char);
    printf(C_RESET);

    printf(C_YELLOW "%*s%s%*s" C_RESET, inner_padding, "", title, inner_padding, "");

    printf(C_CYAN);
    for (int i = 0; i < side_width; i++) putchar(border_char);
    printf(C_RESET "\n");

    // Bottom line
    printf(C_CYAN);
    for (int i = 0; i < total_width; i++) putchar(border_char);
    printf(C_RESET "\n");
}

#define C_RESET   "\033[0m"
#define C_BOLD    "\033[1m"
#define C_BLUE    "\033[38;5;39m"
#define C_CYAN    "\033[36m"
#define C_YELLOW  "\033[1;33m"
#define C_WHITE   "\033[97m"
#define C_GRAY    "\033[90m"

void print_banner(void) {
    printf("\n");
    printf(C_CYAN "//////////////////////////////////////////////////////////////////////////////////////////////////////////\n");
    printf(C_CYAN "////" C_RESET C_BOLD C_YELLOW "               Brownian Dynamics of a Classic Polymer Chain                 " C_RESET C_CYAN "////\n");
    printf(C_CYAN "//////////////////////////////////////////////////////////////////////////////////////////////////////////\n" C_RESET);

    printf(C_WHITE "\n");
    printf(C_BOLD "Code developed at " C_BLUE "Laboratoire de Physique Theorique (LPT - CNRS, Toulouse)" C_RESET "\n");
    printf(C_BOLD "Author: " C_YELLOW "Erwan Le Floch" C_RESET "\n");
    printf(C_BOLD "Last updated: " C_GRAY __DATE__ " - " __TIME__ C_RESET "\n");
    printf(C_BOLD "Compiler: " C_GRAY "gcc / g++ with OpenMP and HPC optimizations" C_RESET "\n");
    printf(C_BOLD "Simulation: " C_WHITE "Classic polymer chain (coarse-grained Brownian dynamics)" C_RESET "\n");

    printf("\n" C_CYAN "//////////////////////////////////////////////////////////////////////////////////////////////////////////" C_RESET "\n\n");
}
