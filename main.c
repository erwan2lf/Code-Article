#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "include/config.h"
#include "include/simulation.h"
#include "basic_functions.h"
#include "file.h"
#include "msd.h"
#include "mt19937ar.h"


#define DIM 3

double timespec_to_sec(struct timespec t) {
    return t.tv_sec + t.tv_nsec / 1e9;
}

int main(int argc, char*argv[])
{

    print_banner();
    setvbuf(stdout, NULL, _IONBF, 0);

    static SimVars sv = {0};
    static Files f = {0};

    fflush(stdout);

    // 1) Read config from the arguments
    Config cfg = {0};
    cfg = parse_config(argc, argv);

    int t_start = 0;

    // 2) Allocate buffers WITH the correct cfg dimensions
    init_sim_vars(&sv, &cfg);



    init_genrand(cfg.seed);
    printf(" Fresh start\n");


    // for(int i = 0; i < cfg.N; i++){
    //     printf("R[%d][0] = %lf R[%d][1] = %lf R[%d][2] = %lf \n", i, sv.R[i][0], i, sv.R[i][1], i, sv.R[i][2]);
    // }



    printf("First random number: %.10f\n", genrand_real2());


    printf("DEBUG period_record = %d\n", cfg.period_record);
    open_simulation_files(&cfg, &f);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////              Chromatin Creation           ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // create_straight_polymer(cfg.N, cfg.a, 2, sv.R);
    create_random_polymer(cfg.N, cfg.a, cfg.r_conf, sv.R);
    // create_fractal_globule(N, a, spacing, R);

    if (cfg.resume_from_checkpoint == 0)
    {
        // char file_path[512];

        // #ifdef CLUSTER
        //     snprintf(
        //         file_path,
        //         sizeof(file_path),
        //         "/home/elefloch/Simulation/Start_simu/N_1000/Simulations/simulation_seed_/brownian_LJ.lammpstrj",
        //         cfg.seed
        //     );
        // #else
        //     snprintf(
        //         file_path,
        //         sizeof(file_path),
        //         "/Users/erwan/Documents/These/Cluster/Start/simulation_seed_%lu/brownian_LJ.lammpstrj",
        //         cfg.seed
        //     );
        // #endif

        // printf(" Opening file: %s\n", file_path);

        // double** R_matrix = get_last_structure(file_path, cfg.N);
        // if (R_matrix == NULL)
        // {
        //     fprintf(stderr, "Error: Could not read the structure from the file.\n");
        //     return 1;
        // }
        // for (int i = 0; i < cfg.N; i++)
        // {
        //     for (int j = 0; j < 3; j++)
        //     {
        //         sv.R[i][j] = R_matrix[i][j];
        //     }
        // }
        // free_matrix_if_allocated(&R_matrix, cfg.N);

    }

    run_simulation(&cfg, &sv, &f, t_start);

    cleanup_sim_vars(&sv, &cfg);
    close_simulation_files(&cfg, &f);
    FILE *fin = fopen(".FINISHED", "w");
    if (fin) fprintf(fin, "DONE\n");
    fclose(fin);
    printf("=== SIMULATION TERMINATED ===\n");

    return 0;
}
