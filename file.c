#include "file.h"
#include"config.h"
#include"traj_binary.h"
#include <stdlib.h>

#include <sys/stat.h>   // mkdir
#include <sys/types.h>  // mode_t
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * 📂 open_simulation_files
 * --------------------------------------------------------------------------------------------------
 * Creates a "Results/" folder (if it doesn't already exist) and opens all the simulation's
 * output files inside that folder.
 *
 * Each file is opened in write mode ("w") in the directory:
 *    ./Results/
 * Example: ./Results/brownian_LJ.lammpstrj
 *
 * If a file fails to open, the function prints an explicit error and terminates the program.
 * --------------------------------------------------------------------------------------------------
 */
void open_simulation_files(const Config *cfg, Files *f)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  Creating the "Results" directory  ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    printf("\n--- DEBUG CONFIG open_simulation_files() ---\n");
    printf("cfg->file_name ptr = %p\n", (void*)cfg->file_name);
    printf("cfg->file_name = '%s'\n\n", cfg->file_name);

    const char *result_dir = "Results";

    struct stat st = {0};

    if (stat(result_dir, &st) == -1) {
        if (mkdir(result_dir, 0777) != 0) {
            perror("Error: unable to create the 'Results' folder");
            exit(EXIT_FAILURE);
        } else {
            printf("Folder 'Results' created successfully.\n");
        }
    } else {
        printf("Folder 'Results' already exists.\n");
    }



    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  ✏️ Open mode: "w" or "a"  ////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    // If resuming from a checkpoint -> APPEND to the existing files
    const char *mode = (cfg->resume_from_checkpoint ? "a" : "w");

    char path[512];

    // === Main files ===

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name);
    f->file = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_equilibrium);
    f->file_equilibrium = fopen(path, mode);

    // === Tests and monitoring ===
    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_test);
    f->test2 = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_test2);
    f->test = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_center_of_mass);
    f->center_of_mass = fopen(path, mode);

    // === Forces ===
    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_force);
    f->file_force = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_force_thermal);
    f->file_force_thermal = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_force_lj);
    f->file_force_lj = fopen(path, mode);

    // === End-to-end, neighbors, correlations ===
    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_endtoend_segment);
    f->file_endtoend_segment = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_endtoend_before);
    f->file_endtoend_before = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_endtoend_after);
    f->file_endtoend_after = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_endtoend);
    f->file_endtoend = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_neighbor);
    f->file_neighbor = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/%s", result_dir, cfg->file_name_correl_segment);
    f->file_correl_segment = fopen(path, mode);

    // === Parameters ===
    snprintf(path, sizeof(path), "%s/param.txt", result_dir);
    f->param = fopen(path, mode);

    snprintf(path, sizeof(path), "%s/trajectory.bin", result_dir);
    f->traj_bin = traj_open(path, cfg);



    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Open check  ////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

#define CHECK_FILE(ptr, name) \
    if (!(ptr)) { fprintf(stderr, "Error: unable to open '%s/%s'\n", result_dir, name); perror("fopen"); exit(EXIT_FAILURE); }

    CHECK_FILE(f->file, cfg->file_name);
    CHECK_FILE(f->file_equilibrium, cfg->file_name_equilibrium);
    CHECK_FILE(f->test, cfg->file_name_test2);
    CHECK_FILE(f->test2, cfg->file_name_test);
    CHECK_FILE(f->center_of_mass, cfg->file_name_center_of_mass);
    CHECK_FILE(f->file_force, cfg->file_name_force);
    CHECK_FILE(f->file_force_thermal, cfg->file_name_force_thermal);
    CHECK_FILE(f->file_force_lj, cfg->file_name_force_lj);
    CHECK_FILE(f->file_endtoend_segment, cfg->file_name_endtoend_segment);
    CHECK_FILE(f->file_endtoend_before, cfg->file_name_endtoend_before);
    CHECK_FILE(f->file_endtoend_after, cfg->file_name_endtoend_after);
    CHECK_FILE(f->file_endtoend, cfg->file_name_endtoend);
    CHECK_FILE(f->file_neighbor, cfg->file_name_neighbor);
    CHECK_FILE(f->file_correl_segment, cfg->file_name_correl_segment);
    CHECK_FILE(f->param, "param.txt");
    CHECK_FILE(f->traj_bin, "trajectory.bin");

#undef CHECK_FILE

    printf("All files were opened in ./Results/ (mode = %s)\n", mode);
}

void close_simulation_files(const Config *cfg, Files *f){

    fclose(f->file);
    fclose(f->file_equilibrium);

    fclose(f->test);
    fclose(f->test2);

    fclose(f->file_force);
    fclose(f->center_of_mass);
    fclose(f->file_force_thermal);
    fclose(f->file_force_lj);

    fclose(f->file_endtoend);
    fclose(f->file_endtoend_before);
    fclose(f->file_endtoend_segment);
    fclose(f->file_endtoend_after);
    fclose(f->file_neighbor);
    fclose(f->file_correl_segment);

    fclose(f->param);

    traj_close(f->traj_bin);

}

