#include "config.h"
#include "time.h"
#include "simulation.h"
#include "basic_functions.h"
#include "movement.h"
#include "forces.h"
#include "file.h"
#include "traj_binary.h"
#include "mt19937ar.h"
#include "neighborlist.h"
#include "vmd_export.h"


#include <stdlib.h>  // for malloc, calloc
#include <math.h>    // for log10()
#include <string.h>  // memset()

#define DEBUG_TIMING 0

typedef struct {
    const char *name;
    double wall_total;   // cumulated wall time
    double cpu_total;    // cumulated CPU time
    long calls;          // number of calls
} DebugTimer;

static double timespec_diff_s(struct timespec a, struct timespec b)
{
    return (double)(b.tv_sec - a.tv_sec) + 1e-9 * (double)(b.tv_nsec - a.tv_nsec);
}

static void timer_add(DebugTimer *tmr, struct timespec w0, struct timespec w1,
                      clock_t c0, clock_t c1)
{
    tmr->wall_total += timespec_diff_s(w0, w1);
    tmr->cpu_total  += (double)(c1 - c0) / CLOCKS_PER_SEC;
    tmr->calls++;
}

static void print_timer_summary(const DebugTimer *timers, int n, const char *title)
{
    printf("\n================ %s ================\n", title);
    printf("%-28s %12s %12s %10s %12s %12s\n",
           "Function", "Wall(s)", "CPU(s)", "Calls", "Wall/call", "CPU/call");

    for (int i = 0; i < n; i++) {
        double wall_avg = (timers[i].calls > 0) ? timers[i].wall_total / timers[i].calls : 0.0;
        double cpu_avg  = (timers[i].calls > 0) ? timers[i].cpu_total  / timers[i].calls : 0.0;

        printf("%-28s %12.6f %12.6f %10ld %12.6e %12.6e\n",
               timers[i].name,
               timers[i].wall_total,
               timers[i].cpu_total,
               timers[i].calls,
               wall_avg,
               cpu_avg);
    }
    printf("========================================================\n\n");
}


void init_sim_vars(SimVars *sv, Config *cfg) {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////                 Periods                    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    sv->T_record = cfg->T_record;
    sv->T_msd = cfg->T_msd;
    sv->T_correlation = cfg->T_correlation;
    sv->T_endtoend = cfg->T_endtoend;
    sv->T_center_of_mass = cfg->T_center_of_mass;
    sv->T_force = cfg->T_force;
    sv->T_neighbor = cfg->T_neighbor;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////          Array Initialization           ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    sv->bending_forces = allocate_matrix(cfg->N, 3);
    sv->R = allocate_matrix(cfg->N, 3);
    sv->R_new = allocate_matrix(cfg->N, 3);
    sv->t_link = allocate_matrix(cfg->N-1,3);
    sv->monomer_list = (int*)malloc(cfg->Nm * sizeof(int));
    for(int i = 0; i < cfg->Nm; i++) {
        sv->monomer_list[i] = i;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////           Correlations and Time            ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->correlation_history         = calloc(sv->T_correlation, sizeof(double));
    sv->correlation_segment_history = calloc(sv->T_correlation, sizeof(double));
    sv->time                      = malloc(sv->T_correlation * sizeof(double));
    sv->log10_time                = malloc(sv->T_correlation * sizeof(double));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////              End-to-end Storage              ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->stock         = allocate_matrix(sv->T_endtoend, 2);
    sv->R_endtoend_segment   = allocate_matrix(sv->T_endtoend, 3);
    sv->R_endtoend_before     = allocate_matrix(sv->T_endtoend, 3);
    sv->R_endtoend_after     = allocate_matrix(sv->T_endtoend, 3);
    sv->neighbor_count= malloc(cfg->N * sizeof(int));
    sv->R_endtoend           = allocate_matrix(sv->T_correlation, 3);
    sv->R_segment     = allocate_matrix(sv->T_correlation, 3);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////     Center of Mass & Gyration Radius    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->R_center_of_mass = allocate_matrix(sv->T_center_of_mass, 3);
    sv->com_history         = allocate_matrix(sv->T_center_of_mass, 3);
    sv->gyration_radius   = malloc(sv->T_center_of_mass * sizeof(double));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////                     MSD                    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->msd_storage = malloc(cfg->total_num_simulations * sizeof(double**));
    for (int i = 0; i < cfg->total_num_simulations; i++) {
        sv->msd_storage[i] = malloc(sv->T_msd * sizeof(double*));
        for (int j = 0; j < sv->T_msd; j++)
            sv->msd_storage[i][j] = calloc(cfg->Nm, sizeof(double));
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////         Monomer Trajectories         ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->monomer_position_arrays = malloc(sv->T_msd * sizeof(double**));
    for (int i = 0; i < cfg->T_msd; i++) {
        sv->monomer_position_arrays[i] = malloc(cfg->Nm * sizeof(double*));
        for (int j = 0; j < cfg->Nm; j++)
            sv->monomer_position_arrays[i][j] = calloc(3, sizeof(double));
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////      time/log10_time Initialization     ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < sv->T_correlation; i++) {
        sv->time[i]       = (double)i;
        sv->log10_time[i] = log10(sv->time[i] + 1e-12); // avoid log10(0)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////                Counters                   ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    sv->large_displacement_counter = 0;
    sv->test_variable_2 = 0;
    sv->test_variable_3 = 0;


    sv->com_conf = malloc(3 * sizeof(double));
    sv->num_large_moves = 0;


    // Chromatin forces
    sv->F_chrom = calloc(cfg->N, sizeof(*sv->F_chrom));
}

void cleanup_sim_vars(SimVars *sv, Config *cfg)
{
    // --- Simple matrices ---
    free_matrix_if_allocated(&sv->bending_forces, cfg->N);
    free_matrix_if_allocated(&sv->R, cfg->N);
    free_matrix_if_allocated(&sv->t_link, cfg->N - 1);
    free_matrix_if_allocated(&sv->R_endtoend, sv->T_correlation);
    free_matrix_if_allocated(&sv->R_segment, sv->T_correlation);
    free_matrix_if_allocated(&sv->R_endtoend_segment, sv->T_endtoend);
    free_matrix_if_allocated(&sv->R_endtoend_before, sv->T_endtoend);
    free_matrix_if_allocated(&sv->R_endtoend_after, sv->T_endtoend);
    free_matrix_if_allocated(&sv->R_center_of_mass, sv->T_center_of_mass);
    free_matrix_if_allocated(&sv->com_history, sv->T_center_of_mass);
    free_matrix_if_allocated(&sv->stock, sv->T_endtoend);

    // --- Simple arrays ---
    free_if_allocated((void**)&sv->neighbor_count);
    free_if_allocated((void**)&sv->gyration_radius);
    free_if_allocated((void**)&sv->correlation_history);
    free_if_allocated((void**)&sv->correlation_segment_history);
    free_if_allocated((void**)&sv->time);
    free_if_allocated((void**)&sv->log10_time);
    free_if_allocated((void**)&sv->monomer_list);

    // --- MSD ---
    free_matrix_cube_if_allocated(sv->msd_storage, cfg->total_num_simulations, sv->T_msd);

    // --- monomer_position_arrays ---
    if (sv->monomer_position_arrays)
    {
        for (int i = 0; i < sv->T_msd; i++)
        {
            if (sv->monomer_position_arrays[i])
            {
                for (int j = 0; j < cfg->Nm; j++)
                {
                    free_if_allocated((void**)&sv->monomer_position_arrays[i][j]);
                }
                free(sv->monomer_position_arrays[i]);
            }
        }
        free(sv->monomer_position_arrays);
    }

    // Forces
    free(sv->F_chrom);
}


void record_simulation_data(SimVars *sv, const Config *cfg, const Files *f, int t)
{

    if(t%cfg->period_lammps == 0){
        // Trajectory in LAMMPS format readable by VMD
        export_frame_vmd(f->file, sv->R, cfg->N, cfg->T_eq + t);
    }

    // Binary MSD — every period_record frames
    if (t % cfg->period_record == 0) {
        traj_write_frame(f->traj_bin, t, sv->R, cfg);
    }
}



#define DEBUG_TIMING   0
#define DEBUG_PROGRESS 1

#if DEBUG_TIMING
    #define TIMER_START(w0,c0) \
        clock_gettime(CLOCK_MONOTONIC, &(w0)); \
        (c0) = clock()

    #define TIMER_END(w0,w1,c0,c1,tmr) \
        clock_gettime(CLOCK_MONOTONIC, &(w1)); \
        (c1) = clock(); \
        timer_add(&(tmr), (w0), (w1), (c0), (c1))
#else
    #define TIMER_START(w0,c0) do {} while(0)
    #define TIMER_END(w0,w1,c0,c1,tmr) do {} while(0)
#endif

void run_calculation(SimVars *sv, const Config *cfg, const Files *f,
            NeighborList *neighbor_lists,
            int t_start)
{
    // ------------------------------------------------------------------
    // Global timers
    // ------------------------------------------------------------------
    struct timespec chrono_start, last, now;
    struct timespec loop_w0, loop_w1, sec_w0, sec_w1;
    clock_t loop_c0, loop_c1, sec_c0, sec_c1;

    double interval = 60.0;              // debug display every 60 s
    double checkpoint_limit_h = 47.0;    // 47 h
    double checkpoint_limit_s = checkpoint_limit_h * 3600.0;

    double total_loop_wall = 0.0;
    double total_loop_cpu  = 0.0;

    int print_stride = (cfg->T >= 10) ? (cfg->T / 10) : 1;

    // ------------------------------------------------------------------
    // Per-function / per-section timers
    // ------------------------------------------------------------------
    enum {
        TM_BUILD_NL = 0,
        TM_CONFINEMENT,
        TM_POLYMER_BM,
        TM_LJ,
        TM_LARGE_DISP_COUNT,
        TM_RECORD,
        TM_COMPUTE_MEASUREMENTS,
        TM_SAVE_CHECKPOINT,
        TM_COUNT
    };

    DebugTimer timers[TM_COUNT] = {
        [TM_BUILD_NL]       = {"build_neighbor_list",        0, 0, 0},
        [TM_CONFINEMENT]    = {"confinement_sphere",         0, 0, 0},
        [TM_POLYMER_BM]     = {"polymer_brownian_motion",   0, 0, 0},
        [TM_LJ]             = {"lennard_jones_forces",       0, 0, 0},
        [TM_LARGE_DISP_COUNT]    = {"count_large_displacements",       0, 0, 0},
        [TM_RECORD] = {"record_simulation_data",        0, 0, 0},
        [TM_COMPUTE_MEASUREMENTS] = {"compute_measurements",             0, 0, 0},
        [TM_SAVE_CHECKPOINT]= {"save_checkpoint",            0, 0, 0}
    };

    // ------------------------------------------------------------------
    // Initialization
    // ------------------------------------------------------------------
    clock_gettime(CLOCK_MONOTONIC, &chrono_start);
    clock_gettime(CLOCK_MONOTONIC, &last);

    int time_depart = (t_start != 0) ? t_start : 0;
    printf("▶️  Simulation started at t=%d with a limit of %.1f h (%.0f s)\n",
           time_depart, checkpoint_limit_h, checkpoint_limit_s);
    fflush(stdout);

    for (int t = time_depart; t < cfg->T; t++) {

        // --------------------------------------------------------------
        // Start loop timer
        // --------------------------------------------------------------
        clock_gettime(CLOCK_MONOTONIC, &loop_w0);
        loop_c0 = clock();

        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = timespec_diff_s(last, now);

        // --------------------------------------------------------------
        // Periodic DEBUG display
        // --------------------------------------------------------------
        #if DEBUG_TIMING
                if (elapsed >= interval) {
                    double total_elapsed_dbg = timespec_diff_s(chrono_start, now);
                    printf("Iteration %d (%.1f s since start)\n", t, total_elapsed_dbg);
                    fflush(stdout);
                    last = now;

                    print_timer_summary(timers, TM_COUNT, "Partial summary");
                }
        #endif

        // --------------------------------------------------------------
        // Checkpoint limit check
        // --------------------------------------------------------------
        double total_elapsed = timespec_diff_s(chrono_start, now);

        if (total_elapsed >= checkpoint_limit_s) {
            printf("\n Limit of %.1f h reached (%.2f h elapsed)\n",
                   checkpoint_limit_h, total_elapsed / 3600.0);
            printf("   -> Automatic checkpoint save and clean exit.\n");

            TIMER_START(sec_w0, sec_c0);
            save_checkpoint(sv, cfg, t);
            TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_SAVE_CHECKPOINT]);

            printf("✅ Checkpoint saved at t=%d\n", t);
            fflush(stdout);

            clock_gettime(CLOCK_MONOTONIC, &loop_w1);
            loop_c1 = clock();
            total_loop_wall += timespec_diff_s(loop_w0, loop_w1);
            total_loop_cpu  += (double)(loop_c1 - loop_c0) / CLOCKS_PER_SEC;

            #if DEBUG_TIMING
                        print_timer_summary(timers, TM_COUNT, "Final summary");
                        printf("Cumulated loop time: wall = %.6f s | cpu = %.6f s\n",
                            total_loop_wall, total_loop_cpu);
            #endif
                        exit(0);
                    }

        // --------------------------------------------------------------
        // Neighbor lists
        // --------------------------------------------------------------
        if (t % 1000 == 0) {
            TIMER_START(sec_w0, sec_c0);
            build_neighbor_list(sv->R, neighbor_lists, cfg->N, 2, 0);
            TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_BUILD_NL]);

            for (int i = 0; i < cfg->N; i++) {
                fprintf(f->file_neighbor, "%d ", neighbor_lists[i].count);
            }
            fprintf(f->file_neighbor, "\n");
        }


        memset(sv->F_chrom, 0, cfg->N * sizeof(*sv->F_chrom)); // Reset forces to zero

        for (int i = 0; i < cfg->N; i++) {
            if (isnan(sv->F_chrom[i][0]) || isinf(sv->F_chrom[i][0])) {
                printf(" NaN/Inf force on F_chrom[%d] at t=%d\n", i, t);
                exit(1);
            }
        }
        // --------------------------------------------------------------
        // Polymer dynamics
        // --------------------------------------------------------------
        TIMER_START(sec_w0, sec_c0);
        accumulate_spring_forces(sv->R, sv->F_chrom, cfg->K, cfg->N);
        TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_POLYMER_BM]);

        // --------------------------------------------------------------
        // Lennard-Jones
        // --------------------------------------------------------------
        TIMER_START(sec_w0, sec_c0);
        accumulate_lj_forces(sv->R, sv->F_chrom, neighbor_lists, cfg->N,
                    0.001 * cfg->epsilon, cfg->sigma6, cfg->sigma12,
                    cfg->fixed_ends);
        TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_LJ]);

        // --------------------------------------------------------------
        // Confinement
        // --------------------------------------------------------------
        if (cfg->confinement == 1) {
            TIMER_START(sec_w0, sec_c0);
            accumulate_conf_forces(sv->R, sv->F_chrom, cfg->N,
                sv->com_conf, cfg->r_conf,
                cfg->sigma_conf, cfg->epsilon_conf);
            TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_CONFINEMENT]);
        }

        // --------------------------------------------------------------
        // Large displacement counter
        // --------------------------------------------------------------
        double **R_temp = allocate_matrix(cfg->N, 3);

        for (int i = 0; i<cfg->N; i++){
            R_temp[i][0] = sv->R[i][0];
            R_temp[i][1] = sv->R[i][1];
            R_temp[i][2] = sv->R[i][2];
        }

        // --------------------------------------------------------------
        // Update positions (Euler-Maruyama)
        // --------------------------------------------------------------
        euler_maruyama_update(sv->R, sv->F_chrom, cfg->N, cfg->Delta,
            cfg->temperature, cfg->fixed_ends, cfg->plan, cfg->gamma_fric);

        // --------------------------------------------------------------
        // Large displacement counter
        // --------------------------------------------------------------

        TIMER_START(sec_w0, sec_c0);
        count_large_displacements(cfg->N, cfg->T, sv->R, R_temp,
                                     &sv->large_displacement_counter);
        TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_LARGE_DISP_COUNT]);
        free_matrix_if_allocated(&R_temp, cfg->N);


        // --------------------------------------------------------------
        // Miscellaneous
        // --------------------------------------------------------------

        TIMER_START(sec_w0, sec_c0);
        record_simulation_data(sv, cfg, f, t);
        TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_RECORD]);

        // --------------------------------------------------------------
        // End loop timer
        // --------------------------------------------------------------
        clock_gettime(CLOCK_MONOTONIC, &loop_w1);
        loop_c1 = clock();

        double loop_wall = timespec_diff_s(loop_w0, loop_w1);
        double loop_cpu  = (double)(loop_c1 - loop_c0) / CLOCKS_PER_SEC;

        total_loop_wall += loop_wall;
        total_loop_cpu  += loop_cpu;

        #if DEBUG_PROGRESS
                if (t % print_stride == 0) {
                    Measurements measurements = compute_measurements(sv->R, cfg->N);

        #if DEBUG_TIMING
                    TIMER_START(sec_w0, sec_c0);
                    // we don't want to time compute_measurements twice, so it can be omitted here
                    TIMER_END(sec_w0, sec_w1, sec_c0, sec_c1, timers[TM_COMPUTE_MEASUREMENTS]);
        #endif

            double duration_min = (int)(total_loop_wall / 60.0);
            double duration_sec = total_loop_wall - (duration_min * 60.0);
            double time_remaining = loop_wall * (cfg->T - t - 1) / 60.0;

            printf("%d/%d | wall_loop=%.6f s | cpu_loop=%.6f s | cum=%.0f:%.2f | ETA=%.2f min | std=%.10f | mean=%.10f\n",
                   t, cfg->T,
                   loop_wall, loop_cpu,
                   duration_min, duration_sec,
                   time_remaining,
                   measurements.std, measurements.mean);
            fflush(stdout);
        }
        #endif
    }

    // ------------------------------------------------------------------
    // Final summary
    // ------------------------------------------------------------------
#if DEBUG_TIMING
    print_timer_summary(timers, TM_COUNT, "Final summary");
#endif

    printf("Total time spent in the loop:\n");
    printf("  - wall      : %.6f s (%.3f h)\n", total_loop_wall, total_loop_wall / 3600.0);
    printf("  - cpu       : %.6f s (%.3f h)\n", total_loop_cpu,  total_loop_cpu  / 3600.0);
    printf("Number of large displacements: %d Percentage of large displacements: %lf\n",sv->large_displacement_counter, (double)sv->large_displacement_counter/(cfg->N*cfg->T));
    fflush(stdout);
}

void f_equilibrate(SimVars *sv, const Config *cfg, const Files *f, NeighborList *neighbor_lists)
{
    printf("//// Before the simulation //// \n");

    clock_t start, end;
    double loop_duration, total_duration = 0, time_remaining;

    for (int t = 0; t < cfg->T_eq; t++){

        start = clock();

        if (t%1000==0){
            build_neighbor_list(sv->R, neighbor_lists, cfg->N, 2, 0);

            for(int i = 0; i < cfg->N; i++){
                fprintf(f->file_neighbor,"%d ",neighbor_lists[i].count);
            }
            fprintf(f->file_neighbor, "\n");
        }

        polymer_brownian_motion(
            sv->R, cfg->K, cfg->Delta, cfg->N, cfg->K_bend, sv->bending_forces,
            cfg->fixed_ends, cfg->plan, t, f->test, cfg->bending, sv->unused_legacy_param,
            cfg->T, f->file_force, cfg->period_force,
            f->file_force_thermal, cfg->temperature);

        //update_link_vectors(R, t_link, N);
        //f_bending_forces(R, t_link, bending_forces, K_bend, N, t);

        end = clock();
        loop_duration = (double)(end - start)/CLOCKS_PER_SEC;
        total_duration += loop_duration ;

        if (t%(cfg->T_eq/10) == 0){
            Measurements measurements = compute_measurements (sv->R, cfg->N);
            double duration_min = (int)(total_duration / 60);
            double duration_sec = total_duration - (duration_min * 60);
            time_remaining = loop_duration * (cfg->T-t-1) / 60;
            printf("%d/%.d %.f:%.f %.2fmin std %.10f mean %.10f\n",t, cfg->T_eq, duration_min, duration_sec, time_remaining, measurements.std, measurements.mean);
        }

        if(t%cfg->period_record == 0){
            export_frame_vmd(f->file_equilibrium, sv->R, cfg->N, t);
        }
    }
}



void save_checkpoint(SimVars *sv, const Config *cfg, int t)
{
    unsigned long state[624];
    int index;
    get_mt_state(state, &index);

    char current[256], backup[256];
    snprintf(current, sizeof(current), "./Results/checkpoint_last.dat");
    snprintf(backup, sizeof(backup), "./Results/checkpoint_prev.dat");

    rename(current, backup);

    FILE *f = fopen(current, "wb");
    if (!f) {
        perror(" fopen checkpoint");
        return;
    }

    // --- META DATA ---
    fwrite(&t,          sizeof(int), 1, f);

    // --- POLYMER ---
    for (int i = 0; i < cfg->N; i++)
        fwrite(sv->R[i], sizeof(double), 3, f);

    // --- RNG ---
    fwrite(state, sizeof(state), 1, f);
    fwrite(&index, sizeof(int), 1, f);

    fclose(f);

    printf(" Checkpoint saved at t = %d\n", t);
}



void confinement_sphere(const Config *cfg, SimVars *sv, int t)
{
    // --- CENTER OF MASS INITIALIZATION AT t = 0 ---
    if (t == 0)
    {
        sv->com_conf[0] = 0.0;
        sv->com_conf[1] = 0.0;
        sv->com_conf[2] = 0.0;

        for (int i = 0; i < cfg->N; i++)
        {
            sv->com_conf[0] += sv->R[i][0];
            sv->com_conf[1] += sv->R[i][1];
            sv->com_conf[2] += sv->R[i][2];
        }
        sv->com_conf[0] /= cfg->N;
        sv->com_conf[1] /= cfg->N;
        sv->com_conf[2] /= cfg->N;
    }

    // --- PARAMETERS ---
    double R = cfg->r_conf;
    double sigma = cfg->sigma_conf;
    double eps = cfg->epsilon_conf;

    double sigma6 = pow(sigma, 6);
    double sigma12 = sigma6 * sigma6;
    (void)eps; (void)sigma12;

    // displacement threshold (important)
    double threshold = 20 * cfg->a;

    // local counter
    int displacement_too_large = 0;


    // --- LOOP OVER MONOMERS ---
    for (int i = 0; i < cfg->N; i++)
    {
        double dx = sv->R[i][0] - sv->com_conf[0];
        double dy = sv->R[i][1] - sv->com_conf[1];
        double dz = sv->R[i][2] - sv->com_conf[2];

        double d = sqrt(dx*dx + dy*dy + dz*dz);
        double wall_dist = R - d;

        if (wall_dist < sigma)
        {
            double inv = sigma / (sigma - wall_dist);
            double inv6 = pow(inv, 6);
            double inv12 = inv6 * inv6;

            double delta = 0.02 * (2.0 * inv12 - inv6);

            // clamp to avoid blow-ups
            double delta_max = 1e-4 * cfg->a;    // 5% of a bead's size
            if (delta > delta_max)
                delta = delta_max;

            // direction toward the center
            double nx = -dx / d;
            double ny = -dy / d;
            double nz = -dz / d;

            // --- Compute the displacement ---
            double dx_move = delta * nx;
            double dy_move = delta * ny;
            double dz_move = delta * nz;

            double move_norm = sqrt(dx_move*dx_move + dy_move*dy_move + dz_move*dz_move);

            // --- Check the displacement ---
            if (move_norm > threshold)
                displacement_too_large++;

            // --- Apply the displacement ---
            sv->R[i][0] += dx_move;
            sv->R[i][1] += dy_move;
            sv->R[i][2] += dz_move;

            // Safety: clamp onto the sphere
            double d2 = sqrt(
                (sv->R[i][0] - sv->com_conf[0]) * (sv->R[i][0] - sv->com_conf[0]) +
                (sv->R[i][1] - sv->com_conf[1]) * (sv->R[i][1] - sv->com_conf[1]) +
                (sv->R[i][2] - sv->com_conf[2]) * (sv->R[i][2] - sv->com_conf[2])
            );

            if (d2 > R)
            {
                double k = R / d2;
                sv->R[i][0] = sv->com_conf[0] + k * (sv->R[i][0] - sv->com_conf[0]);
                sv->R[i][1] = sv->com_conf[1] + k * (sv->R[i][1] - sv->com_conf[1]);
                sv->R[i][2] = sv->com_conf[2] + k * (sv->R[i][2] - sv->com_conf[2]);
            }
        }
    }

    // --- Store the counter in sv ---
    sv->num_large_moves += displacement_too_large;

    // // You can print this periodically
    // if (t % 1000 == 0 && t > 0)
    //     printf("[CONF] Cumulated large corrections = %d\n", sv->num_large_moves);
}


/*******************************************************************************************************
 * 🧬 run_simulation
 * -----------------------------------------------------------------------------------------------------
 * • General description:
 *   Main function for the Brownian dynamics simulation of a classic polymer chain (chromatin)
 *   subject to Lennard-Jones (LJ) interactions and spherical confinement.
 *
 *   It is the central engine of the simulation:
 *     - Initializes the neighbor structures (NeighborList),
 *     - Equilibrates the system (relaxation phase),
 *     - Runs the main time loop of the dynamics,
 *     - Handles recording (positions, metrics),
 *     - Finalizes the simulation and frees memory resources.
 *
 * • Inputs:
 *   - cfg  : physical and numerical parameters of the simulation (Config struct)
 *   - sv   : dynamic state of the system (positions, forces, etc.)
 *   - f    : struct grouping the output files
 *   - t_start : resume time step (0 for a fresh start)
 *
 *******************************************************************************************************/
void run_simulation(const Config *cfg, SimVars *sv, const Files *f, int t_start)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////   Simulation Initialization   ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    // --- Neighbor list allocation ---
    NeighborList *neighbor_lists = malloc(cfg->N * sizeof(NeighborList));

    if (!neighbor_lists) {
        fprintf(stderr, " Allocation error: neighbor_lists\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < cfg->N; i++) {
        neighbor_lists[i].neighbors = malloc(10 * sizeof(int));
        neighbor_lists[i].capacity  = 10;
        neighbor_lists[i].count     = 0;
        if (!neighbor_lists[i].neighbors) {
            fprintf(stderr, " Allocation error: neighbor_lists[%d].neighbors\n", i);
            exit(EXIT_FAILURE);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////    Equilibration   //////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    if (cfg->resume_from_checkpoint == 0 && cfg->equilibrate)
    {
        f_equilibrate(sv, cfg, f, neighbor_lists);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////    Main Loop   //////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    run_calculation(sv, cfg, f, neighbor_lists, t_start);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////    Finalization and Exit   //////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < cfg->N; i++)
        free(neighbor_lists[i].neighbors);
    free(neighbor_lists);

    printf("🏁 Simulation completed successfully.\n\n");
}
