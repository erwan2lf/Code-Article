#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include <stdio.h>
#include <stdbool.h>
#include <time.h>

/*******************************************************************************************************
 * 🧰 basic_functions.h
 * -------------------------------------------------------------------------------------------------
 * Generic toolbox: vector algebra, matrix allocation/deallocation, statistics,
 * correlations, MSD (legacy), Gaussian random number generation, gnuplot plotting, etc.
 * Nothing here is specific to RNAP.
 *******************************************************************************************************/

typedef struct {
    double std;
    double mean;
} Measurements;

// --- Random / algebra ---
double randn(void);
double distance(double* p1, double* p2);
double distance2(double* p1, double* p2);
double norm(double* vec);
double norm2(double* vec);
double dot_product(double *vec1, double *vec2);
double calculate_angle(double *vec1, double *vec2);

// --- Allocation ---
double** allocate_matrix(size_t rows, size_t cols);
bool is_allocated(void *ptr);
void free_if_allocated(void **ptr);
void free_matrix_if_allocated(double ***matrix, int rows);
void free_matrix_cube_if_allocated(double ***matrix, int rows, int cols);
void free_matrix_4D(double ****matrix, int Nm, int T_msd);

// --- Initial structures ---
void create_random_polymer(int N, double a, double r_conf, double **R);

// --- Divers / mesures ---
void print_structure(double** R_matrix, int N);
int* particles_outside(double** R, double* origin_point, double radius, double thickness, int N);
Measurements compute_measurements(double ** R, int N);
void legacy_record_frame(FILE* file, double** R_matrix, int N, int t);
double calculate_mean_for_correlation(double **R, int i, int j);
void calculate_autocorrelation(double** R_endtoend, double *correl_history, int T_correlation, int N);
void plot_columns_from_file(const char *filename, const char *title, const char *xlabel, const char *ylabel);
void plot_autocorrelation(double *time, double *autocorrelation, int length, const char *filename,
                          const char *title, const char *xlabel, const char *ylabel);
void plot_msd(double **msd_storage, int Tp, int Nm, int *monomer_list,
             const char *base_filename, const char *base_log_filename);
void save_center_of_mass(double** R_center_of_mass, int N, int T_center_of_mass,
                          int sim_index, int period_center_of_mass);
bool is_in_list(int *list, int size, int num);
void generate_unique_random_numbers(int *list, int N, int count);
double sum_array(double* data, int size);
void linear_regression(double* x, double* y, int size, double* alpha, double* logA);
void print_time_remaining(int current_iteration, int total_iterations, time_t start_time);
void plot_cols(double **array, int rows, int cols, int column_index, const char *filename,
              const char *title, const char *xlabel, const char *ylabel);
void update_link_vectors(double **r_new, double **t_link, int N);
void count_large_displacements(int N, int T, double **R, double **r_new, int *compteur);
void create_histogram(double *data, int data_size, int *bins, int bin_count, double min_value, double max_value);
bool is_present(int *list, int size, int value);
void swap_two_particles(double **R, int i, int j);
void reorganize_positions(double **R, int N);
void four1(double* data, int n, int Isign);
void realft(double* vec_double, double* data, int Isign, int data_size);
void calculate_gyration_radius(double** R_center_of_mass, double* gyration_radius, double** R, int N,
                               int T_center_of_mass, double** com_history);
void save_gyration(double* gyration_radius, int sim_index, int T_center_of_mass, int period_center_of_mass);

#endif // BASIC_FUNCTIONS_H
