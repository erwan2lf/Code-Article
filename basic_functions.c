#include "basic_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include "potentials.h"
#include "mt19937ar.h"

#include <string.h>
#include <sys/stat.h> // Pour mkdir()
#include <errno.h> 
#include <omp.h>



#define FILENAME_SIZE 100

#define PI 3.14159265358979323846


double randn() {
    //unsigned long seed=97714454678; init_genrand(seed);
    static int hasSpare = 0;
    static double spare;
    if (hasSpare) {
        hasSpare = 0;
        return spare;
    }
    hasSpare = 1;
    double u, v, s;
    do {
        u = genrand_real2() * 2 - 1;
        v = genrand_real2() * 2 - 1;
        //u = (rand() / ((double) RAND_MAX)) * 2.0 - 1;
        //v = (rand() / ((double) RAND_MAX)) * 2.0 - 1;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

// Computes the distance between two points in 3D
double distance(double* p1, double* p2) {
    return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}

// Computes the squared distance between two points in 3D
double distance2(double* p1, double* p2) {
    return pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2);
}

double norm(double* vec) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}
double norm2(double* vec) {
    return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}


double** allocate_matrix(size_t rows, size_t cols) {
    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for matrix row\n");
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

bool is_allocated(void *ptr) {
    return ptr != NULL;
}

void free_if_allocated(void **ptr) {
    if (is_allocated(*ptr)) {
        free(*ptr);
        *ptr = NULL;
    }
}

void free_matrix_if_allocated(double ***matrix, int rows) {
    if (is_allocated(*matrix)) {
        for (int i = 0; i < rows; i++) {
            if (is_allocated((*matrix)[i])) {
                free((*matrix)[i]);
                (*matrix)[i] = NULL; // Prevent double free
            }
        }
        free(*matrix);
        *matrix = NULL; // Prevent double free
    }
}

void free_matrix_cube_if_allocated(double ***matrix, int rows, int cols) {
    if (is_allocated(matrix)) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (is_allocated(matrix[i][j])) {
                    free(matrix[i][j]);
                    matrix[i][j] = NULL; // Prevent double free
                }
            }
            free(matrix[i]);
            matrix[i] = NULL; // Prevent double free
        }
        free(matrix);
        matrix = NULL; // Prevent double free
    }
}


//typedef struct {double mean;double std;} Measurements;
Measurements compute_measurements (double ** R, int N){
    double distance_list[N-1]; double mean = 0.0 ; double std = 0.0;
    for (int part=0; part<N-1; part ++){
    double d = distance(R[part], R[part+1]);
    distance_list[part]= d;
    mean += d; }
    mean = mean / (N-1) ;
    for (int part=0; part<N-1; part ++){std += (mean - distance_list[part]) * (mean - distance_list[part]);}
    std = sqrt(std/(N-1));

    //double mesures[2] ; mesures [0] =mean ; mesures[1] = std ;
    //return mean, std ;return mesures ;
    Measurements result;
    result.mean = mean;
    result.std = std;
    //printf("%f %f \n", mean, result.mean);
    return result;
    }

/*Mesures_temp calcul_mesures_temp (double ** R, int N,  int t){
    
}*/

void legacy_record_frame(FILE* file, double** R_matrix, int N, int t) {
    double TT = 1e+3;
    double cdm[3] = {0, 0, 0};
    fprintf(file, "ITEM: TIMESTEP\n");
    fprintf(file, "%d\n", t);
    fprintf(file, "ITEM: NUMBER OF ATOMS\n");
    fprintf(file, "%d\n", N);
    fprintf(file, "ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(file, "%lf %lf\n", -TT, TT);
    fprintf(file, "%lf %lf\n", -TT, TT);
    fprintf(file, "%lf %lf\n", -TT, TT);
    fprintf(file, "ITEM: ATOMS id type xs ys zs\n");

    for(int i = 0; i < N; i++){
        cdm[0] += R_matrix[i][0];
        cdm[1] += R_matrix[i][1];
        cdm[2] += R_matrix[i][2];
    }
    cdm[0] /= N;
    cdm[1] /= N;
    cdm[2] /= N;

    for (int particle = 0; particle < N; particle++) {
        fprintf(file, "%d 1 %lf %lf %lf -2\n", particle, (R_matrix[particle][0]) / 1000, (R_matrix[particle][1]) / 1000, (R_matrix[particle][2]) / 1000);
    }
}

double dot_product(double *vec1, double *vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void plot_columns_from_file(const char *filename, const char *title, const char *xlabel, const char *ylabel) {
    // Open the file for reading
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s for reading.\n", filename);
        return;
    }

    // Create a temporary file to store the data for Gnuplot
    FILE *tempFile = fopen("temp_data.txt", "w");
    if (tempFile == NULL) {
        fprintf(stderr, "Error: Could not create temporary file.\n");
        fclose(file);
        return;
    }

    // Read data from the input file and write to the temporary file
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        fprintf(tempFile, "%s", line);
    }

    fclose(file);
    fclose(tempFile);

    // Open a pipe to Gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title '%s'\n", title);
        fprintf(gnuplotPipe, "set xlabel '%s'\n", xlabel);
        fprintf(gnuplotPipe, "set ylabel '%s'\n", ylabel);
        fprintf(gnuplotPipe, "plot for [col=2:*] 'temp_data.txt' using 1:col with linespoints title columnheader\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        fprintf(stderr, "Error: Could not open gnuplot pipe.\n");
    }
    // Delete the temporary file
    if (remove("temp_data.txt") != 0) {
        fprintf(stderr, "Error: Could not delete temporary file.\n");
    }
}

void calculate_msd(double ****monomer_position_arrays,
                   double ***msd_storage,
                   int Tp, int Nm,
                   int *monomer_list,
                   int sim_index, int total_num_simulations)
{
    const char *folder_name = "MSD";

    // Create the folder if needed
    if (mkdir(folder_name, 0777) == -1) {
        if (errno == EEXIST)
            printf("Folder '%s' already exists.\n", folder_name);
        else
            perror("Error while creating folder");
    } else {
        printf("Folder '%s' created successfully.\n", folder_name);
    }

    // Prepare the output file
    char file_path[256];
    snprintf(file_path, sizeof(file_path), "%s/msd_file_%d.txt", folder_name, sim_index);
    FILE *msd_file = fopen(file_path, "w");
    if (!msd_file) {
        fprintf(stderr, "Error: unable to open %s\n", file_path);
        return;
    }

    printf("[Simulation %d] Computing MSD with OpenMP (%d threads)\n",
           sim_index, omp_get_max_threads());

    // --- Main loop ---
    for (int i = 0; i < Tp; i++) {
        fprintf(msd_file, "%d ", i);

        double msd_total = 0.0;

        // Parallelized over monomers
        #pragma omp parallel for reduction(+:msd_total) schedule(static)
        for (int k = 0; k < Nm; k++) {

            double msd = 0.0;

            for (int j = 0; j < Tp - i; j++) {
                msd += distance2(
                    monomer_position_arrays[sim_index][i + j][k],
                    monomer_position_arrays[sim_index][j][k]);
            }

            msd /= (Tp - i);   // time mean

            // Thread-safe storage: each thread writes to its own slot [k]
            msd_storage[sim_index][i][k] = msd;

            #pragma omp critical  // sequential write to the file
            fprintf(msd_file, "%f ", msd);

            msd_total += msd;
        }

        // Average over all monomers
        fprintf(msd_file, "%f\n", msd_total / Nm);
    }

    fclose(msd_file);
    printf("✅ File saved: %s\n", file_path);
}

// Function to check if a number is already in the list
bool is_in_list(int *list, int size, int num) {
    for (int i = 0; i < size; i++) {
        if (list[i] == num) {
            return true;
        }
    }
    return false;
}



void calculate_msd_serial(double ****R, double ***msd_storage,
                          int Tp, int Nm, int sim_index)
{
    for (int i = 0; i < Tp; i++) {
        double msd_total = 0.0;
        for (int k = 0; k < Nm; k++) {
            double msd = 0.0;
            for (int j = 0; j < Tp - i; j++) {
                msd += distance2(R[sim_index][i + j][k], R[sim_index][j][k]);
            }
            msd /= (Tp - i);
            msd_storage[sim_index][i][k] = msd;
            msd_total += msd;
        }
        msd_storage[sim_index][i][0] = msd_total / Nm;
    }
}


void calculate_msd_parallel(double ****R, double ***msd_storage,
                            int Tp, int Nm, int sim_index)
{
    for (int i = 0; i < Tp; i++) {
        double msd_total = 0.0;

        #pragma omp parallel for reduction(+:msd_total) schedule(static)
        for (int k = 0; k < Nm; k++) {
            double msd = 0.0;
            for (int j = 0; j < Tp - i; j++) {
                msd += distance2(R[sim_index][i + j][k], R[sim_index][j][k]);
            }
            msd /= (Tp - i);
            msd_storage[sim_index][i][k] = msd;
            msd_total += msd;
        }

        msd_storage[sim_index][i][0] = msd_total / Nm;
    }
}


// Computes the sum of an array
double sum_array(double* data, int size) {
    double sum = 0.0;
    for (int i = 1; i < size; i++) {
        sum += data[i];
    }
    
    return sum;
}

double calculate_angle(double *vec1, double *vec2){
    double dot_product_value = dot_product(vec1, vec2);
    double norm_v1 = norm(vec1);
    double norm_v2 = norm(vec2);
    
    return acos(dot_product_value / (norm_v1 * norm_v2));
}

void update_link_vectors(double **r_new,double **t_link, int N){
     for (int i = 1; i < N; i++){
        for(int j = 0; j < 3; j++){
            t_link[i-1][j] = r_new[i][j] - r_new[i-1][j];
        }
    }
}

void count_large_displacements(int N, int T, double **R, double **r_new, int *compteur){
    for(int i = 0; i < N; i++){
        if(distance(R[i], r_new[i]) > 0.1){
            (*compteur)++;
        }
    } 
}

//////////// // Replaces data with its discrete Fourier transform if Isign = 1, and the inverse if Isign = -1 ////////////////
void four1(double* data, int n, int Isign){
    int nn, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
    if(n<2 || n&(n-1)) {
        fprintf(stderr, "Number of points n must be a power of 2 in four1\n");
        exit(EXIT_FAILURE);
    }
    nn = n << 1;
    j = 1;
    for(i = 1; i < nn; i += 2) {
        if(j > i) {
            tempr = data[j]; data[j] = data[i]; data[i] = tempr;
            tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
        }
        m = n;
        while(m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while(nn > mmax) {
        istep = mmax << 1;
        theta = Isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for(m = 1; m < mmax; m += 2) {
            for(i = m; i <= nn; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

void create_random_polymer(int N, double a, double r_conf, double **R)
{
    double r_init = (r_conf > 0) ? r_conf * 0.5 : 50.0;

    // Initial position at the center of the sphere (± r_init/2)
    R[0][0] = (genrand_real2() - 0.5) * r_init;
    R[0][1] = (genrand_real2() - 0.5) * r_init;
    R[0][2] = (genrand_real2() - 0.5) * r_init;

    for (int i = 1; i < N; i++) {
        int max_tries = 1000;
        int ok = 0;

        for (int t = 0; t < max_tries; t++) {
            double phi   = genrand_real2() * 2 * PI;
            double theta = genrand_real2() * PI;

            double x = R[i-1][0] + a * sin(theta) * cos(phi);
            double y = R[i-1][1] + a * sin(theta) * sin(phi);
            double z = R[i-1][2] + a * cos(theta);

            // Accept only if inside the sphere
            if (r_conf > 0) {
                double d2 = x*x + y*y + z*z;
                if (d2 > r_conf * r_conf) continue;
            }

            R[i][0] = x;
            R[i][1] = y;
            R[i][2] = z;
            ok = 1;
            break;
        }

        // If no solution found, fall back toward the center
        if (!ok) {
            double phi   = genrand_real2() * 2 * PI;
            double theta = genrand_real2() * PI;
            double r_back = r_conf * 0.3;
            R[i][0] = r_back * sin(theta) * cos(phi);
            R[i][1] = r_back * sin(theta) * sin(phi);
            R[i][2] = r_back * cos(theta);
        }
    }
}
