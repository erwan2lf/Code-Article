#include "initial_structures.h"
#include "basic_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846


// Function to create a polymer
double** create_polymer(int N, double a, double spacing) {
    double** R = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        R[i] = (double*)malloc(3 * sizeof(double));}
    for (int i = 0; i < N; i++) {
        R[i][1] = spacing;
        R[i][2] = 0;}
    for (int i = 0; i < N; i++) {R[i][0] = -N * a / 2 + i * a;}

    return R; }

int find_max_timestep(const char* file_name) {
    FILE* file = fopen(file_name, "r");
    if (file == NULL) {
        printf("Error: Could not open file %s\n", file_name);
        return -1;
    }

    int max_timestep = -1;
    int stock = -1;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), file) != NULL) {

        if (strstr(buffer, "ITEM: TIMESTEP") != NULL) {
            int timestep;
            fscanf(file, "%d", &timestep);
            if (timestep > max_timestep) {
                stock = max_timestep;
                max_timestep = timestep;
            }
        }
    }

    fclose(file);
    printf("%d %d \n", max_timestep, stock);
    return stock;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Reads the last atomic structure from a LAMMPS file (.lammpstrj)
 * and returns a matrix [N][3] containing the positions (in nm).
 *
 * - Automatically handles the last frame of the file
 * - Ignores extra columns after zs
 * - Tolerates variable line formats
 */
double** get_last_structure(const char* file_name, int N) {
    FILE* file = fopen(file_name, "r");
    if (!file) {
        printf(" Error: unable to open file %s\n", file_name);
        return NULL;
    }

    // Allocate the position matrix
    double** R_matrix = malloc(N * sizeof(double*));
    if (!R_matrix) {
        printf(" Error: failed to allocate R_matrix\n");
        fclose(file);
        return NULL;
    }
    for (int i = 0; i < N; i++) {
        R_matrix[i] = malloc(3 * sizeof(double));
        if (!R_matrix[i]) {
            printf(" Error: failed to allocate R_matrix[%d]\n", i);
            for (int j = 0; j < i; j++) free(R_matrix[j]);
            free(R_matrix);
            fclose(file);
            return NULL;
        }
    }

    char buffer[512];
    long last_atoms_pos = -1; // position of the last "ITEM: ATOMS" block

    // Look for the last ATOMS section
    while (fgets(buffer, sizeof(buffer), file)) {
        if (strstr(buffer, "ITEM: ATOMS id type xs ys zs") != NULL) {
            // Save the position of the start of the block
            last_atoms_pos = ftell(file);
        }
    }

    if (last_atoms_pos == -1) {
        printf(" Error: no ATOMS section found in %s\n", file_name);
        for (int i = 0; i < N; i++) free(R_matrix[i]);
        free(R_matrix);
        fclose(file);
        return NULL;
    }

    // Rewind to the last ATOMS data block
    fseek(file, last_atoms_pos, SEEK_SET);

    // Read the N coordinate lines
    for (int i = 0; i < N; i++) {
        if (!fgets(buffer, sizeof(buffer), file)) {
            printf(" Error: unexpected end of file at particle %d\n", i);
            goto error_exit;
        }

        int id, type;
        double xs, ys, zs;
        char line_remainder[256]; // to ignore the extra columns

        // Read the first 5 columns, ignore the line_remainder
        int n = sscanf(buffer, "%d %d %lf %lf %lf %[^\n]", &id, &type, &xs, &ys, &zs, line_remainder);
        if (n < 5) {
            printf(" Error: invalid line at particle %d: %s\n", i, buffer);
            goto error_exit;
        }

        // Convert to nm
        R_matrix[i][0] = xs * 1000.0;
        R_matrix[i][1] = ys * 1000.0;
        R_matrix[i][2] = zs * 1000.0;
    }

    fclose(file);
    return R_matrix;

error_exit:
    for (int i = 0; i < N; i++) free(R_matrix[i]);
    free(R_matrix);
    fclose(file);
    return NULL;
}

void create_solenoid_polymer(int N, double a, double spacing, double thickness, double** R) {
    double end_to_end_distance = 30*20;// 30 beads of size 20
    double step_spacing = end_to_end_distance/N ;
    double th = 2 * PI / 300.0; // 10 points on the circle
    double theta = 3 * PI / 2.0;
    double radius = a / 2.0 / sin(th / 2.0);
    for (int i = 0; i < N; ++i) {
        R[i][0] = i * step_spacing + spacing;
        R[i][1] = radius * cos(theta);
        R[i][2] = radius * sin(theta);
        theta += th;
    }
    // Adjust Z for thickness and spacing
    double min_Z = R[0][2];
    for (int i = 1; i < N; ++i) {
        if (R[i][2] < min_Z) {
            min_Z = R[i][2];
        }
    }
    for (int i = 0; i < N; ++i) {
        R[i][2] = R[i][2] - min_Z - thickness + spacing;
    }

    // Center X
    double max_X = R[0][0];
    for (int i = 1; i < N; ++i) {
        if (R[i][0] > max_X) {
            max_X = R[i][0];
        }
    }
    for (int i = 0; i < N; ++i) {
        R[i][0] = R[i][0] - max_X / 2.0;
    }
}

void create_straight_polymer(int N, double a, double spacing, double** R){
    for ( int i = 0; i < N; i ++){
        R[i][0] = 0;
        R[i][1] = 0;
        R[i][2] = i * a + spacing;
    }
}

void create_fractal_globule(int N, double a, double spacing, double** R){
   for (int i = 0; i < N; i++){
    /*if (i<(int)N/3){
        R[i][0]=0;
        R[i][1]=0;
        R[i][2]=i*a+spacing;
    }
    if(i>=(int)N/3 && i<=(int)2*N/3){
        R[i][0]=0;
        R[i][1]=(i-N/3)*a+spacing;
        R[i][2]=(N/3)*a+spacing;
    }
    if ( i > (int)2*N/3){
        R[i][0]=0;
        R[i][1]=(N/3)*a+spacing;
        R[i][2]=-(i-N+1)*a+spacing;
    }*/
   if (i<(int)N/2){
        R[i][0]=0;
        R[i][1]=0;
        R[i][2]=i*a+spacing;
    }
    if ( i >= (int)N/2){
        R[i][0]=0;
        R[i][1]=a;
        R[i][2]=-(i-N+1)*a+spacing;
    }
   }

   }

void create_knot_structure(int N, double a, double **R){
    for(int i = 0; i < N; i++){
        if(i <= (N-24)/2){
            R[i][0] = 0; 
            R[i][1] = 0;
            R[i][2] = i*a;
        }
        
        if( i >= (N+24)/2){
            R[i][0] = 1; 
            R[i][1] = 0; 
            R[i][2] = -(i-N+1)*a;
        }
    }

    R[(N-24)/2+1][0] = 0;
    R[(N-24)/2+1][1] = 0;
    R[(N-24)/2+1][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+2][0] = 0;
    R[(N-24)/2+2][1] = 1;
    R[(N-24)/2+2][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+3][0] = 0;
    R[(N-24)/2+3][1] = 2;
    R[(N-24)/2+3][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+4][0] = 1;
    R[(N-24)/2+4][1] = 2;
    R[(N-24)/2+4][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+5][0] = 2;
    R[(N-24)/2+5][1] = 2;
    R[(N-24)/2+5][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+6][0] = 3;
    R[(N-24)/2+6][1] = 2;
    R[(N-24)/2+6][2] = R[(N-24)/2][2] +1;

    R[(N-24)/2+7][0] = 3;
    R[(N-24)/2+7][1] = 2;
    R[(N-24)/2+7][2] = R[(N-24)/2][2];

    R[(N-24)/2+8][0] = 3;
    R[(N-24)/2+8][1] = 2;
    R[(N-24)/2+8][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+9][0] = 3;
    R[(N-24)/2+9][1] = 1;
    R[(N-24)/2+9][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+10][0] = 2;
    R[(N-24)/2+10][1] = 1;
    R[(N-24)/2+10][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+11][0] = 1;
    R[(N-24)/2+11][1] = 1;
    R[(N-24)/2+11][2] = R[(N-24)/2][2] - 1;

    R[(N-24)/2+12][0] = 1;
    R[(N-24)/2+12][1] = 1;
    R[(N-24)/2+12][2] = R[(N-24)/2][2];

    R[(N-24)/2+13][0] = 1;
    R[(N-24)/2+13][1] = 1;
    R[(N-24)/2+13][2] = R[(N-24)/2][2] + 1;

    R[(N-24)/2+14][0] = 1;
    R[(N-24)/2+14][1] = 1;
    R[(N-24)/2+14][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+15][0] = 1;
    R[(N-24)/2+15][1] = 2;
    R[(N-24)/2+15][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+16][0] = 1;
    R[(N-24)/2+16][1] = 3;
    R[(N-24)/2+16][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+17][0] = 2;
    R[(N-24)/2+17][1] = 3;
    R[(N-24)/2+17][2] = R[(N-24)/2][2] + 2;

    R[(N-24)/2+18][0] = 2;
    R[(N-24)/2+18][1] = 3;
    R[(N-24)/2+18][2] = R[(N-24)/2][2] + 1;

    R[(N-24)/2+19][0] = 2;
    R[(N-24)/2+19][1] = 3;
    R[(N-24)/2+19][2] = R[(N-24)/2][2];

    R[(N-24)/2+20][0] = 2;
    R[(N-24)/2+20][1] = 2;
    R[(N-24)/2+20][2] = R[(N-24)/2][2];

    R[(N-24)/2+21][0] = 2;
    R[(N-24)/2+21][1] = 1;
    R[(N-24)/2+21][2] = R[(N-24)/2][2];

    R[(N-24)/2+22][0] = 2;
    R[(N-24)/2+22][1] = 0;
    R[(N-24)/2+22][2] = R[(N-24)/2][2];

    R[(N-24)/2+23][0] = 1;
    R[(N-24)/2+23][1] = 0;
    R[(N-24)/2+23][2] = R[(N-24)/2][2];

}
