#ifndef INITIAL_STRUCTURES_H
#define INITIAL_STRUCTURES_H

/*******************************************************************************************************
 * 🧬 initial_structures.h
 * -------------------------------------------------------------------------------------------------
 * Generation of the chain's initial configuration (straight, random, solenoid, fractal globule,
 * knot) and re-reading of a last structure from a .lammpstrj file.
 *******************************************************************************************************/

double** create_polymer(int N, double a, double spacing);

int find_max_timestep(const char* file_name);

double** get_last_structure(const char* file_name, int N);

void create_solenoid_polymer(int N, double a, double spacing, double thickness, double** R);

void create_straight_polymer(int N, double a, double spacing, double** R);

void create_fractal_globule(int N, double a, double spacing, double** R);

void create_knot_structure(int N, double a, double **R);

#endif // INITIAL_STRUCTURES_H
