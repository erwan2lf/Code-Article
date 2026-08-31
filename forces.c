#include "forces.h"
#include "basic_functions.h"   /* randn() */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ================================================================
 * CHROMATIN — harmonic spring
 * ----------------------------------------------------------------
 * Parallelizable version (OpenMP): reformulated as a loop PER MONOMER
 * (gather) rather than PER BOND (scatter). Each iteration i only reads
 * R[i-1], R[i], R[i+1] and only writes to F[i]: no write is shared
 * between threads, so no race condition, no need for atomics or
 * private buffers to reduce.
 *
 * The old version (one iteration per bond i<->i+1, which increments
 * F[i] AND F[i+1]) cannot be parallelized as-is: two threads handling
 * consecutive bonds both write to the same F[i+1], causing a "lost
 * update" (race condition).
 *
 * The "if(N > SPRING_PARALLEL_THRESHOLD)" clause automatically disables
 * parallelization (sequential execution, no thread overhead) for small
 * chains: at N~1000, the cost of OpenMP fork-join at EVERY time step
 * far outweighs the gain, and in practice it slows the simulation
 * down. The threshold below only becomes useful for much longer
 * chains (N gtrsim a few tens of thousands).
 * ================================================================ */
#define SPRING_PARALLEL_THRESHOLD 5000

void accumulate_spring_forces(double **R, double (*F)[3],
                              double K, int N)
{
    const double eps2 = 1e-16;

    #pragma omp parallel for schedule(static) if(N > SPRING_PARALLEL_THRESHOLD)
    for (int i = 0; i < N; i++) {
        double fx = 0.0, fy = 0.0, fz = 0.0;

        // Contribution from the bond to the next monomer (i, i+1)
        if (i < N - 1) {
            double dx = R[i+1][0] - R[i][0];
            double dy = R[i+1][1] - R[i][1];
            double dz = R[i+1][2] - R[i][2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > eps2) {
                double r    = sqrt(r2);
                double coef = K * (r - 1.0) / r;
                fx += coef * dx;
                fy += coef * dy;
                fz += coef * dz;
            }
        }

        // Contribution from the bond to the previous monomer (i-1, i)
        // (reaction opposite to the force this bond exerts on i-1)
        if (i > 0) {
            double dx = R[i][0] - R[i-1][0];
            double dy = R[i][1] - R[i-1][1];
            double dz = R[i][2] - R[i-1][2];
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > eps2) {
                double r    = sqrt(r2);
                double coef = K * (r - 1.0) / r;
                fx -= coef * dx;
                fy -= coef * dy;
                fz -= coef * dz;
            }
        }

        F[i][0] += fx;
        F[i][1] += fy;
        F[i][2] += fz;
    }
}

/* ================================================================
 * CHROMATIN — Lennard-Jones chrom-chrom
 * ================================================================ */
void accumulate_lj_forces(double **R, double (*F)[3],
                          NeighborList *neighbor_lists, int N,
                          double epsilon, double sigma6, double sigma12,
                          int fixed_ends)
{
    const double fmax = 300.0;
    const double c12  = 48.0 * epsilon * sigma12;
    const double c6   = 24.0 * epsilon * sigma6;

    for (int i = 1; i < N - 1; i++) {
        double *Ri = R[i];

        for (int k = 0; k < neighbor_lists[i].count; k++) {
            int j = neighbor_lists[i].neighbors[k];
            if (j <= i) continue;
            double *Rj = R[j];

            double dx = Ri[0] - Rj[0];
            double dy = Ri[1] - Rj[1];
            double dz = Ri[2] - Rj[2];

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 == 0.0) continue;

            double inv_r2  = 1.0 / r2;
            double inv_r6  = inv_r2 * inv_r2 * inv_r2;
            double inv_r8  = inv_r6 * inv_r2;
            double inv_r14 = inv_r8 * inv_r6;

            double f = c12 * inv_r14 - c6 * inv_r8;
            if (f >  fmax) f =  fmax;
            if (f < -fmax) f = -fmax;

            double fx = f * dx;
            double fy = f * dy;
            double fz = f * dz;

            if (fixed_ends == 1) {
                if (i != 0 && i != N-1) {
                    F[i][0] += fx; F[i][1] += fy; F[i][2] += fz;
                }
                if (j != 0 && j != N-1) {
                    F[j][0] -= fx; F[j][1] -= fy; F[j][2] -= fz;
                }
            } else {
                F[i][0] += fx; F[i][1] += fy; F[i][2] += fz;
                F[j][0] -= fx; F[j][1] -= fy; F[j][2] -= fz;
            }
        }
    }
}

/* ================================================================
 * CHROMATIN — spherical confinement
 * ================================================================ */
void accumulate_conf_forces(double **R, double (*F)[3], int N,
                            const double cdm[3],
                            double r_conf, double sigma_conf,
                            double epsilon_conf)
{
    const double fmax = 300.0;

    for (int i = 0; i < N; i++) {
        double dx = R[i][0] - cdm[0];
        double dy = R[i][1] - cdm[1];
        double dz = R[i][2] - cdm[2];

        double d = sqrt(dx*dx + dy*dy + dz*dz);
        if (d < 1e-12) continue;

        double d_wall = r_conf - d;
        if (d_wall >= sigma_conf || d_wall <= 0.0) continue;

        double inv   = sigma_conf / d_wall;
        double inv6  = inv*inv*inv*inv*inv*inv;
        double inv12 = inv6 * inv6;

        double f_mag = epsilon_conf * (12.0 * inv12 - 6.0 * inv6) / d_wall;
        if (f_mag >  fmax) f_mag =  fmax;
        if (f_mag < -fmax) f_mag = -fmax;

        double nx = -dx / d;
        double ny = -dy / d;
        double nz = -dz / d;

        F[i][0] += f_mag * nx;
        F[i][1] += f_mag * ny;
        F[i][2] += f_mag * nz;
    }
}

/* ================================================================
 * MISE À JOUR EULER-MARUYAMA
 * ================================================================ */
void euler_maruyama_update(double **R, double (*F)[3], int N,
                           double Delta, int temperature,
                           int fixed_ends, int plan, double gamma_fric)
{

    const double sigma_th = temperature ? sqrt(2.0 * (Delta / gamma_fric)) : 0.0;

    for (int i = 0; i < N; i++) {

        if (fixed_ends == 1 && (i == 0 || i == N-1))
            continue;

        R[i][0] += (Delta / gamma_fric) * F[i][0];
        R[i][1] += (Delta / gamma_fric) * F[i][1];
        R[i][2] += (Delta / gamma_fric) * F[i][2];

        if (temperature) {
            R[i][0] += sigma_th * randn();
            R[i][1] += sigma_th * randn();
            R[i][2] += sigma_th * randn();
        }

        if (plan && i > 0 && i < N-1 && R[i][2] < 0.0)
            R[i][2] = -R[i][2];
    }
}
