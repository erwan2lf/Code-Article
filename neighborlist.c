#include "simulation.h"
#include "config.h"
#include "basic_functions.h"
#include "movement.h"
#include "neighborlist.h"
#include <assert.h>


/**
 * @brief Builds the neighbor list for each particle based on a cutoff radius.
 *
 * For every particle i, this function finds all other particles j such that
 * the squared distance r_ij² < (RCUT + SKIN)², and stores their indices
 * in neighbor_lists[i]. The list expands dynamically as needed.
 *
 * @param R                Positions of all particles [N][3]
 * @param neighbor_lists   Array of NeighborList structures, one per particle
 * @param N                Total number of particles
 * @param RCUT             Cutoff radius for neighbor search
 * @param SKIN             Additional buffer distance for Verlet list updates
 *
 * @note This function reallocates memory dynamically if capacity is exceeded.
 *       Make sure neighbor_lists[i].neighbors and capacity are properly initialized.
 */
void build_neighbor_list(double **R, NeighborList *neighbor_lists, int N, double RCUT, double SKIN) {

    double rcut_sq = (RCUT + SKIN) * (RCUT + SKIN);

    for (int i = 0; i < N; i++) {
        neighbor_lists[i].count = 0;

        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            double dx = R[i][0] - R[j][0];
            double dy = R[i][1] - R[j][1];
            double dz = R[i][2] - R[j][2];
            double dist_sq = dx * dx + dy * dy + dz * dz;

            if (dist_sq < rcut_sq) {
                if (neighbor_lists[i].count >= neighbor_lists[i].capacity) {
                    neighbor_lists[i].capacity *= 2;
                    int *tmp = realloc(neighbor_lists[i].neighbors, neighbor_lists[i].capacity * sizeof(int));
                    if (!tmp) {
                        fprintf(stderr, " Error: realloc failed for neighbor list of particle %d\n", i);
                        exit(EXIT_FAILURE);
                    }
                    neighbor_lists[i].neighbors = tmp;
                }
                neighbor_lists[i].neighbors[neighbor_lists[i].count++] = j;
            }
        }
    }
}
