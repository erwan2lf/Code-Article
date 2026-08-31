#include "vmd_export.h"

#include <stdio.h>

void export_frame_vmd(FILE *file, double **R, int N, int t)
{
    if (!file) return;

    // --- Bounding box: chain, 10% margin ---
    double xmin = R[0][0], xmax = R[0][0];
    double ymin = R[0][1], ymax = R[0][1];
    double zmin = R[0][2], zmax = R[0][2];

    for (int i = 1; i < N; i++) {
        if (R[i][0] < xmin) xmin = R[i][0];
        if (R[i][0] > xmax) xmax = R[i][0];
        if (R[i][1] < ymin) ymin = R[i][1];
        if (R[i][1] > ymax) ymax = R[i][1];
        if (R[i][2] < zmin) zmin = R[i][2];
        if (R[i][2] > zmax) zmax = R[i][2];
    }

    double mx = (xmax - xmin) * 0.1 + 1.0;
    double my = (ymax - ymin) * 0.1 + 1.0;
    double mz = (zmax - zmin) * 0.1 + 1.0;
    xmin -= mx; xmax += mx;
    ymin -= my; ymax += my;
    zmin -= mz; zmax += mz;

    // --- LAMMPS header (REAL coordinates, "ff" = non-periodic box) ---
    fprintf(file, "ITEM: TIMESTEP\n%d\n", t);
    fprintf(file, "ITEM: NUMBER OF ATOMS\n%d\n", N);
    fprintf(file, "ITEM: BOX BOUNDS ff ff ff\n");
    fprintf(file, "%lf %lf\n", xmin, xmax);
    fprintf(file, "%lf %lf\n", ymin, ymax);
    fprintf(file, "%lf %lf\n", zmin, zmax);
    fprintf(file, "ITEM: ATOMS id type x y z\n");

    // --- Chromatin: type 1 ---
    for (int i = 0; i < N; i++) {
        fprintf(file, "%d 1 %lf %lf %lf\n", i, R[i][0], R[i][1], R[i][2]);
    }

    fflush(file);
}
