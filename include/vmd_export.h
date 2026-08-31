#ifndef VMD_EXPORT_H
#define VMD_EXPORT_H

#include <stdio.h>

/*******************************************************************************************************
 * 🎥 vmd_export.h
 * -------------------------------------------------------------------------------------------------
 * Writes one frame in LAMMPS trajectory format (.lammpstrj), directly readable in VMD:
 *   File > New Molecule... > Browse (choose the file) > Determine file type: "LAMMPS Trajectory"
 *   (or leave it on automatic, the .lammpstrj extension is usually recognized).
 *
 * Unlike the old legacy_record_frame() in basic_functions.c:
 *   - REAL coordinates (no division by 1000, "x y z" columns instead of "xs ys zs" — the old
 *     code advertised "scaled" coordinates (xs/ys/zs, i.e. fractions of the box) while actually
 *     writing plain x/1000, which is not a box fraction at all: a source of confusion for a
 *     strict reader of the LAMMPS format)
 *   - bounding box recomputed at every frame (with a 10% margin) instead of a fixed, arbitrary
 *     ±1000 box, for correctly framed rendering in VMD.
 *
 * "Classic chain" version: a single atom type (1 = chromatin monomer), always cfg->N atoms
 * per frame. No RNAP here.
 *******************************************************************************************************/

void export_frame_vmd(FILE *file, double **R, int N, int t);

#endif // VMD_EXPORT_H
