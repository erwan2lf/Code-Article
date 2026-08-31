#ifndef TRAJ_BINARY_H
#define TRAJ_BINARY_H

#include <stdio.h>
#include "config.h"

/*******************************************************************************************************
 * 💾 traj_binary.h
 * -------------------------------------------------------------------------------------------------
 * Compact binary trajectory format (segment [segment_start, segment_end) of the chain),
 * used by msd.c to compute the MSD via FFT. The format no longer contains any RNAP field.
 *
 * Header (in order):
 *   int32  magic
 *   int32  version
 *   int32  N
 *   int32  N_segment
 *   int32  segment_start
 *   int32  segment_end
 *   f64    Delta
 *   f64    a
 *   int32  period_record
 *
 * Then, for each frame:
 *   int32  timestep
 *   float32[N_segment][3]  positions
 *******************************************************************************************************/

#define TRAJ_MAGIC   0x54524A31   // "TRJ1"
#define TRAJ_VERSION 2            // v2 = RNAP-free format

FILE *traj_open(const char *path, const Config *cfg);
void traj_write_frame(FILE *f, int timestep, double **R, const Config *cfg);
void traj_close(FILE *f);

#endif // TRAJ_BINARY_H
