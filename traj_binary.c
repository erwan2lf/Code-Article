#include "traj_binary.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* ============================================================
 * Internal helpers
 * ============================================================ */

/* Writes an int32 in little-endian (portable). */
static void write_i32(FILE *f, int32_t v)
{
    fwrite(&v, sizeof(int32_t), 1, f);
}

/* Writes a float64 (double). */
static void write_f64(FILE *f, double v)
{
    uint8_t buf[8];
    memcpy(buf, &v, 8);
    fwrite(buf, 1, 8, f);
}

/* Writes a float32.
 * Positions don't need double precision for the MSD. */
static void write_f32(FILE *f, double v)
{
    float fv = (float)v;
    fwrite(&fv, sizeof(float), 1, f);
}


/* ============================================================
 * traj_open
 * ============================================================ */
FILE *traj_open(const char *path, const Config *cfg)
{
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, " traj_open: unable to open '%s'\n", path);
        return NULL;
    }

    int N_segment = cfg->segment_end - cfg->segment_start;

    /* --- Header --- */
    write_i32(f, (int32_t)TRAJ_MAGIC);
    write_i32(f, (int32_t)TRAJ_VERSION);
    write_i32(f, (int32_t)cfg->N);
    write_i32(f, (int32_t)N_segment);
    write_i32(f, (int32_t)cfg->segment_start);
    write_i32(f, (int32_t)cfg->segment_end);
    write_f64(f, cfg->Delta);
    write_f64(f, cfg->a);
    write_i32(f, (int32_t)cfg->period_record);

    fflush(f);

    printf(" traj_binary opened: %s\n", path);
    printf("   N_segment=%d  [%d -> %d]\n",
           N_segment,
           cfg->segment_start, cfg->segment_end);

    /* Estimate the final size */
    long bytes_per_frame =
        sizeof(int32_t)                                          /* timestep  */
        + (long)N_segment * 3 * sizeof(float);                   /* chromatin */

    long n_frames = cfg->T_record;
    double size_mb = (double)(bytes_per_frame * n_frames) / (1024.0 * 1024.0);

    printf("   ~%.1f MB estimated (%ld frames x %ld bytes/frame)\n\n",
           size_mb, n_frames, bytes_per_frame);

    return f;
}


/* ============================================================
 * traj_write_frame
 * ============================================================ */
void traj_write_frame(FILE *f, int timestep,
                      double **R,
                      const Config *cfg)
{
    if (!f) return;

    /* --- Timestep --- */
    write_i32(f, (int32_t)timestep);

    /* --- Chromatin segment [segment_start, segment_end) --- */
    for (int i = cfg->segment_start; i < cfg->segment_end; i++) {
        write_f32(f, R[i][0]);
        write_f32(f, R[i][1]);
        write_f32(f, R[i][2]);
    }
}


/* ============================================================
 * traj_close
 * ============================================================ */
void traj_close(FILE *f)
{
    if (f) {
        fflush(f);
        fclose(f);
        printf(" traj_binary closed.\n");
    }
}
