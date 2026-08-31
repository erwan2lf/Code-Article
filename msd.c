#include "msd.h"
#include "traj_binary.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <fftw3.h>

/* ============================================================
 * Internal structures
 * ============================================================ */

typedef struct {
    int32_t magic;
    int32_t version;
    int32_t N;
    int32_t N_segment;
    int32_t segment_start;
    int32_t segment_end;
    double  Delta;
    double  a;
    int32_t period;
} BinHeader;

/* ============================================================
 * Header reading
 * ============================================================ */
static int read_header(FILE *f, BinHeader *h)
{
    if (fread(&h->magic,         sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->version,       sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->N,             sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->N_segment,     sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->segment_start, sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->segment_end,   sizeof(int32_t), 1, f) != 1) return -1;
    if (fread(&h->Delta,         sizeof(double),  1, f) != 1) return -1;
    if (fread(&h->a,             sizeof(double),  1, f) != 1) return -1;
    if (fread(&h->period,        sizeof(int32_t), 1, f) != 1) return -1;

    if (h->magic != (int32_t)TRAJ_MAGIC) {
        fprintf(stderr, "msd: incorrect magic (0x%X expected, 0x%X read)\n",
                TRAJ_MAGIC, h->magic);
        return -1;
    }
    return 0;
}

/* ============================================================
 * Frame counting
 * ============================================================ */
static int count_frames(FILE *f, const BinHeader *h)
{
    long frame_bytes =
        sizeof(int32_t)
        + (long)h->N_segment * 3 * sizeof(float);

    long pos_start = ftell(f);
    fseek(f, 0, SEEK_END);
    long pos_end = ftell(f);
    fseek(f, pos_start, SEEK_SET);

    long data_bytes = pos_end - pos_start;
    int n_frames = (int)(data_bytes / frame_bytes);

    printf("Binary file: %d frames detected (%.1f MB)\n",
           n_frames, (double)(pos_end) / (1024.0 * 1024.0));

    return n_frames;
}

/* ============================================================
 * SEQUENTIAL read of all trajectories
 *
 * Reads the file ONLY ONCE from start to end — no fseek.
 * Returns traj[mono][frame][3] allocated by the function.
 *
 * RAM: nseg x n_frames x 3 x sizeof(float)
 *       For nseg=100, n_frames=100000: ~115 MB
 * ============================================================ */
static float ***read_all_trajs_sequential(FILE *f, const BinHeader *h,
                                          int n_frames)
{
    int  nseg      = h->N_segment;

    /* Allocate traj[mono][frame][3] */
    float ***traj = malloc((size_t)nseg * sizeof(float **));
    if (!traj) { perror("malloc traj"); return NULL; }

    for (int i = 0; i < nseg; i++) {
        traj[i] = malloc((size_t)n_frames * sizeof(float *));
        if (!traj[i]) { perror("malloc traj[i]"); goto error; }
        for (int fr = 0; fr < n_frames; fr++) {
            traj[i][fr] = malloc(3 * sizeof(float));
            if (!traj[i][fr]) { perror("malloc traj[i][fr]"); goto error; }
        }
    }

    /* Sequential read — ONE frame at a time, no fseek */
    for (int fr = 0; fr < n_frames; fr++) {

        /* Timestep — ignored */
        int32_t ts;
        if (fread(&ts, sizeof(int32_t), 1, f) != 1) {
            fprintf(stderr, " fread timestep frame %d\n", fr);
            goto error;
        }

        /* Chromatin positions — stored */
        for (int i = 0; i < nseg; i++) {
            if (fread(traj[i][fr], sizeof(float), 3, f) != 3) {
                fprintf(stderr, " fread chrom frame %d mono %d\n", fr, i);
                goto error;
            }
        }

        if (fr % 10000 == 0)
            printf("   reading frame %d/%d\n", fr, n_frames);
    }

    return traj;

error:
    for (int i = 0; i < nseg; i++) {
        if (!traj[i]) break;
        for (int fr = 0; fr < n_frames; fr++)
            free(traj[i][fr]);
        free(traj[i]);
    }
    free(traj);
    return NULL;
}

/* ============================================================
 * Freeing the trajectories
 * ============================================================ */
static void free_all_trajs(float ***traj, int nseg, int n_frames)
{
    if (!traj) return;
    for (int i = 0; i < nseg; i++) {
        if (!traj[i]) continue;
        for (int fr = 0; fr < n_frames; fr++)
            free(traj[i][fr]);
        free(traj[i]);
    }
    free(traj);
}

/* ============================================================
 * MSD computation of a 1D trajectory via FFT
 * (Calandrini / Kneller algorithm, O(T log T))
 * ============================================================ */
static void msd_1d_fft(const double *x, int n, int n_lags,
                       double *msd_out, long *count)
{
    int nfft = 1;
    while (nfft < 2 * n) nfft <<= 1;

    fftw_complex *in  = fftw_alloc_complex(nfft);
    fftw_complex *out = fftw_alloc_complex(nfft);

    fftw_plan plan_fwd = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d(nfft, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int i = 0;  i < n;    i++) { in[i][0] = x[i]; in[i][1] = 0.0; }
    for (int i = n; i < nfft; i++) { in[i][0] = 0.0;  in[i][1] = 0.0; }

    fftw_execute(plan_fwd);
    for (int i = 0; i < nfft; i++) {
        double re = out[i][0], im = out[i][1];
        out[i][0] = re*re + im*im;
        out[i][1] = 0.0;
    }
    fftw_execute(plan_inv);

    double *S2 = malloc((size_t)n * sizeof(double));
    for (int m = 0; m < n; m++)
        S2[m] = in[m][0] / (double)nfft;

    double S1 = 0.0;
    for (int i = 0; i < n; i++) S1 += 2.0 * x[i] * x[i];

    for (int m = 0; m < n_lags; m++) {
        if (m > 0) S1 -= x[m-1]*x[m-1] + x[n-m]*x[n-m];

        int denom = n - m;
        if (denom <= 0) break;

        double msd_m = (S1 - 2.0 * S2[m]) / (double)denom;
        msd_out[m] += msd_m;
        count[m]   += 1;
    }

    free(S2);
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_inv);
    fftw_free(in);
    fftw_free(out);
}

/* ============================================================
 * msd_compute_from_file — public entry point
 * ============================================================ */
int msd_compute_from_file(const char *bin_path,
                          const char *out_path,
                          int         n_lags)
{
    FILE *f = fopen(bin_path, "rb");
    if (!f) {
        fprintf(stderr, " msd: unable to open '%s'\n", bin_path);
        return -1;
    }

    BinHeader h;
    if (read_header(f, &h) != 0) { fclose(f); return -1; }

    printf(" Header read:\n");
    printf("   N=%d  segment=[%d,%d)  N_segment=%d\n",
           h.N, h.segment_start, h.segment_end, h.N_segment);
    printf("   Delta=%.2e  period=%d\n", h.Delta, h.period);

    int n_frames = count_frames(f, &h);
    if (n_frames <= 0) {
        fprintf(stderr, " msd: no frame in the file\n");
        fclose(f);
        return -1;
    }

    if (n_lags <= 0 || n_lags > n_frames)
        n_lags = n_frames;

    printf("   n_frames=%d  n_lags=%d\n", n_frames, n_lags);

    double dt_frame = (double)h.period;
    int    nseg     = h.N_segment;

    double ram_mb = (double)nseg * n_frames * 3 * sizeof(float) / (1024.0 * 1024.0);
    printf("   Trajectory RAM: %.1f MB\n\n", ram_mb);

    /* ── Sequential read — A SINGLE PASS over the file ─── */
    printf(" Sequential reading...\n");
    float ***traj = read_all_trajs_sequential(f, &h, n_frames);
    fclose(f);

    if (!traj) {
        fprintf(stderr, " msd: read failed\n");
        return -1;
    }
    printf(" Reading complete.\n\n");

    /* ── Accumulator allocation ────────────────────────────────── */
    double **msd_per_mono = malloc((size_t)nseg * sizeof(double *));
    long   **count_mono   = malloc((size_t)nseg * sizeof(long *));
    for (int i = 0; i < nseg; i++) {
        msd_per_mono[i] = calloc((size_t)n_lags, sizeof(double));
        count_mono[i]   = calloc((size_t)n_lags, sizeof(long));
        if (!msd_per_mono[i] || !count_mono[i]) {
            perror("calloc");
            free_all_trajs(traj, nseg, n_frames);
            return -1;
        }
    }

    /* ── MSD computation ──────────────────────────────────────────────── */
    printf(" Computing MSD...\n");
    double *x = malloc((size_t)n_frames * sizeof(double));
    if (!x) { perror("malloc x"); free_all_trajs(traj, nseg, n_frames); return -1; }

    for (int i = 0; i < nseg; i++) {
        if (i % 10 == 0)
            printf("   monomer %d/%d\n", h.segment_start + i, h.segment_end);

        for (int dim = 0; dim < 3; dim++) {
            for (int fr = 0; fr < n_frames; fr++)
                x[fr] = (double)traj[i][fr][dim];
            msd_1d_fft(x, n_frames, n_lags, msd_per_mono[i], count_mono[i]);
        }
    }
    free(x);
    free_all_trajs(traj, nseg, n_frames);

    /* ── Mean and standard deviation ───────────────────────────────────── */
    double *msd_mean = calloc((size_t)n_lags, sizeof(double));
    double *msd_std  = calloc((size_t)n_lags, sizeof(double));

    for (int m = 0; m < n_lags; m++) {
        for (int i = 0; i < nseg; i++)
            if (count_mono[i][m] > 0)
                msd_per_mono[i][m] /= (double)count_mono[i][m];

        double sum = 0.0;
        for (int i = 0; i < nseg; i++) sum += msd_per_mono[i][m];
        msd_mean[m] = sum / (double)nseg;

        double sum2 = 0.0;
        for (int i = 0; i < nseg; i++) {
            double d = msd_per_mono[i][m] - msd_mean[m];
            sum2 += d * d;
        }
        msd_std[m] = sqrt(sum2 / (double)nseg);
    }

    /* ── Writing ────────────────────────────────────────────────── */
    FILE *fout = fopen(out_path, "w");
    if (!fout) {
        fprintf(stderr, " msd: unable to open '%s'\n", out_path);
        return -1;
    }

    fprintf(fout, "# MSD — segment [%d, %d)\n", h.segment_start, h.segment_end);
    fprintf(fout, "# n_monomers=%d  n_frames=%d  Delta=%.6e  period=%d\n",
            nseg, n_frames, h.Delta, h.period);
    fprintf(fout, "# columns: lag_index  lag_time  msd_mean  msd_std\n");
    fprintf(fout, "# [adim]    [timestep]    [a^2]     [a^2]\n");

    for (int m = 0; m < n_lags; m++) {
        double lag_time = (double)m * dt_frame;
        fprintf(fout, "%d\t%.6e\t%.6e\t%.6e\n",
                m, lag_time, msd_mean[m], msd_std[m]);
    }

    fclose(fout);
    printf(" MSD written to '%s'\n", out_path);

    /* ── Freeing ──────────────────────────────────────────────── */
    for (int i = 0; i < nseg; i++) {
        free(msd_per_mono[i]);
        free(count_mono[i]);
    }
    free(msd_per_mono);
    free(count_mono);
    free(msd_mean);
    free(msd_std);

    return 0;
}
