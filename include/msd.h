#ifndef MSD_H
#define MSD_H

/*******************************************************************************************************
 * 📈 msd.h
 * -------------------------------------------------------------------------------------------------
 * Computes the mean squared displacement (MSD) per monomer from the binary trajectory,
 * via an O(T log T) algorithm based on the FFT (Calandrini / Kneller).
 *******************************************************************************************************/

int msd_compute_from_file(const char *bin_path, const char *out_path, int n_lags);

#endif // MSD_H
