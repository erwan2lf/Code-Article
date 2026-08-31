"""
Reader and viewer for the binary trajectory format.

Computes the end-to-end vector between two beads i and j:
    r(t) = R[j](t) - R[i](t)
and plots |r(t)|^2 as well as its x, y, z components.

Binary format (classic chain, no RNAP):
    Header: magic(i32) version(i32) N(i32) N_segment(i32)
            segment_start(i32) segment_end(i32) Delta(f64) a(f64)
            period_record(i32)
    Per frame: timestep(i32) then N_segment*3 float32 (positions)
"""

import sys
import struct
import numpy as np
import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────
# Reading
# ──────────────────────────────────────────────────────────────

def read_header(f):
    def ri32(): return struct.unpack("<i", f.read(4))[0]
    def rf64(): return struct.unpack("<d", f.read(8))[0]

    hdr = {}
    hdr["magic"]          = ri32()
    hdr["version"]        = ri32()
    hdr["N"]               = ri32()
    hdr["N_segment"]      = ri32()
    hdr["segment_start"]  = ri32()
    hdr["segment_end"]    = ri32()
    hdr["Delta"]           = rf64()
    hdr["a"]                = rf64()
    hdr["period_record"]  = ri32()
    return hdr


def load_trajectory(filepath: str):
    with open(filepath, "rb") as f:
        hdr = read_header(f)

        N_seg = hdr["N_segment"]

        frame_bytes = 4 + N_seg * 3 * 4
        raw = f.read()

    n_frames, remainder = divmod(len(raw), frame_bytes)
    if remainder != 0:
        print(f"⚠️  {remainder} leftover bytes ignored (incomplete frame).")

    print(f"[load] Header: N={hdr['N']}, N_segment={N_seg}, "
          f"start={hdr['segment_start']}, end={hdr['segment_end']}, "
          f"Delta={hdr['Delta']}, a={hdr['a']}")
    print(f"[load] {n_frames} frames found in '{filepath}'")

    timesteps = np.empty(n_frames, dtype=np.int32)
    chrom     = np.empty((n_frames, N_seg, 3), dtype=np.float32)

    offset = 0
    for k in range(n_frames):
        timesteps[k] = struct.unpack_from("<i", raw, offset)[0]
        offset += 4

        n_f = N_seg * 3
        chrom[k] = np.frombuffer(raw, dtype="<f4", count=n_f, offset=offset).reshape(N_seg, 3)
        offset += n_f * 4

    return hdr, timesteps, chrom


# ──────────────────────────────────────────────────────────────
# Visualization
# ──────────────────────────────────────────────────────────────

def plot_end_to_end(filepath: str, bead_i: int = 300, bead_j: int = 400):
    """
    Plots the end-to-end vector r(t) = R[bead_j](t) - R[bead_i](t):
      1. |r(t)|^2
      2. Components rx(t), ry(t), rz(t)

    bead_i and bead_j are GLOBAL indices in the chain
    (the script subtracts segment_start to access the right index in the array).
    """
    hdr, timesteps, chrom = load_trajectory(filepath)

    start = hdr["segment_start"]
    N_seg = hdr["N_segment"]

    # Convert global indices -> local indices in chrom
    li = bead_i - start
    lj = bead_j - start

    if not (0 <= li < N_seg and 0 <= lj < N_seg):
        raise ValueError(
            f"Beads {bead_i} and {bead_j} must be within the recorded segment "
            f"[{start}, {start + N_seg - 1}]."
        )

    # End-to-end vector
    r  = chrom[:, lj, :] - chrom[:, li, :]   # (n_frames, 3)
    rx, ry, rz = r[:, 0], r[:, 1], r[:, 2]
    r2 = rx**2 + ry**2 + rz**2
    t  = timesteps.astype(float)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    label = f"beads {bead_i}->{bead_j}"

    # --- |r(t)|^2 ---
    axes[0].plot(t, r2, color="steelblue", linewidth=1.2)
    axes[0].set_xlabel("Timestep")
    axes[0].set_ylabel(r"$|\mathbf{r}_{" + str(bead_i) + r"\to" + str(bead_j) + r"}(t)|^2$")
    axes[0].set_title(f"Squared end-to-end vector — {label}")
    axes[0].grid(True, alpha=0.3)

    # --- Components ---
    axes[1].plot(t, rx, label="x", color="tomato",       linewidth=1.0)
    axes[1].plot(t, ry, label="y", color="seagreen",     linewidth=1.0)
    axes[1].plot(t, rz, label="z", color="mediumpurple", linewidth=1.0)
    axes[1].set_xlabel("Timestep")
    axes[1].set_ylabel("Component")
    axes[1].set_title(f"End-to-end vector components — {label}")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out = filepath.replace(".bin", f"_r2_{bead_i}_{bead_j}.png")
    plt.savefig(out, dpi=150)
    plt.show()
    print(f"[plot] Figure saved to '{out}'")

    return hdr, timesteps, chrom, r2


# ──────────────────────────────────────────────────────────────
# Gaussian statistical analysis
# ──────────────────────────────────────────────────────────────

def plot_gaussian_analysis(filepath: str, bead_i: int = None, bead_j: int = None):
    """
    Statistical analysis of the end-to-end vector r = R[j] - R[i]:

      1. Histograms of rx, ry, rz with a Gaussian fit
         -> checks that each component ~ N(0, sigma^2)
      2. Histogram of |r|^2 with a chi-squared fit (3 degrees of freedom)
         -> for a Gaussian chain: |r|^2/sigma^2 ~ chi^2(3)
      3. Q-Q plot of each component
         -> visual deviation from the Gaussian line
      4. Prints statistics: <|r|^2>, sigma^2, skewness, kurtosis

    Parameters
    ----------
    filepath : path to the binary file
    bead_i, bead_j : global indices (default: segment start and end)
    """
    from scipy import stats

    hdr, timesteps, chrom = load_trajectory(filepath)

    start = hdr["segment_start"]
    N_seg = hdr["N_segment"]

    if bead_i is None: bead_i = start
    if bead_j is None: bead_j = start + N_seg - 1

    li = bead_i - start
    lj = bead_j - start

    if not (0 <= li < N_seg and 0 <= lj < N_seg):
        raise ValueError(
            f"Beads {bead_i} and {bead_j} must be within [{start}, {start + N_seg - 1}]."
        )

    r  = chrom[:, lj, :] - chrom[:, li, :]
    rx, ry, rz = r[:, 0], r[:, 1], r[:, 2]
    r2 = rx**2 + ry**2 + rz**2

    # ── Text statistics ──────────────────────────────────────
    print(f"\n{'─'*50}")
    print(f"  Gaussian analysis: beads {bead_i} -> {bead_j}")
    print(f"{'─'*50}")
    print(f"  <|r|^2>      = {r2.mean():.4f}")
    print(f"  sigma^2  (via rx) = {rx.var():.4f}  (expected: <|r|^2>/3 = {r2.mean()/3:.4f})")
    for name, comp in [("rx", rx), ("ry", ry), ("rz", rz)]:
        sk = stats.skew(comp)
        ku = stats.kurtosis(comp)   # excess kurtosis (0 = Gaussian)
        _, p = stats.normaltest(comp)
        print(f"  {name} : skew={sk:+.3f}  kurt_exc={ku:+.3f}  p_normaltest={p:.3e}"
              + ("  ✅" if p > 0.05 else "  ❌"))
    print(f"{'─'*50}\n")

    # ── Figure ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f"Gaussian analysis — beads {bead_i}->{bead_j}", fontsize=13)

    components = [("rx", rx, "tomato"), ("ry", ry, "seagreen"), ("rz", rz, "mediumpurple")]

    # --- Row 1: component histograms with Gaussian fit ---
    for col, (name, comp, color) in enumerate(components):
        ax = fig.add_subplot(3, 3, col + 1)
        mu, sigma = comp.mean(), comp.std()
        ax.hist(comp, bins=60, density=True, color=color, alpha=0.6, label="data")
        xs = np.linspace(comp.min(), comp.max(), 300)
        ax.plot(xs, stats.norm.pdf(xs, mu, sigma), "k--", linewidth=1.5, label=f"N({mu:.2f}, {sigma:.2f}^2)")
        ax.set_title(f"Histogram {name}")
        ax.set_xlabel(name)
        ax.set_ylabel("Density")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # --- Row 2: Q-Q plots ---
    for col, (name, comp, color) in enumerate(components):
        ax = fig.add_subplot(3, 3, col + 4)
        (osm, osr), (slope, intercept, _) = stats.probplot(comp, dist="norm")
        ax.scatter(osm, osr, s=1, color=color, alpha=0.4)
        xs = np.array([osm[0], osm[-1]])
        ax.plot(xs, slope * xs + intercept, "k--", linewidth=1.5)
        ax.set_title(f"Q-Q plot {name}")
        ax.set_xlabel("Theoretical quantiles")
        ax.set_ylabel("Observed quantiles")
        ax.grid(True, alpha=0.3)

    # --- Row 3: histogram of |r|^2 with chi^2(3)*sigma^2 fit ---
    ax_r2 = fig.add_subplot(3, 3, 7)
    sigma2 = rx.var()   # estimated from rx (assuming isotropy)
    r2_norm = r2 / sigma2  # should follow chi^2(3)
    ax_r2.hist(r2_norm, bins=60, density=True, color="steelblue", alpha=0.6, label=r"$|r|^2/\sigma^2$")
    xs = np.linspace(0, r2_norm.max(), 300)
    ax_r2.plot(xs, stats.chi2.pdf(xs, df=3), "k--", linewidth=1.5, label=r"$\chi^2(3)$")
    ax_r2.set_title(r"Histogram $|r|^2/\sigma^2$")
    ax_r2.set_xlabel(r"$|r|^2 / \sigma^2$")
    ax_r2.set_ylabel("Density")
    ax_r2.legend(fontsize=8)
    ax_r2.grid(True, alpha=0.3)

    # --- Time evolution of |r|^2 ---
    ax_t = fig.add_subplot(3, 3, (8, 9))
    t = timesteps.astype(float)
    ax_t.plot(t, r2, color="steelblue", linewidth=0.6, alpha=0.7)
    ax_t.axhline(r2.mean(), color="k", linestyle="--", linewidth=1.2, label=f"<|r|^2> = {r2.mean():.2f}")
    ax_t.set_xlabel("Timestep")
    ax_t.set_ylabel(r"$|r|^2$")
    ax_t.set_title(r"$|r(t)|^2$ over time")
    ax_t.legend(fontsize=8)
    ax_t.grid(True, alpha=0.3)

    plt.tight_layout()
    out = filepath.replace(".bin", f"_gauss_{bead_i}_{bead_j}.png")
    plt.savefig(out, dpi=150)
    plt.show()
    print(f"[plot] Analysis saved to '{out}'")

    return r, r2


# ──────────────────────────────────────────────────────────────
# R^2 of the whole chain + autocorrelation -> relaxation time tau_r
# ──────────────────────────────────────────────────────────────

def _autocorr_fft(x: np.ndarray) -> np.ndarray:
    """
    Non-normalized autocorrelation of a 1D signal via FFT (Wiener-Khinchin
    theorem), with correction for the number of overlapping points at each
    lag (same as msd_1d_fft on the C side).
    """
    n = len(x)
    nfft = 1
    while nfft < 2 * n:
        nfft <<= 1
    Xf = np.fft.rfft(x, nfft)
    ac = np.fft.irfft(Xf * np.conjugate(Xf), nfft)[:n]
    counts = n - np.arange(n)
    return ac / counts


def compute_end_to_end_autocorrelation(r: np.ndarray, dt_frame: float):
    """
    r : (n_frames, 3) end-to-end vector over time.

    Returns (lag_time, C, tau_r) where:
      C(tau) = <R(0).R(tau)> / <R(0).R(0)>   (time average, centered components)
      tau_r is estimated by linear regression of ln(C(tau)) vs tau, on the
      part of the curve where C in [e^-2, 1] (avoids the noisy tail at large tau).
    """
    n = r.shape[0]
    r_centered = r - r.mean(axis=0, keepdims=True)

    C = np.zeros(n)
    for d in range(3):
        C += _autocorr_fft(r_centered[:, d])
    C /= C[0]  # normalization: C(0) = 1

    lag_time = np.arange(n) * dt_frame

    cutoff = np.exp(-2.0)
    mask = (lag_time > 0) & (C > cutoff)
    if mask.sum() >= 2:
        slope, _ = np.polyfit(lag_time[mask], np.log(C[mask]), 1)
        tau_r = -1.0 / slope if slope < 0 else np.nan
    else:
        tau_r = np.nan

    return lag_time, C, tau_r


def plot_chain_r2_and_autocorrelation(filepath: str):
    """
    1. Plots R^2(t) = |R_end - R_start|^2 for the WHOLE segment recorded
       in the binary trajectory (by default [segment_start, segment_end)).
    2. Plots the autocorrelation function C(tau) of the end-to-end vector, with
       an exponential fit exp(-tau/tau_r) to estimate the chain's relaxation
       time tau_r.

    ⚠️  The binary file only contains the segment [segment_start,
    segment_end) (configured in config.c), not necessarily the full N-monomer
    chain. "Whole chain" here refers to the entire recorded segment.
    To compute tau_r on the full N-monomer chain, regenerate the trajectory
    with segment_start=0 and segment_end=N in config.c before rerunning the
    simulation.
    """
    hdr, timesteps, chrom = load_trajectory(filepath)

    start = hdr["segment_start"]
    end   = hdr["segment_end"]

    # End-to-end vector of the whole recorded segment: last bead - first bead
    r  = chrom[:, -1, :] - chrom[:, 0, :]
    r2 = (r ** 2).sum(axis=1)
    t  = timesteps.astype(float)

    dt_frame = hdr["Delta"] * hdr["period_record"]
    lag_time, C, tau_r = compute_end_to_end_autocorrelation(r, dt_frame)

    # Autocorrelation display window: up to ~10 tau_r (or everything if tau_r unknown)
    if not np.isnan(tau_r) and tau_r > 0:
        n_show = min(len(lag_time), int(10 * tau_r / dt_frame) + 2)
    else:
        n_show = len(lag_time)

    fig, axes = plt.subplots(2, 1, figsize=(11, 9))

    # --- Panel 1: R^2(t) of the whole chain (recorded segment) ---
    axes[0].plot(t, r2, color="darkorange", linewidth=0.8)
    axes[0].axhline(r2.mean(), color="k", linestyle="--", linewidth=1.2,
                    label=f"<R^2> = {r2.mean():.4g}")
    axes[0].set_xlabel("Timestep")
    axes[0].set_ylabel(r"$R^2(t) = |\mathbf{R}_{" + str(end-1) + r"} - \mathbf{R}_{" + str(start) + r"}|^2$")
    axes[0].set_title(f"End-to-end R^2 — recorded segment [{start}, {end-1}]")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # --- Panel 2: autocorrelation C(tau) + exponential fit ---
    axes[1].plot(lag_time[:n_show], C[:n_show], color="steelblue", linewidth=1.2, label=r"measured $C(\tau)$")
    if not np.isnan(tau_r):
        axes[1].plot(lag_time[:n_show], np.exp(-lag_time[:n_show] / tau_r), "k--", linewidth=1.5,
                     label=fr"fit $e^{{-\tau/\tau_r}}$,  $\tau_r$ = {tau_r:.4g}")
    axes[1].axhline(np.exp(-1), color="gray", linestyle=":", linewidth=1, label=r"$C=1/e$")
    axes[1].set_xlabel(r"$\tau$ (reduced time)")
    axes[1].set_ylabel(r"$C(\tau) = \langle \mathbf{R}(0)\cdot\mathbf{R}(\tau)\rangle / \langle R^2\rangle$")
    axes[1].set_title("End-to-end vector autocorrelation function")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out = filepath.replace(".bin", "_R2_chain_autocorr.png")
    plt.savefig(out, dpi=150)
    plt.show()
    print(f"[plot] Whole-chain R^2 + autocorrelation saved to '{out}'")

    if not np.isnan(tau_r):
        print(f"[tau_r] Estimated relaxation time: tau_r ~= {tau_r:.6g} "
              f"(reduced time units = steps x Delta)")
    else:
        print("[tau_r] ❌ Unable to estimate tau_r: the autocorrelation doesn't decay "
              "enough / signal too noisy (try increasing T).")

    return hdr, r2, lag_time, C, tau_r


# ──────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 plot_trajectory.py <trajectory.bin> [bead_i] [bead_j]")
        sys.exit(1)

    path = sys.argv[1]

    with open(path, "rb") as f:
        hdr_preview = read_header(f)
    start = hdr_preview["segment_start"]
    end   = hdr_preview["segment_end"] - 1

    bead_i = int(sys.argv[2]) if len(sys.argv) > 2 else start
    bead_j = int(sys.argv[3]) if len(sys.argv) > 3 else end

    print(f"[info] End-to-end vector: bead {bead_i} -> bead {bead_j}")

    # Trajectory + components
    plot_end_to_end(path, bead_i=bead_i, bead_j=bead_j)

    # Gaussian statistical analysis
    plot_gaussian_analysis(path, bead_i=bead_i, bead_j=bead_j)

    # R^2 of the whole chain (recorded segment) + autocorrelation -> tau_r
    plot_chain_r2_and_autocorrelation(path)
