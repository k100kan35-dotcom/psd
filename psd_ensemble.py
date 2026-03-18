"""
PSD Ensemble Generator with PCA-based Sampling
================================================
Constructs a PSD ensemble from multiple road surface profile files,
performs PCA decomposition, and generates physically plausible C(q)
samples for Monte-Carlo friction simulations.

Depends on PSDComputer from psd_generator.py for individual PSD computation.
This module is GUI-free (no tkinter dependency at runtime).
"""

import sys
import types
import os
import csv
import numpy as np

# numpy >= 2.0 renamed trapz → trapezoid
_trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _import_psd_computer():
    """Import PSDComputer from psd_generator.py, handling tkinter dependency.

    psd_generator.py imports tkinter and sets matplotlib to TkAgg at module
    level (for its GUI).  In headless environments we must:
      1. Stub out tkinter and its submodules
      2. Pre-set matplotlib backend to 'Agg' so the TkAgg call is a no-op
    Only PSDComputer (pure computation, no GUI) is used here.
    """
    try:
        from psd_generator import PSDComputer
        return PSDComputer
    except (ImportError, ModuleNotFoundError):
        pass

    # Build comprehensive tkinter + matplotlib backend stubs.
    # psd_generator.py imports:
    #   - tkinter (tk, ttk, filedialog, messagebox)
    #   - matplotlib.use('TkAgg')
    #   - matplotlib.backends.backend_tkagg (FigureCanvasTkAgg, NavigationToolbar2Tk)
    # All of these transitively touch many tkinter.* submodules.

    class _Stub:
        """Universal stub: callable, subscriptable, attribute-accessible."""
        def __init__(self, *a, **kw): pass
        def __call__(self, *a, **kw): return _Stub()
        def __getattr__(self, name): return _Stub()
        def __bool__(self): return False

    class _StubModule(types.ModuleType):
        """Module stub that returns _Stub for any attribute access."""
        def __getattr__(self, name):
            return _Stub()

    # Register stubs for tkinter and ALL submodules that matplotlib's
    # _backend_tk might try to import
    _tk_submodules = [
        'tkinter', 'tkinter.ttk', 'tkinter.filedialog',
        'tkinter.messagebox', 'tkinter.font', 'tkinter.constants',
        'tkinter.simpledialog', 'tkinter.colorchooser',
    ]
    for mod_name in _tk_submodules:
        sys.modules[mod_name] = _StubModule(mod_name)

    # Also stub the matplotlib TkAgg backend so its top-level import
    # in psd_generator doesn't chain-fail
    for mod_name in [
        'matplotlib.backends.backend_tkagg',
        'matplotlib.backends._backend_tk',
    ]:
        sys.modules[mod_name] = _StubModule(mod_name)

    # Patch matplotlib.use to silently ignore TkAgg
    _orig_use = matplotlib.use
    matplotlib.use = lambda backend, **kw: (
        None if backend.lower() == 'tkagg' else _orig_use(backend, **kw)
    )

    try:
        from psd_generator import PSDComputer
    finally:
        matplotlib.use = _orig_use

    return PSDComputer


PSDComputer = _import_psd_computer()


class PSDEnsemble:
    """PCA-based PSD ensemble generator.

    Workflow:
        1. load_profiles()   — compute C(q) for each profile, interpolate to common grid
        2. fit_pca()         — SVD in log-space, select principal components
        3. generate_samples() — draw N(0, lambda_k) coefficients, filter for physics
        4. plot_ensemble()   — visualize original + generated PSDs
        5. export_samples()  — save generated samples to disk
    """

    def __init__(self, psd_params=None):
        """Initialise ensemble generator.

        Parameters
        ----------
        psd_params : dict or None
            Keyword arguments forwarded to ``PSDComputer.compute_psd()``.
            If None, the calibrated defaults from psd_generator are used.
        """
        self.psd_params = psd_params if psd_params is not None else {}

        # Populated by load_profiles
        self.q_grid = None          # shape (N_q,)
        self.C_matrix = None        # shape (n_profiles, N_q)  — linear C(q)
        self.n_profiles = 0

        # Populated by fit_pca
        self.Y_mean = None          # shape (N_q,)
        self.eigenvalues = None     # shape (n_components,)
        self.eigenvectors = None    # shape (n_components, N_q)
        self.K = 0

        # Populated by generate_samples
        self.samples = None         # shape (n_samples, N_q)

    # ------------------------------------------------------------------
    # 1. load_profiles
    # ------------------------------------------------------------------

    def load_profiles(self, filepath_list):
        """Load multiple IDADA-format profiles and build a common-q PSD matrix.

        For each file the full PSD pipeline (detrend → FFT → 1D→2D → binning)
        is executed via ``PSDComputer``.  All resulting C(q) curves are then
        interpolated onto a shared log-spaced q grid (intersection of all
        individual q ranges, 100 points).

        Parameters
        ----------
        filepath_list : list of str
            Paths to IDADA-format profile files.

        Raises
        ------
        ValueError
            If fewer than 2 profiles are provided.
        """
        if len(filepath_list) < 2:
            raise ValueError("At least 2 profiles are required for ensemble analysis")

        q_list = []
        C_list = []

        for fpath in filepath_list:
            comp = PSDComputer()
            comp.load_profile(fpath)
            q_bin, C2D_bin, *_ = comp.compute_psd(**self.psd_params)
            q_list.append(q_bin)
            C_list.append(C2D_bin)
            print(f"  Loaded: {os.path.basename(fpath)}  "
                  f"({len(q_bin)} bins, q=[{q_bin[0]:.1f}, {q_bin[-1]:.2e}])")

        # Common q range: intersection of all individual ranges
        q_min_common = max(q[0] for q in q_list)
        q_max_common = min(q[-1] for q in q_list)
        if q_min_common >= q_max_common:
            raise ValueError("No overlapping q range among the provided profiles")

        self.q_grid = np.logspace(np.log10(q_min_common),
                                  np.log10(q_max_common), 100)

        # Interpolate each C(q) onto common grid (log-log linear interpolation)
        self.C_matrix = np.zeros((len(filepath_list), len(self.q_grid)))
        for i, (q_i, C_i) in enumerate(zip(q_list, C_list)):
            log_C_interp = np.interp(np.log10(self.q_grid),
                                     np.log10(q_i), np.log10(C_i))
            self.C_matrix[i] = 10.0 ** log_C_interp

        self.n_profiles = len(filepath_list)
        print(f"\nEnsemble: {self.n_profiles} profiles on {len(self.q_grid)}-point "
              f"q grid [{self.q_grid[0]:.1f}, {self.q_grid[-1]:.2e}] 1/m")

    # ------------------------------------------------------------------
    # 2. fit_pca
    # ------------------------------------------------------------------

    def fit_pca(self, var_threshold=0.90):
        """Perform PCA (via SVD) on the log-transformed PSD matrix.

        Steps:
            1. Y = ln(C_matrix)
            2. Centre: dY = Y - mean(Y)
            3. SVD:  dY = U S Vt
            4. eigenvalues = S² / (n-1)
            5. Select K components capturing ≥ var_threshold of total variance

        Parameters
        ----------
        var_threshold : float
            Cumulative explained-variance fraction for component selection
            (default 0.90 = 90 %).
        """
        if self.C_matrix is None:
            raise RuntimeError("Call load_profiles() first")

        Y = np.log(self.C_matrix)                   # (n, N_q)
        self.Y_mean = np.mean(Y, axis=0)             # (N_q,)
        dY = Y - self.Y_mean                          # (n, N_q)

        U, S, Vt = np.linalg.svd(dY, full_matrices=False)

        self.eigenvalues = S ** 2 / (self.n_profiles - 1)
        self.eigenvectors = Vt                         # (n_components, N_q)

        # Cumulative variance ratio
        total_var = np.sum(self.eigenvalues)
        cum_var = np.cumsum(self.eigenvalues) / total_var

        self.K = int(np.searchsorted(cum_var, var_threshold) + 1)
        self.K = min(self.K, len(self.eigenvalues))

        print(f"\nPCA results (threshold={var_threshold*100:.0f}%):")
        print(f"  Total components: {len(self.eigenvalues)}")
        print(f"  Selected K = {self.K}")
        for k in range(min(self.K + 2, len(self.eigenvalues))):
            marker = " <--" if k < self.K else ""
            print(f"  PC{k+1}: {self.eigenvalues[k]:.4e}  "
                  f"({self.eigenvalues[k]/total_var*100:5.1f}% / "
                  f"cum {cum_var[k]*100:5.1f}%){marker}")

    # ------------------------------------------------------------------
    # 3. generate_samples
    # ------------------------------------------------------------------

    def generate_samples(self, n_samples=1000, random_seed=None):
        """Generate physically plausible C(q) samples via PCA reconstruction.

        For each sample:
            z_k ~ N(0, eigenvalue_k),   k = 0 … K-1
            Y_new = Y_mean + Σ z_k · eigenvector_k
            C_new = exp(Y_new)

        Rejection filters (re-draw on failure):
            * RMS roughness outside mean ± 3σ of the originals
            * C(q) not monotonically decreasing on average

        Parameters
        ----------
        n_samples : int
            Number of valid samples to generate (default 1000).
        random_seed : int or None
            Seed for reproducibility.

        Returns
        -------
        np.ndarray, shape (n_samples, N_q)
            Generated C(q) values (linear scale).
        """
        if self.eigenvalues is None:
            raise RuntimeError("Call fit_pca() first")

        # Reference RMS statistics from original profiles
        rms_orig = np.array([
            np.sqrt(_trapz(self.C_matrix[i], self.q_grid))
            for i in range(self.n_profiles)
        ])
        rms_mean = np.mean(rms_orig)
        rms_std = np.std(rms_orig, ddof=0)
        rms_lo = rms_mean - 3.0 * rms_std
        rms_hi = rms_mean + 3.0 * rms_std

        rng = np.random.default_rng(random_seed)
        std_k = np.sqrt(self.eigenvalues[:self.K])  # shape (K,)
        ev_k = self.eigenvectors[:self.K]            # shape (K, N_q)

        samples = []
        n_rejected = 0
        max_attempts = n_samples * 50  # safety cap
        attempts = 0

        while len(samples) < n_samples and attempts < max_attempts:
            attempts += 1

            # Draw random coefficients
            z = rng.normal(0.0, std_k)               # shape (K,)

            # Reconstruct in log-space
            Y_new = self.Y_mean + z @ ev_k            # shape (N_q,)
            C_new = np.exp(Y_new)                      # shape (N_q,)

            # --- Filter 1: RMS roughness ---
            rms_new = np.sqrt(_trapz(C_new, self.q_grid))
            if rms_new < rms_lo or rms_new > rms_hi:
                n_rejected += 1
                continue

            # --- Filter 2: Monotonic decrease (on average) ---
            dC = np.diff(C_new)
            if np.mean(dC) > 0:
                n_rejected += 1
                continue

            samples.append(C_new)

        if len(samples) < n_samples:
            print(f"WARNING: only {len(samples)}/{n_samples} valid samples generated "
                  f"after {max_attempts} attempts")

        self.samples = np.array(samples)
        print(f"\nSample generation: {len(self.samples)} accepted, "
              f"{n_rejected} rejected ({attempts} total attempts)")
        print(f"  RMS filter range: [{rms_lo:.4e}, {rms_hi:.4e}]")
        return self.samples

    # ------------------------------------------------------------------
    # 4. plot_ensemble
    # ------------------------------------------------------------------

    def plot_ensemble(self, n_show=100, save_path=None):
        """Visualise original profiles, generated samples, and PC directions.

        Left subplot : log10(q) vs log10(C)
            — originals (grey), generated (blue, translucent), mean (red)
        Right subplot: PC1–PC3 eigenvectors with explained variance

        Parameters
        ----------
        n_show : int
            Number of generated samples to display (default 100).
        save_path : str or None
            If given, save figure as PNG.
        """
        if self.samples is None:
            raise RuntimeError("Call generate_samples() first")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        log_q = np.log10(self.q_grid)

        # --- Left: PSD ensemble ---
        # Generated samples (blue, translucent)
        n_plot = min(n_show, len(self.samples))
        for i in range(n_plot):
            ax1.plot(log_q, np.log10(self.samples[i]),
                     color='steelblue', alpha=0.08, linewidth=0.5)

        # Originals (grey, thicker)
        for i in range(self.n_profiles):
            label = 'Original profiles' if i == 0 else None
            ax1.plot(log_q, np.log10(self.C_matrix[i]),
                     color='grey', alpha=0.7, linewidth=1.5, label=label)

        # Mean (red)
        C_mean = np.exp(self.Y_mean)
        ax1.plot(log_q, np.log10(C_mean),
                 'r-', linewidth=2.5, label='Mean (log-space)')

        # Dummy handle for generated
        ax1.plot([], [], color='steelblue', alpha=0.5, linewidth=1,
                 label=f'Generated ({n_plot} shown)')

        ax1.set_xlabel('log10(q)  [q in 1/m]')
        ax1.set_ylabel('log10(C(q))  [C in m^4]')
        ax1.set_title('PSD Ensemble')
        ax1.legend(loc='upper right', fontsize=9)
        ax1.grid(True, alpha=0.3)

        # --- Right: Principal Components ---
        total_var = np.sum(self.eigenvalues)
        n_pc_show = min(3, len(self.eigenvalues))
        colors = ['tab:blue', 'tab:orange', 'tab:green']

        for k in range(n_pc_show):
            pct = self.eigenvalues[k] / total_var * 100
            ax2.plot(log_q, self.eigenvectors[k],
                     color=colors[k], linewidth=1.5,
                     label=f'PC{k+1} ({pct:.1f}%)')

        ax2.axhline(0, color='grey', linewidth=0.5, linestyle='--')
        ax2.set_xlabel('log10(q)  [q in 1/m]')
        ax2.set_ylabel('Eigenvector amplitude')
        ax2.set_title('Principal Components (log-space)')
        ax2.legend(loc='best', fontsize=9)
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=200, bbox_inches='tight')
            print(f"Figure saved: {save_path}")

        plt.close(fig)
        return fig

    # ------------------------------------------------------------------
    # 5. export_samples
    # ------------------------------------------------------------------

    def export_samples(self, save_dir, format='persson'):
        """Export generated C(q) samples to individual CSV files.

        Parameters
        ----------
        save_dir : str
            Directory to write files into (created if needed).
        format : str
            'persson' — two columns: log10(q), log10(C)  (no header)
            'csv'     — two columns: q (1/m), C_2D (m^4) (with header)
        """
        if self.samples is None:
            raise RuntimeError("Call generate_samples() first")

        os.makedirs(save_dir, exist_ok=True)

        for idx in range(len(self.samples)):
            fname = os.path.join(save_dir, f"sample_{idx+1:04d}.csv")

            with open(fname, 'w', newline='') as f:
                writer = csv.writer(f)

                if format == 'csv':
                    writer.writerow(['q (1/m)', 'C_2D (m^4)'])
                    for qi, ci in zip(self.q_grid, self.samples[idx]):
                        writer.writerow([f'{qi:.6E}', f'{ci:.6E}'])
                else:  # persson
                    for qi, ci in zip(self.q_grid, self.samples[idx]):
                        writer.writerow([f'{np.log10(qi):.2E}',
                                         f'{np.log10(ci):.2E}'])

        print(f"Exported {len(self.samples)} samples to {save_dir}/ "
              f"(format={format})")
