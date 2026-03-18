"""
PSD Generator for Persson Friction Modeling
============================================
GUI application for computing Power Spectral Density (PSD) from road surface profiles.
Generates C(q) data compatible with Persson's contact mechanics / RubberFriction programs.

Input format (IDADA road profile):
  Line 1: Profile name
  Line 2: Flag
  Line 3: Number of data points
  Line 4: dx (spacing in meters)
  Line 5: Unit conversion factor (data_unit -> meters)
  Line 6+: x_value  h_value (in data units, e.g. micrometers)

Output: log10(q) vs log10(C_2D) matching Persson program format
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import detrend as scipy_detrend
from scipy.special import gamma as gamma_func
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import os
import csv


# ==============================================================================
# Core PSD Computation Engine
# ==============================================================================

class PSDComputer:
    """Computes Power Spectral Density from 1D road surface profiles.

    Algorithm (optimized to match Persson program output):
    1. Load profile, convert to meters
    2. Linear detrend (remove mean + slope)
    3. (Optional) Top profile extraction for asymmetric surfaces
    4. (Optional) Window function
    5. FFT -> two-sided 1D PSD: C_1D(q) = dx/(2*pi*N) * |FFT|^2
    6. 1D->2D conversion: C_2D(q) = C_1D(q) / (pi*q) * correction_factor
    7. Arithmetic mean binning in log-spaced q bins
    """

    def __init__(self):
        self.profile_name = ""
        self.h_raw = None       # raw height array (meters)
        self.dx = None           # spacing (meters)
        self.N = None            # number of points
        self.L = None            # total scan length (meters)
        self.unit_factor = 1e-6  # default: micrometers -> meters

    def load_profile(self, filepath):
        """Load road profile from IDADA-format file."""
        with open(filepath, 'r') as f:
            lines = f.readlines()

        self.profile_name = lines[0].strip()
        n_points_header = int(lines[2].strip())
        self.dx = float(lines[3].strip())         # meters
        self.unit_factor = float(lines[4].strip()) # conversion to meters

        h_list = []
        for line in lines[5:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                h_list.append(float(parts[1]))

        self.h_raw = np.array(h_list) * self.unit_factor  # convert to meters
        self.N = len(self.h_raw)
        self.L = self.N * self.dx

        info = {
            'name': self.profile_name,
            'n_points': self.N,
            'dx_um': self.dx * 1e6,
            'L_mm': self.L * 1e3,
            'h_rms_um': np.std(self.h_raw) * 1e6,
            'h_min_um': np.min(self.h_raw) * 1e6,
            'h_max_um': np.max(self.h_raw) * 1e6,
            'q_min': 2 * np.pi / self.L,
            'q_max': np.pi / self.dx,
        }
        return info

    def compute_psd(self, detrend='linear', window='none', use_top_psd=False,
                    conversion_method='standard', hurst=0.8,
                    correction_factor=1.1615, n_bins=88):
        """
        Full PSD computation pipeline.

        Parameters:
        -----------
        detrend : str
            'none', 'mean', 'linear', 'quadratic'
        window : str
            'none', 'hanning', 'hamming', 'blackman'
        use_top_psd : bool
            If True, zero out negative heights (for asymmetric surfaces)
        conversion_method : str
            'standard': C_2D = C_1D / (pi*q) * correction_factor
            'gamma': C_2D = C_1D / q * Gamma(1+H)/(sqrt(pi)*Gamma(H+0.5))
            'sqrt': C_2D = C_1D / q * sqrt(1+3H)
        hurst : float
            Hurst exponent (used by gamma and sqrt methods)
        correction_factor : float
            Multiplicative correction for standard method (default 1.1615)
        n_bins : int
            Number of logarithmic bins

        Returns: (q_binned, C2D_binned, q_raw, C1D_raw, C2D_raw)
        """
        if self.h_raw is None:
            raise ValueError("No profile loaded")

        h = self.h_raw.copy()

        # Step 1: Detrend
        h = self._detrend(h, method=detrend)

        # Step 2: Top PSD (optional)
        if use_top_psd:
            h = self._top_profile(h)

        # Step 3: Window
        h_w = self._apply_window(h, window_type=window)

        # Step 4: FFT -> two-sided 1D PSD
        q_pos, C1D = self._compute_1d_psd(h_w)

        # Step 5: 1D -> 2D conversion
        C2D = self._convert_1d_to_2d(q_pos, C1D, method=conversion_method,
                                      H=hurst, corr=correction_factor)

        # Step 6: Log-space binning (arithmetic mean)
        q_bin, C2D_bin = self._log_bin(q_pos, C2D, n_bins=n_bins)

        return q_bin, C2D_bin, q_pos, C1D, C2D

    def _detrend(self, h, method='linear'):
        """Remove trend from profile."""
        if method == 'mean':
            return h - np.mean(h)
        elif method == 'linear':
            return scipy_detrend(h, type='linear')
        elif method == 'quadratic':
            x_idx = np.arange(len(h), dtype=np.float64)
            coeffs = np.polyfit(x_idx, h, 2)
            return h - np.polyval(coeffs, x_idx)
        return h  # 'none'

    def _top_profile(self, h):
        """Extract top profile: keep only heights above mean, zero out the rest."""
        h_centered = h - np.mean(h)
        h_top = h_centered.copy()
        h_top[h_top < 0] = 0.0
        return h_top

    def _apply_window(self, h, window_type='none'):
        """Apply window function with energy correction."""
        N = len(h)
        if window_type == 'none':
            return h.copy()

        if window_type == 'hanning':
            w = np.hanning(N)
        elif window_type == 'hamming':
            w = np.hamming(N)
        elif window_type == 'blackman':
            w = np.blackman(N)
        else:
            return h.copy()

        # Energy correction factor
        correction = np.sqrt(N / np.sum(w ** 2))
        return h * w * correction

    def _compute_1d_psd(self, h):
        """
        Compute two-sided 1D PSD using FFT.

        C_1D(q) = dx / (2*pi*N) * |FFT(h)|^2

        Returns only positive q values.
        """
        N = self.N
        dx = self.dx

        H_fft = fft(h)
        freqs = fftfreq(N, d=dx)
        q = 2.0 * np.pi * freqs

        # Two-sided PSD
        C1D = (dx / (2.0 * np.pi * N)) * np.abs(H_fft) ** 2

        mask = q > 0
        return q[mask], C1D[mask]

    def _convert_1d_to_2d(self, q, C1D, method='standard', H=0.8, corr=1.1615):
        """
        Convert 1D PSD to 2D isotropic PSD.

        Methods:
        - 'standard': C_2D = C_1D / (pi*q) * correction_factor
          The base formula C_1D/(pi*q) is the simplest isotropic conversion.
          The correction_factor (default 1.1615) is calibrated to match
          Persson's program output.
        - 'gamma': C_2D = C_1D / q * Gamma(1+H) / (sqrt(pi)*Gamma(H+0.5))
          Exact analytical conversion for self-affine fractal surfaces.
        - 'sqrt': C_2D = C_1D / q * sqrt(1+3H)
          Literature approximation from reference code.
        """
        if method == 'standard':
            C2D = C1D / (np.pi * q) * corr
        elif method == 'gamma':
            factor = gamma_func(1.0 + H) / (np.sqrt(np.pi) * gamma_func(H + 0.5))
            C2D = (C1D / q) * factor
        elif method == 'sqrt':
            C2D = (C1D / q) * np.sqrt(1.0 + 3.0 * H)
        else:
            C2D = C1D / (np.pi * q) * corr

        return C2D

    def _log_bin(self, q, C, n_bins=88):
        """Logarithmically bin PSD using arithmetic mean within each bin."""
        log_q = np.log10(q)
        log_q_min = log_q.min()
        log_q_max = log_q.max()

        edges = np.linspace(log_q_min, log_q_max, n_bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2.0

        C_binned = np.full(n_bins, np.nan)
        for i in range(n_bins):
            mask = (log_q >= edges[i]) & (log_q < edges[i + 1])
            if np.any(mask):
                C_binned[i] = np.mean(C[mask])  # arithmetic mean

        valid = ~np.isnan(C_binned) & (C_binned > 0)
        return 10.0 ** centers[valid], C_binned[valid]


# ==============================================================================
# GUI Application
# ==============================================================================

class PSDGeneratorApp:
    """Tkinter GUI for PSD generation and visualization."""

    DEFAULT_PROFILE = 'IDADA_road_profile.csv'
    DEFAULT_REFERENCE = 'persson_real_program_psd.csv'

    def __init__(self, root):
        self.root = root
        self.root.title("PSD Generator - Persson Friction Model")
        self.root.geometry("1500x950")

        self.computer = PSDComputer()
        self.profile_loaded = False
        self.reference_data = None   # (log10_q, log10_C) arrays
        self.current_psd = None      # (q, C2D) linear arrays
        self.current_psd_raw = None  # (q_raw, C1D, C2D_raw) for diagnostics

        self._build_gui()
        self._auto_load_data()

    # ------------------------------------------------------------------
    # GUI Construction
    # ------------------------------------------------------------------

    def _build_gui(self):
        """Build complete GUI layout."""
        main = ttk.Frame(self.root)
        main.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel: controls (fixed width)
        ctrl = ttk.Frame(main, width=340)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        ctrl.pack_propagate(False)
        self._build_controls(ctrl)

        # Right panel: plots
        plot_frame = ttk.Frame(main)
        plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._build_plots(plot_frame)

    def _build_controls(self, parent):
        """Build the control panel."""
        canvas = tk.Canvas(parent)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scroll_frame = ttk.Frame(canvas)
        scroll_frame.bind("<Configure>",
                          lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        p = scroll_frame  # shorthand

        # --- Profile Info ---
        info_lf = ttk.LabelFrame(p, text="Profile Info")
        info_lf.pack(fill=tk.X, padx=5, pady=5)

        self.info_text = tk.Text(info_lf, height=9, width=38, state=tk.DISABLED,
                                 font=('Consolas', 9), bg='#f5f5f5')
        self.info_text.pack(fill=tk.X, padx=3, pady=3)

        btn_row = ttk.Frame(p)
        btn_row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Button(btn_row, text="Load Profile",
                   command=self._load_profile_dialog).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row, text="Load Reference",
                   command=self._load_reference_dialog).pack(side=tk.LEFT, padx=2)

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- Preprocessing ---
        pre_lf = ttk.LabelFrame(p, text="Preprocessing")
        pre_lf.pack(fill=tk.X, padx=5, pady=5)

        self._add_combo_row(pre_lf, "Detrend:", 'detrend_var', 'linear',
                            ['none', 'mean', 'linear', 'quadratic'])
        self._add_combo_row(pre_lf, "Window:", 'window_var', 'none',
                            ['none', 'hanning', 'hamming', 'blackman'])

        row = ttk.Frame(pre_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="Top PSD:").pack(side=tk.LEFT)
        self.top_psd_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row, variable=self.top_psd_var,
                        text="(for asymmetric surfaces)").pack(side=tk.RIGHT)

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- 1D -> 2D Conversion ---
        conv_lf = ttk.LabelFrame(p, text="1D -> 2D Conversion")
        conv_lf.pack(fill=tk.X, padx=5, pady=5)

        self._add_combo_row(conv_lf, "Method:", 'conv_method_var', 'standard',
                            ['standard', 'gamma', 'sqrt'])

        self.method_desc = ttk.Label(conv_lf, text="", wraplength=300,
                                     font=('Consolas', 8), foreground='gray')
        self.method_desc.pack(fill=tk.X, padx=5, pady=2)
        self.conv_method_var.trace_add('write', self._update_method_desc)

        # Correction factor (for standard method)
        row = ttk.Frame(conv_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="Correction:").pack(side=tk.LEFT)
        self.corr_var = tk.DoubleVar(value=1.1615)
        self.corr_entry = ttk.Entry(row, textvariable=self.corr_var, width=10)
        self.corr_entry.pack(side=tk.RIGHT)

        # Hurst exponent
        row = ttk.Frame(conv_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="Hurst H:").pack(side=tk.LEFT)
        self.hurst_var = tk.DoubleVar(value=0.80)
        self.hurst_label = ttk.Label(row, text="0.80", width=5)
        self.hurst_label.pack(side=tk.RIGHT)
        ttk.Scale(row, from_=0.1, to=1.0, variable=self.hurst_var,
                  orient=tk.HORIZONTAL,
                  command=lambda v: self.hurst_label.config(
                      text=f"{float(v):.2f}")).pack(
            side=tk.RIGHT, fill=tk.X, expand=True, padx=5)

        self._update_method_desc()

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- Binning ---
        bin_lf = ttk.LabelFrame(p, text="Output Binning")
        bin_lf.pack(fill=tk.X, padx=5, pady=5)

        row = ttk.Frame(bin_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="Log bins:").pack(side=tk.LEFT)
        self.nbins_var = tk.IntVar(value=88)
        self.nbins_entry = ttk.Entry(row, textvariable=self.nbins_var, width=6)
        self.nbins_entry.pack(side=tk.RIGHT)

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- Action Buttons ---
        act_frame = ttk.Frame(p)
        act_frame.pack(fill=tk.X, padx=5, pady=5)

        style = ttk.Style()
        style.configure('Compute.TButton', font=('Helvetica', 10, 'bold'))

        ttk.Button(act_frame, text="COMPUTE PSD", style='Compute.TButton',
                   command=self._compute).pack(fill=tk.X, pady=3)
        ttk.Button(act_frame, text="Export CSV (linear)",
                   command=self._export_csv).pack(fill=tk.X, pady=2)
        ttk.Button(act_frame, text="Export Persson Format (log10)",
                   command=self._export_persson_format).pack(fill=tk.X, pady=2)

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- Statistics / Status ---
        stat_lf = ttk.LabelFrame(p, text="Comparison Statistics")
        stat_lf.pack(fill=tk.X, padx=5, pady=5)
        self.stat_text = tk.Text(stat_lf, height=6, width=38, state=tk.DISABLED,
                                 font=('Consolas', 8), bg='#f5f5f5')
        self.stat_text.pack(fill=tk.X, padx=3, pady=3)

        self.status_var = tk.StringVar(value="Ready. Load a profile to begin.")
        ttk.Label(p, textvariable=self.status_var, wraplength=320,
                  font=('Consolas', 8), foreground='blue').pack(
            fill=tk.X, padx=5, pady=5)

    def _add_combo_row(self, parent, label, var_name, default, values):
        """Helper to add a label + combobox row."""
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text=label).pack(side=tk.LEFT)
        var = tk.StringVar(value=default)
        setattr(self, var_name, var)
        ttk.Combobox(row, textvariable=var, width=14,
                     values=values, state='readonly').pack(side=tk.RIGHT)
        return var

    def _build_plots(self, parent):
        """Build matplotlib plot area."""
        self.fig = Figure(figsize=(10, 8), dpi=100)

        # Profile subplot (top)
        self.ax_profile = self.fig.add_subplot(2, 1, 1)
        self.ax_profile.set_xlabel('Position (mm)')
        self.ax_profile.set_ylabel('Height (um)')
        self.ax_profile.set_title('Road Surface Profile')
        self.ax_profile.grid(True, alpha=0.3)

        # PSD subplot (bottom)
        self.ax_psd = self.fig.add_subplot(2, 1, 2)
        self.ax_psd.set_xlabel('log10(q)  [q in 1/m]')
        self.ax_psd.set_ylabel('log10(C(q))  [C in m^4]')
        self.ax_psd.set_title('2D Power Spectral Density C(q)')
        self.ax_psd.grid(True, alpha=0.3)

        self.fig.tight_layout(pad=3.0)

        # Toolbar
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.draw()

    # ------------------------------------------------------------------
    # Data Loading
    # ------------------------------------------------------------------

    def _auto_load_data(self):
        """Auto-load profile and reference if available in script directory."""
        script_dir = os.path.dirname(os.path.abspath(__file__))

        profile_path = os.path.join(script_dir, self.DEFAULT_PROFILE)
        if os.path.exists(profile_path):
            self._load_profile(profile_path)

        ref_path = os.path.join(script_dir, self.DEFAULT_REFERENCE)
        if os.path.exists(ref_path):
            self._load_reference(ref_path)

    def _load_profile_dialog(self):
        path = filedialog.askopenfilename(
            title="Select Road Profile",
            filetypes=[("CSV/Text", "*.csv *.txt"), ("All", "*.*")]
        )
        if path:
            self._load_profile(path)

    def _load_profile(self, filepath):
        try:
            info = self.computer.load_profile(filepath)
            self.profile_loaded = True

            self.info_text.config(state=tk.NORMAL)
            self.info_text.delete('1.0', tk.END)
            self.info_text.insert(tk.END,
                f"Name:     {info['name']}\n"
                f"Points:   {info['n_points']:,}\n"
                f"dx:       {info['dx_um']:.4f} um\n"
                f"Length:   {info['L_mm']:.2f} mm\n"
                f"h_rms:    {info['h_rms_um']:.2f} um\n"
                f"h_min:    {info['h_min_um']:.1f} um\n"
                f"h_max:    {info['h_max_um']:.1f} um\n"
                f"q_min:    {info['q_min']:.2f} 1/m\n"
                f"q_max:    {info['q_max']:.2e} 1/m"
            )
            self.info_text.config(state=tk.DISABLED)

            self._plot_profile()
            self.status_var.set(f"Profile loaded: {info['name']}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load profile:\n{e}")

    def _load_reference_dialog(self):
        path = filedialog.askopenfilename(
            title="Select Reference PSD (log10 format)",
            filetypes=[("CSV/Text", "*.csv *.txt"), ("All", "*.*")]
        )
        if path:
            self._load_reference(path)

    def _load_reference(self, filepath):
        try:
            data = np.loadtxt(filepath, delimiter=',')
            if data.ndim == 2 and data.shape[1] >= 2:
                self.reference_data = (data[:, 0], data[:, 1])
                self._update_psd_plot()
                self.status_var.set("Reference PSD loaded.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load reference:\n{e}")

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _plot_profile(self):
        self.ax_profile.clear()
        h = self.computer.h_raw
        x_mm = np.arange(len(h)) * self.computer.dx * 1e3
        h_um = h * 1e6

        # Downsample for responsive plotting
        if len(h) > 5000:
            step = max(1, len(h) // 5000)
            self.ax_profile.plot(x_mm[::step], h_um[::step], 'b-', linewidth=0.4)
        else:
            self.ax_profile.plot(x_mm, h_um, 'b-', linewidth=0.5)

        self.ax_profile.set_xlabel('Position (mm)')
        self.ax_profile.set_ylabel('Height (um)')
        self.ax_profile.set_title(f'Road Surface Profile: {self.computer.profile_name}')
        self.ax_profile.grid(True, alpha=0.3)
        self.canvas.draw()

    def _update_psd_plot(self):
        """Redraw PSD plot with current data and reference."""
        self.ax_psd.clear()

        if self.current_psd is not None:
            q, C = self.current_psd
            log_q = np.log10(q)
            log_C = np.log10(C)
            self.ax_psd.plot(log_q, log_C, 'b-o', linewidth=1.5, markersize=2,
                            label='Computed C(q)', zorder=3)

        if self.reference_data is not None:
            rq, rC = self.reference_data
            self.ax_psd.plot(rq, rC, 'r-s', linewidth=1.5, markersize=3,
                            label='Persson Reference', alpha=0.8, zorder=2)

        self.ax_psd.set_xlabel('log10(q)  [q in 1/m]')
        self.ax_psd.set_ylabel('log10(C(q))  [C in m^4]')
        self.ax_psd.set_title('2D Power Spectral Density C(q)')
        self.ax_psd.legend(loc='upper right')
        self.ax_psd.grid(True, alpha=0.3)
        self.canvas.draw()

    # ------------------------------------------------------------------
    # Computation
    # ------------------------------------------------------------------

    def _compute(self):
        if not self.profile_loaded:
            messagebox.showwarning("Warning", "No profile loaded!")
            return

        self.status_var.set("Computing PSD...")
        self.root.update_idletasks()

        try:
            params = {
                'detrend': self.detrend_var.get(),
                'window': self.window_var.get(),
                'use_top_psd': self.top_psd_var.get(),
                'conversion_method': self.conv_method_var.get(),
                'hurst': self.hurst_var.get(),
                'correction_factor': self.corr_var.get(),
                'n_bins': int(self.nbins_var.get()),
            }

            q_bin, C2D_bin, q_raw, C1D_raw, C2D_raw = self.computer.compute_psd(**params)

            self.current_psd = (q_bin, C2D_bin)
            self.current_psd_raw = (q_raw, C1D_raw, C2D_raw)

            self._update_psd_plot()
            self._update_statistics(q_bin, C2D_bin)

            self.status_var.set(
                f"Done: {len(q_bin)} bins, "
                f"q=[{q_bin[0]:.1f}, {q_bin[-1]:.2e}] 1/m"
            )
        except Exception as e:
            messagebox.showerror("Error", f"Computation failed:\n{e}")
            import traceback
            traceback.print_exc()

    def _update_statistics(self, q, C):
        """Compute and display comparison statistics."""
        self.stat_text.config(state=tk.NORMAL)
        self.stat_text.delete('1.0', tk.END)

        log_q = np.log10(q)
        log_C = np.log10(C)

        # Basic stats
        self.stat_text.insert(tk.END,
            f"q range: {q[0]:.1f} - {q[-1]:.2e} 1/m\n"
            f"C range: {C[0]:.2e} - {C[-1]:.2e} m^4\n"
        )

        # Comparison with reference
        if self.reference_data is not None:
            ref_lq, ref_lC = self.reference_data
            try:
                fn = interp1d(log_q, log_C, kind='linear',
                              bounds_error=False, fill_value=np.nan)
                c_at_ref = fn(ref_lq)
                v = ~np.isnan(c_at_ref)
                if np.any(v):
                    diff = c_at_ref[v] - ref_lC[v]
                    rmse = np.sqrt(np.mean(diff ** 2))
                    max_err = np.max(np.abs(diff))
                    bias = np.mean(diff)
                    self.stat_text.insert(tk.END,
                        f"\n--- vs Reference ---\n"
                        f"RMSE(log10):   {rmse:.4f}\n"
                        f"MaxErr(log10): {max_err:.4f}\n"
                        f"Bias(log10):   {bias:+.4f}"
                    )
            except Exception:
                pass

        self.stat_text.config(state=tk.DISABLED)

    def _update_method_desc(self, *args):
        m = self.conv_method_var.get()
        descs = {
            'standard': 'C2D = C1D/(pi*q) * corr_factor\n'
                        'Default corr=1.1615 (calibrated to Persson)',
            'gamma': 'C2D = C1D/q * G(1+H)/(sqrt(pi)*G(H+0.5))\n'
                     'Exact for self-affine fractal surfaces',
            'sqrt': 'C2D = C1D/q * sqrt(1+3H)\n'
                    'Literature approximation',
        }
        self.method_desc.config(text=descs.get(m, ''))

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export_csv(self):
        """Export PSD as CSV with headers (linear q and C values)."""
        if self.current_psd is None:
            messagebox.showwarning("Warning", "Compute PSD first!")
            return

        path = filedialog.asksaveasfilename(
            title="Export PSD (Linear CSV)",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")]
        )
        if path:
            q, C = self.current_psd
            with open(path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['q (1/m)', 'C_2D (m^4)'])
                for qi, ci in zip(q, C):
                    writer.writerow([f'{qi:.6E}', f'{ci:.6E}'])
            self.status_var.set(f"Exported: {os.path.basename(path)}")

    def _export_persson_format(self):
        """Export in Persson format: log10(q), log10(C) without headers."""
        if self.current_psd is None:
            messagebox.showwarning("Warning", "Compute PSD first!")
            return

        path = filedialog.asksaveasfilename(
            title="Export PSD (Persson log10 Format)",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")]
        )
        if path:
            q, C = self.current_psd
            with open(path, 'w', newline='') as f:
                writer = csv.writer(f)
                for qi, ci in zip(q, C):
                    writer.writerow([f'{np.log10(qi):.2E}', f'{np.log10(ci):.2E}'])
            self.status_var.set(f"Exported Persson format: {os.path.basename(path)}")


# ==============================================================================
# Entry Point
# ==============================================================================

def main():
    root = tk.Tk()
    app = PSDGeneratorApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
