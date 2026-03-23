"""
PSD Generator for Persson Friction Modeling
============================================
Tab-based GUI application:
  Tab 1 — Single PSD:   Compute C(q) from one road profile
  Tab 2 — Ensemble:     Load multiple profiles, PCA-based C(q) generation

Input format (IDADA road profile):
  Line 1: Profile name
  Line 2: Flag
  Line 3: Number of data points
  Line 4: dx (spacing in meters)
  Line 5: Unit conversion factor (data_unit -> meters)
  Line 6+: x_value  h_value (in data units, e.g. micrometers)
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

# numpy >= 2.0 renamed trapz -> trapezoid
_trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

# Color palette for multi-profile plots
_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
]


# ==============================================================================
# Core PSD Computation Engine (unchanged)
# ==============================================================================

class PSDComputer:
    """Computes Power Spectral Density from 1D road surface profiles."""

    def __init__(self):
        self.profile_name = ""
        self.h_raw = None
        self.dx = None
        self.N = None
        self.L = None
        self.unit_factor = 1e-6

    def load_profile(self, filepath, x_unit='um', h_unit='um'):
        """Load road profile with auto-format detection.

        Supported formats:
          1. IDADA format (5-line header + x h data)
          2. Generic CSV: x,h two columns (comma/space/tab separated)
             with or without a single header line

        Handles encodings: UTF-8, CP949 (Korean Windows), EUC-KR, Latin-1.

        Parameters
        ----------
        filepath : str
        x_unit : str   Unit of x values for generic CSV ('m', 'mm', 'um')
        h_unit : str   Unit of h values for generic CSV ('m', 'mm', 'um')
        """
        self.profile_name = os.path.basename(filepath)

        lines = self._read_with_encoding(filepath)

        if self._try_load_idada(lines):
            pass  # IDADA loaded successfully
        else:
            self._load_generic_csv(lines, x_unit, h_unit)

        return {
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

    @staticmethod
    def _read_with_encoding(filepath):
        """Read file trying multiple encodings (UTF-8, CP949, EUC-KR, Latin-1)."""
        for enc in ['utf-8-sig', 'utf-8', 'cp949', 'euc-kr', 'latin-1']:
            try:
                with open(filepath, 'r', encoding=enc) as f:
                    return f.readlines()
            except (UnicodeDecodeError, UnicodeError):
                continue
        raise ValueError("Cannot decode file with any supported encoding")

    def _try_load_idada(self, lines):
        """Try parsing as IDADA format. Returns True on success."""
        try:
            if len(lines) < 6:
                return False
            name = lines[0].strip()
            int(lines[1].strip())          # flag — must be int
            n_pts = int(lines[2].strip())  # n_points — must be int
            dx = float(lines[3].strip())   # dx in meters
            uf = float(lines[4].strip())   # unit factor
            # Sanity: dx should be a small positive number (meters)
            if dx <= 0 or dx > 1:
                return False
            # Line 6+ should be "x h" pairs
            parts = lines[5].strip().split()
            if len(parts) < 2:
                return False
            float(parts[0]); float(parts[1])  # must be numeric

            self.profile_name = name
            self.dx = dx
            self.unit_factor = uf
            h_list = []
            for line in lines[5:]:
                p = line.strip().split()
                if len(p) >= 2:
                    h_list.append(float(p[1]))
            self.h_raw = np.array(h_list) * self.unit_factor
            self.N = len(self.h_raw)
            self.L = self.N * self.dx
            return True
        except (ValueError, IndexError):
            return False

    @staticmethod
    def _parse_unit_from_header(header_text):
        """Extract unit from header like 'X(um)', 'Z(mm)', 'height(m)'."""
        import re
        m = re.search(r'\((um|mm|m)\)', header_text, re.IGNORECASE)
        if m:
            u = m.group(1).lower()
            # 'um' also matches end of 'mm' sometimes; be explicit
            return u
        # Also try µm / micrometer patterns
        if 'um' in header_text.lower() or 'µm' in header_text:
            return 'um'
        if 'mm' in header_text.lower():
            return 'mm'
        return None

    def _load_generic_csv(self, lines, x_unit='um', h_unit='um'):
        """Load generic two-column CSV with auto-detection of delimiter,
        header units, and missing values.

        Handles:
          - Values wrapped in quotes: "123.45", "-678.9 "
          - Headers like 'X(um)  Z(um)' or '"X(µm)","Z(µm)"'
          - Rows where the second column is empty / '""' (skipped)
          - Comma / tab / semicolon delimiters
          - Windows line endings, BOM, trailing whitespace
        """
        unit_map = {'m': 1.0, 'mm': 1e-3, 'um': 1e-6}

        def _clean(s):
            """Strip quotes, whitespace, BOM from a cell value."""
            s = s.strip().strip('"').strip("'").strip()
            # Remove BOM if present
            s = s.lstrip('\ufeff')
            return s

        # --- Detect delimiter ---
        # Sample lines from various positions
        sample_indices = set()
        for pos in [0.3, 0.5, 0.7]:
            sample_indices.add(min(int(len(lines) * pos), len(lines) - 1))

        delim = None
        for candidate in [',', '\t', ';']:
            for si in sample_indices:
                parts = lines[si].strip().split(candidate)
                if len(parts) >= 2:
                    v0, v1 = _clean(parts[0]), _clean(parts[1])
                    if v1 != '':
                        try:
                            float(v0); float(v1)
                            delim = candidate; break
                        except ValueError:
                            continue
            if delim is not None:
                break
        # fallback: any whitespace
        if delim is None:
            delim = None

        def _split(line):
            return line.strip().split(delim) if delim else line.strip().split()

        # --- Check first line for header / units ---
        first_parts = _split(lines[0])
        if len(first_parts) >= 2:
            h0, h1 = _clean(first_parts[0]), _clean(first_parts[1])
            try:
                float(h0)
            except ValueError:
                # First line is a header — try extracting units
                # Also check raw text for µm / um patterns
                raw_header = lines[0]
                xu = self._parse_unit_from_header(h0) or self._parse_unit_from_header(raw_header)
                hu = self._parse_unit_from_header(h1) or self._parse_unit_from_header(raw_header)
                if xu:
                    x_unit = xu
                if hu:
                    h_unit = hu

        x_scale = unit_map.get(x_unit, 1e-6)
        h_scale = unit_map.get(h_unit, 1e-6)

        # --- Parse data, skip rows with missing h ---
        x_list, h_list = [], []
        for line in lines:
            parts = _split(line)
            if len(parts) < 2:
                continue
            v0, v1 = _clean(parts[0]), _clean(parts[1])
            # Skip if second column is empty
            if v1 == '':
                continue
            try:
                x_list.append(float(v0))
                h_list.append(float(v1))
            except ValueError:
                continue  # skip header or non-numeric lines

        if len(x_list) < 10:
            raise ValueError(f"Too few valid data points ({len(x_list)})")

        x_arr = np.array(x_list) * x_scale   # -> meters
        h_arr = np.array(h_list) * h_scale    # -> meters

        # Compute dx from x values (assume approximately uniform spacing)
        dx_vals = np.diff(x_arr)
        self.dx = np.median(np.abs(dx_vals))
        if self.dx <= 0:
            raise ValueError(f"Invalid spacing dx={self.dx}")

        self.h_raw = h_arr
        self.N = len(self.h_raw)
        self.L = self.N * self.dx

    def compute_psd(self, detrend='linear', window='multitaper', use_top_psd=False,
                    conversion_method='standard', hurst=0.8,
                    correction_factor=1.1615, n_bins=88,
                    welch_nperseg=None, welch_overlap=0.5):
        """Full PSD pipeline: detrend -> window -> FFT -> 1D->2D -> bin.

        window : str
            'none', 'hanning', 'hamming', 'blackman' — apply window then
            single FFT over the full profile.
            'welch' — Welch method: split into overlapping segments, apply
            Hanning window to each, average their PSDs.
        welch_nperseg : int or None
            Segment length for Welch. None = N//8 (clamped to [64, N]).
        welch_overlap : float
            Overlap fraction for Welch segments (0–0.9). Default 0.5.
        """
        if self.h_raw is None:
            raise ValueError("No profile loaded")
        h = self.h_raw.copy()
        h = self._detrend(h, method=detrend)
        if use_top_psd:
            h = self._top_profile(h)

        if window == 'welch':
            q_pos, C1D = self._compute_1d_psd_welch(
                h, nperseg=welch_nperseg, overlap=welch_overlap)
        elif window == 'multitaper':
            q_pos, C1D = self._compute_1d_psd_multitaper(h)
        else:
            h_w = self._apply_window(h, window_type=window)
            q_pos, C1D = self._compute_1d_psd(h_w)

        C2D = self._convert_1d_to_2d(q_pos, C1D, method=conversion_method,
                                      H=hurst, corr=correction_factor)
        q_bin, C2D_bin = self._log_bin(q_pos, C2D, n_bins=n_bins)
        return q_bin, C2D_bin, q_pos, C1D, C2D

    def _detrend(self, h, method='linear'):
        if method == 'mean':
            return h - np.mean(h)
        elif method == 'linear':
            return scipy_detrend(h, type='linear')
        elif method == 'quadratic':
            x_idx = np.arange(len(h), dtype=np.float64)
            coeffs = np.polyfit(x_idx, h, 2)
            return h - np.polyval(coeffs, x_idx)
        return h

    def _top_profile(self, h):
        h_centered = h - np.mean(h)
        h_top = h_centered.copy()
        h_top[h_top < 0] = 0.0
        return h_top

    def _apply_window(self, h, window_type='none'):
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
        correction = np.sqrt(N / np.sum(w ** 2))
        return h * w * correction

    def _compute_1d_psd(self, h):
        N = self.N
        dx = self.dx
        H_fft = fft(h)
        freqs = fftfreq(N, d=dx)
        q = 2.0 * np.pi * freqs
        C1D = (dx / (2.0 * np.pi * N)) * np.abs(H_fft) ** 2
        mask = q > 0
        q_pos = q[mask]
        C1D_pos = C1D[mask]
        # Correct for sinc² rolloff from discrete sampling
        C1D_pos = self._sinc2_correction(q_pos, C1D_pos)
        return q_pos, C1D_pos

    @staticmethod
    def _sinc2_correction(q, C1D):
        """Compensate the sinc²(f·dx) attenuation from rectangular sampling.

        The DFT of discretely-sampled data is attenuated by sinc²(f·dx)
        where sinc(x)=sin(πx)/(πx). This becomes significant near the
        Nyquist frequency, causing the PSD to drop off.
        """
        q_max = q.max()  # ≈ π/dx (Nyquist)
        # x = f·dx = q/(2·q_Nyquist), so at Nyquist x=0.5
        x = q / (2.0 * q_max)
        sinc_val = np.where(x < 1e-10, 1.0, np.sin(np.pi * x) / (np.pi * x))
        sinc2 = sinc_val ** 2
        return C1D / sinc2

    def _compute_1d_psd_welch(self, h, nperseg=None, overlap=0.5):
        """Welch method: split h into overlapping Hanning-windowed segments,
        compute PSD of each, and return the average."""
        N = len(h)
        dx = self.dx

        # Determine segment length
        if nperseg is None:
            nperseg = max(64, N // 8)
        nperseg = min(nperseg, N)

        step = max(1, int(nperseg * (1 - overlap)))
        w = np.hanning(nperseg)
        # Window energy normalization
        S1 = np.sum(w ** 2)

        starts = list(range(0, N - nperseg + 1, step))
        if not starts:
            starts = [0]

        C1D_accum = None
        for s in starts:
            seg = h[s:s + nperseg] * w
            H_fft = fft(seg)
            freqs = fftfreq(nperseg, d=dx)
            q = 2.0 * np.pi * freqs
            # Normalize: dx / (2π × S1) to account for window energy
            C1D_seg = (dx / (2.0 * np.pi * S1)) * np.abs(H_fft) ** 2
            if C1D_accum is None:
                C1D_accum = C1D_seg.copy()
            else:
                C1D_accum += C1D_seg
        C1D_accum /= len(starts)

        mask = q > 0
        q_pos = q[mask]
        C1D_pos = C1D_accum[mask]
        C1D_pos = self._sinc2_correction(q_pos, C1D_pos)
        return q_pos, C1D_pos

    def _compute_1d_psd_multitaper(self, h):
        """Multi-window PSD: average FFTs with different window functions.

        Uses rectangular, Hanning, Hamming, and Blackman windows on the
        full-length profile. Averaging across windows reduces spectral
        leakage while preserving full frequency resolution (unlike Welch
        which shortens segments).
        """
        N = len(h)
        dx = self.dx
        norm = dx / (2.0 * np.pi * N)
        freqs = fftfreq(N, d=dx)
        q = 2.0 * np.pi * freqs

        windows = [
            np.ones(N),            # rectangular
            np.hanning(N),
            np.hamming(N),
            np.blackman(N),
        ]

        C1D_sum = np.zeros(N)
        for w in windows:
            energy_corr = np.sqrt(N / np.sum(w ** 2))
            H_fft = fft(h * w * energy_corr)
            C1D_sum += norm * np.abs(H_fft) ** 2
        C1D_avg = C1D_sum / len(windows)

        mask = q > 0
        q_pos = q[mask]
        C1D_pos = C1D_avg[mask]
        C1D_pos = self._sinc2_correction(q_pos, C1D_pos)
        return q_pos, C1D_pos

    def _convert_1d_to_2d(self, q, C1D, method='standard', H=0.8, corr=1.1615):
        if method == 'standard':
            return C1D / (np.pi * q) * corr
        elif method == 'gamma':
            f = gamma_func(1.0 + H) / (np.sqrt(np.pi) * gamma_func(H + 0.5))
            return (C1D / q) * f
        elif method == 'sqrt':
            return (C1D / q) * np.sqrt(1.0 + 3.0 * H)
        return C1D / (np.pi * q) * corr

    def _log_bin(self, q, C, n_bins=88):
        """Persson-style adaptive log-binning with low-q smoothing.

        At low q where the FFT frequency spacing > desired log-bin width,
        each FFT frequency gets its own bin but with local log-space
        smoothing (geometric mean of ±neighbors) to reduce periodogram
        noise.  At higher q where many FFT points fall in each bin,
        log-uniform averaging is used.
        """
        log_q = np.log10(q)
        lq_min, lq_max = log_q.min(), log_q.max()
        target_dlogq = (lq_max - lq_min) / n_bins

        # Sort q to identify FFT spacing
        idx = np.argsort(q)
        q_s = q[idx]
        C_s = C[idx]
        lq_s = np.log10(q_s)
        lC_s = np.log10(np.maximum(C_s, 1e-50))

        # Phase 1: Individual FFT points at low-q where spacing > target bin width
        bin_q = []
        bin_C = []
        i_transition = 0
        for i in range(len(q_s) - 1):
            spacing = lq_s[i + 1] - lq_s[i]
            if spacing > target_dlogq * 0.7:
                if C_s[i] > 0:
                    bin_q.append(q_s[i])
                    bin_C.append(C_s[i])
                i_transition = i + 1
            else:
                break

        # Phase 2: Log-uniform bins for the dense region
        if i_transition < len(q_s):
            lq_start = lq_s[i_transition]
            n_remaining = max(1, n_bins - len(bin_q))
            edges = np.linspace(lq_start, lq_max, n_remaining + 1)
            lq_tail = lq_s[i_transition:]
            C_tail = C_s[i_transition:]
            for j in range(n_remaining):
                mask = (lq_tail >= edges[j]) & (lq_tail < edges[j + 1])
                if j == n_remaining - 1:  # include right edge in last bin
                    mask = (lq_tail >= edges[j]) & (lq_tail <= edges[j + 1])
                if np.any(mask):
                    avg_C = np.mean(C_tail[mask])
                    if avg_C > 0:
                        bin_q.append(10.0 ** ((edges[j] + edges[j + 1]) / 2.0))
                        bin_C.append(avg_C)

        return np.array(bin_q), np.array(bin_C)


# ==============================================================================
# Shared GUI helpers
# ==============================================================================

def _bind_mousewheel(canvas):
    """Bind mousewheel scrolling to a canvas, handling enter/leave focus."""
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    def _on_mousewheel_linux(event):
        if event.num == 4:
            canvas.yview_scroll(-3, "units")
        elif event.num == 5:
            canvas.yview_scroll(3, "units")

    def _bind(event):
        import sys
        if sys.platform == 'linux':
            canvas.bind_all("<Button-4>", _on_mousewheel_linux)
            canvas.bind_all("<Button-5>", _on_mousewheel_linux)
        else:
            canvas.bind_all("<MouseWheel>", _on_mousewheel)

    def _unbind(event):
        import sys
        if sys.platform == 'linux':
            canvas.unbind_all("<Button-4>")
            canvas.unbind_all("<Button-5>")
        else:
            canvas.unbind_all("<MouseWheel>")

    canvas.bind("<Enter>", _bind)
    canvas.bind("<Leave>", _unbind)


def _add_combo_row(parent, label_text, default, values):
    """Add label + combobox row, return StringVar."""
    row = ttk.Frame(parent)
    row.pack(fill=tk.X, padx=5, pady=2)
    ttk.Label(row, text=label_text).pack(side=tk.LEFT)
    var = tk.StringVar(value=default)
    ttk.Combobox(row, textvariable=var, width=14,
                 values=values, state='readonly').pack(side=tk.RIGHT)
    return var


def _build_psd_param_widgets(parent):
    """Build PSD parameter controls, return dict of tk vars."""
    pre_lf = ttk.LabelFrame(parent, text="PSD Parameters")
    pre_lf.pack(fill=tk.X, padx=5, pady=5)

    v = {}
    v['detrend'] = _add_combo_row(pre_lf, "Detrend:", 'linear',
                                  ['none', 'mean', 'linear', 'quadratic'])
    v['window'] = _add_combo_row(pre_lf, "Window:", 'multitaper',
                                 ['multitaper', 'none', 'welch', 'hanning', 'hamming', 'blackman'])

    row = ttk.Frame(pre_lf)
    row.pack(fill=tk.X, padx=5, pady=2)
    ttk.Label(row, text="Top PSD:").pack(side=tk.LEFT)
    v['top_psd'] = tk.BooleanVar(value=False)
    ttk.Checkbutton(row, variable=v['top_psd']).pack(side=tk.RIGHT)

    conv_lf = ttk.LabelFrame(parent, text="1D -> 2D Conversion")
    conv_lf.pack(fill=tk.X, padx=5, pady=5)
    v['conv_method'] = _add_combo_row(conv_lf, "Method:", 'standard',
                                      ['standard', 'gamma', 'sqrt'])

    row = ttk.Frame(conv_lf)
    row.pack(fill=tk.X, padx=5, pady=2)
    ttk.Label(row, text="Correction:").pack(side=tk.LEFT)
    v['corr'] = tk.DoubleVar(value=1.1615)
    ttk.Entry(row, textvariable=v['corr'], width=10).pack(side=tk.RIGHT)

    row = ttk.Frame(conv_lf)
    row.pack(fill=tk.X, padx=5, pady=2)
    ttk.Label(row, text="Hurst H:").pack(side=tk.LEFT)
    v['hurst'] = tk.DoubleVar(value=0.80)
    hlabel = ttk.Label(row, text="0.80", width=5)
    hlabel.pack(side=tk.RIGHT)
    ttk.Scale(row, from_=0.1, to=1.0, variable=v['hurst'],
              orient=tk.HORIZONTAL,
              command=lambda val: hlabel.config(text=f"{float(val):.2f}")
              ).pack(side=tk.RIGHT, fill=tk.X, expand=True, padx=5)

    bin_lf = ttk.LabelFrame(parent, text="Output Binning")
    bin_lf.pack(fill=tk.X, padx=5, pady=5)
    row = ttk.Frame(bin_lf)
    row.pack(fill=tk.X, padx=5, pady=2)
    ttk.Label(row, text="Log bins:").pack(side=tk.LEFT)
    v['nbins'] = tk.IntVar(value=88)
    ttk.Entry(row, textvariable=v['nbins'], width=6).pack(side=tk.RIGHT)

    return v


def _get_psd_params(v):
    """Extract compute_psd kwargs from a param-var dict."""
    return {
        'detrend': v['detrend'].get(),
        'window': v['window'].get(),
        'use_top_psd': v['top_psd'].get(),
        'conversion_method': v['conv_method'].get(),
        'hurst': v['hurst'].get(),
        'correction_factor': v['corr'].get(),
        'n_bins': int(v['nbins'].get()),
    }


# ==============================================================================
# Tab 1  —  Single Profile PSD
# ==============================================================================

class SinglePSDTab:
    """Original single-profile PSD generator (now inside a tab)."""

    DEFAULT_PROFILE = 'IDADA_road_profile.csv'
    DEFAULT_REFERENCE = 'persson_real_program_psd.csv'

    def __init__(self, parent, root):
        self.root = root
        self.computer = PSDComputer()
        self.profile_loaded = False
        self.reference_data = None
        self.current_psd = None

        frame = ttk.Frame(parent)
        frame.pack(fill=tk.BOTH, expand=True)

        # Left: controls
        ctrl = ttk.Frame(frame, width=340)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        ctrl.pack_propagate(False)
        self._build_controls(ctrl)

        # Right: plots
        plot_frame = ttk.Frame(frame)
        plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._build_plots(plot_frame)

        self._auto_load()

    # --- Controls ---
    def _build_controls(self, parent):
        canvas = tk.Canvas(parent, highlightthickness=0)
        sb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        sf = ttk.Frame(canvas)
        sf.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=sf, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # Resize inner frame width to match canvas
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(
            canvas.find_withtag("all")[0], width=e.width) if canvas.find_withtag("all") else None)
        # Mousewheel scrolling
        _bind_mousewheel(canvas)

        # Profile info
        lf = ttk.LabelFrame(sf, text="Profile Info")
        lf.pack(fill=tk.X, padx=5, pady=5)
        self.info_text = tk.Text(lf, height=9, width=38, state=tk.DISABLED,
                                 font=('Consolas', 9), bg='#f5f5f5')
        self.info_text.pack(fill=tk.X, padx=3, pady=3)

        row = ttk.Frame(sf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Button(row, text="Load Profile", command=self._load_profile_dlg).pack(side=tk.LEFT, padx=2)
        ttk.Button(row, text="Load Reference", command=self._load_ref_dlg).pack(side=tk.LEFT, padx=2)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        self.pv = _build_psd_param_widgets(sf)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        act = ttk.Frame(sf)
        act.pack(fill=tk.X, padx=5, pady=5)
        ttk.Button(act, text="COMPUTE PSD", command=self._compute).pack(fill=tk.X, pady=3)
        ttk.Button(act, text="Export CSV (linear)", command=self._export_csv).pack(fill=tk.X, pady=2)
        ttk.Button(act, text="Export Persson (log10)", command=self._export_log).pack(fill=tk.X, pady=2)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        lf2 = ttk.LabelFrame(sf, text="Comparison Statistics")
        lf2.pack(fill=tk.X, padx=5, pady=5)
        self.stat_text = tk.Text(lf2, height=6, width=38, state=tk.DISABLED,
                                 font=('Consolas', 8), bg='#f5f5f5')
        self.stat_text.pack(fill=tk.X, padx=3, pady=3)

        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(sf, textvariable=self.status_var, wraplength=320,
                  font=('Consolas', 8), foreground='blue').pack(fill=tk.X, padx=5, pady=5)

    # --- Plots ---
    def _build_plots(self, parent):
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax_prof = self.fig.add_subplot(2, 1, 1)
        self.ax_psd = self.fig.add_subplot(2, 1, 2)
        for ax, xl, yl, t in [
            (self.ax_prof, 'Position (mm)', 'Height (um)', 'Road Surface Profile'),
            (self.ax_psd, 'log10(q) [1/m]', 'log10(C(q)) [m^4]', 'PSD C(q)'),
        ]:
            ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(t); ax.grid(True, alpha=0.3)
        self.fig.tight_layout(pad=3.0)
        tb = ttk.Frame(parent); tb.pack(side=tk.TOP, fill=tk.X)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        NavigationToolbar2Tk(self.canvas, tb)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.draw()

    # --- Data ---
    def _auto_load(self):
        d = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(d, self.DEFAULT_PROFILE)
        if os.path.exists(p):
            self._load_profile(p)
        r = os.path.join(d, self.DEFAULT_REFERENCE)
        if os.path.exists(r):
            self._load_ref(r)

    def _load_profile_dlg(self):
        p = filedialog.askopenfilename(title="Select Profile",
                                       filetypes=[("CSV/Text", "*.csv *.txt *.dat"), ("All", "*.*")])
        if p: self._load_profile(p)

    def _load_profile(self, path):
        try:
            info = self.computer.load_profile(path)
            self.profile_loaded = True
            self.info_text.config(state=tk.NORMAL)
            self.info_text.delete('1.0', tk.END)
            self.info_text.insert(tk.END,
                f"Name:     {info['name']}\nPoints:   {info['n_points']:,}\n"
                f"dx:       {info['dx_um']:.4f} um\nLength:   {info['L_mm']:.2f} mm\n"
                f"h_rms:    {info['h_rms_um']:.2f} um\nh_min:    {info['h_min_um']:.1f} um\n"
                f"h_max:    {info['h_max_um']:.1f} um\nq_min:    {info['q_min']:.2f} 1/m\n"
                f"q_max:    {info['q_max']:.2e} 1/m")
            self.info_text.config(state=tk.DISABLED)
            self._plot_profile()
            self.status_var.set(f"Loaded: {info['name']}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _load_ref_dlg(self):
        p = filedialog.askopenfilename(title="Select Reference PSD",
                                       filetypes=[("CSV/Text", "*.csv *.txt"), ("All", "*.*")])
        if p: self._load_ref(p)

    def _load_ref(self, path):
        try:
            d = np.loadtxt(path, delimiter=',')
            if d.ndim == 2 and d.shape[1] >= 2:
                self.reference_data = (d[:, 0], d[:, 1])
                self._update_psd_plot()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # --- Plot ---
    def _plot_profile(self):
        self.ax_prof.clear()
        h = self.computer.h_raw
        x = np.arange(len(h)) * self.computer.dx * 1e3
        hu = h * 1e6
        step = max(1, len(h) // 5000)
        self.ax_prof.plot(x[::step], hu[::step], 'b-', lw=0.4)
        self.ax_prof.set_xlabel('Position (mm)'); self.ax_prof.set_ylabel('Height (um)')
        self.ax_prof.set_title(f'Profile: {self.computer.profile_name}')
        self.ax_prof.grid(True, alpha=0.3); self.canvas.draw()

    def _update_psd_plot(self):
        self.ax_psd.clear()
        if self.current_psd:
            q, C = self.current_psd
            self.ax_psd.plot(np.log10(q), np.log10(C), 'b-o', lw=1.5, ms=2, label='Computed')
        if self.reference_data:
            self.ax_psd.plot(*self.reference_data, 'r-s', lw=1.5, ms=3, label='Reference', alpha=0.8)
        self.ax_psd.set_xlabel('log10(q) [1/m]'); self.ax_psd.set_ylabel('log10(C(q)) [m^4]')
        self.ax_psd.set_title('PSD C(q)'); self.ax_psd.legend(loc='upper right')
        self.ax_psd.grid(True, alpha=0.3); self.canvas.draw()

    # --- Compute ---
    def _compute(self):
        if not self.profile_loaded:
            messagebox.showwarning("Warning", "No profile loaded!"); return
        self.status_var.set("Computing..."); self.root.update_idletasks()
        try:
            q, C, *_ = self.computer.compute_psd(**_get_psd_params(self.pv))
            self.current_psd = (q, C)
            self._update_psd_plot(); self._update_stats(q, C)
            self.status_var.set(f"Done: {len(q)} bins")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _update_stats(self, q, C):
        self.stat_text.config(state=tk.NORMAL); self.stat_text.delete('1.0', tk.END)
        self.stat_text.insert(tk.END, f"q: {q[0]:.1f} - {q[-1]:.2e}\nC: {C[0]:.2e} - {C[-1]:.2e}\n")
        if self.reference_data:
            rq, rC = self.reference_data
            fn = interp1d(np.log10(q), np.log10(C), bounds_error=False, fill_value=np.nan)
            c = fn(rq); v = ~np.isnan(c)
            if np.any(v):
                d = c[v] - rC[v]
                self.stat_text.insert(tk.END,
                    f"\n--- vs Reference ---\nRMSE: {np.sqrt(np.mean(d**2)):.4f}\n"
                    f"MaxErr: {np.max(np.abs(d)):.4f}\nBias: {np.mean(d):+.4f}")
        self.stat_text.config(state=tk.DISABLED)

    # --- Export ---
    def _export_csv(self):
        if not self.current_psd: messagebox.showwarning("Warning", "Compute first!"); return
        p = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if p:
            q, C = self.current_psd
            with open(p, 'w', newline='') as f:
                w = csv.writer(f); w.writerow(['q (1/m)', 'C_2D (m^4)'])
                for qi, ci in zip(q, C): w.writerow([f'{qi:.6E}', f'{ci:.6E}'])

    def _export_log(self):
        if not self.current_psd: messagebox.showwarning("Warning", "Compute first!"); return
        p = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if p:
            q, C = self.current_psd
            with open(p, 'w', newline='') as f:
                w = csv.writer(f)
                for qi, ci in zip(q, C): w.writerow([f'{np.log10(qi):.2E}', f'{np.log10(ci):.2E}'])


# ==============================================================================
# Tab 2  —  Ensemble Analysis
# ==============================================================================

class EnsembleTab:
    """Multi-profile PSD ensemble with PCA-based sample generation."""

    def __init__(self, parent, root):
        self.root = root

        # Data storage
        self.file_paths = []          # list of loaded file paths
        self.file_names = []          # display names
        self.computers = []           # PSDComputer per file
        self.psd_results = []         # list of (q, C2D) per file
        self.q_grid = None            # common q grid (100 pts)
        self.C_matrix = None          # (n_profiles, N_q)

        # PCA results
        self.Y_mean = None
        self.eigenvalues = None
        self.eigenvectors = None
        self.K = 0

        # Samples
        self.samples = None

        frame = ttk.Frame(parent)
        frame.pack(fill=tk.BOTH, expand=True)

        # Left: controls
        ctrl = ttk.Frame(frame, width=360)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        ctrl.pack_propagate(False)
        self._build_controls(ctrl)

        # Right: sub-tabs with plots
        self.plot_nb = ttk.Notebook(frame)
        self.plot_nb.pack(fill=tk.BOTH, expand=True)
        self._build_plot_tabs()

    # ------------------------------------------------------------------
    # Controls
    # ------------------------------------------------------------------

    def _build_controls(self, parent):
        canvas = tk.Canvas(parent, highlightthickness=0)
        sb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        sf = ttk.Frame(canvas)
        sf.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=sf, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # Resize inner frame width to match canvas
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(
            canvas.find_withtag("all")[0], width=e.width) if canvas.find_withtag("all") else None)
        # Mousewheel scrolling
        _bind_mousewheel(canvas)

        # --- File List ---
        lf = ttk.LabelFrame(sf, text="Profile Files")
        lf.pack(fill=tk.X, padx=5, pady=5)

        self.file_listbox = tk.Listbox(lf, height=10, font=('Consolas', 8),
                                        selectmode=tk.EXTENDED)
        self.file_listbox.pack(fill=tk.X, padx=3, pady=3)

        btn_row = ttk.Frame(lf)
        btn_row.pack(fill=tk.X, padx=3, pady=2)
        ttk.Button(btn_row, text="Add Files", command=self._add_files).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row, text="Remove", command=self._remove_files).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row, text="Clear All", command=self._clear_files).pack(side=tk.LEFT, padx=2)

        self.file_count_var = tk.StringVar(value="0 files loaded")
        ttk.Label(lf, textvariable=self.file_count_var,
                  font=('Consolas', 8)).pack(padx=5, pady=2)

        # --- Data Units (for generic CSV files) ---
        unit_lf = ttk.LabelFrame(sf, text="Data Units (auto-detect from header)")
        unit_lf.pack(fill=tk.X, padx=5, pady=5)
        row = ttk.Frame(unit_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="X unit:").pack(side=tk.LEFT)
        self.x_unit_var = tk.StringVar(value='um')
        ttk.Combobox(row, textvariable=self.x_unit_var, width=6,
                     values=['um', 'mm', 'm'], state='readonly').pack(side=tk.RIGHT)
        row = ttk.Frame(unit_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="H unit:").pack(side=tk.LEFT)
        self.h_unit_var = tk.StringVar(value='um')
        ttk.Combobox(row, textvariable=self.h_unit_var, width=6,
                     values=['um', 'mm', 'm'], state='readonly').pack(side=tk.RIGHT)
        ttk.Label(unit_lf, text="(IDADA format auto-detected, units ignored)",
                  font=('Consolas', 7), foreground='gray').pack(padx=5, pady=1)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- PSD Params ---
        self.pv = _build_psd_param_widgets(sf)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- PCA / Generation ---
        pca_lf = ttk.LabelFrame(sf, text="PCA & Generation")
        pca_lf.pack(fill=tk.X, padx=5, pady=5)

        row = ttk.Frame(pca_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="Var threshold:").pack(side=tk.LEFT)
        self.var_thresh_var = tk.DoubleVar(value=0.90)
        ttk.Entry(row, textvariable=self.var_thresh_var, width=6).pack(side=tk.RIGHT)

        row = ttk.Frame(pca_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="N samples:").pack(side=tk.LEFT)
        self.nsamp_var = tk.IntVar(value=1000)
        ttk.Entry(row, textvariable=self.nsamp_var, width=8).pack(side=tk.RIGHT)

        row = ttk.Frame(pca_lf)
        row.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(row, text="Random seed:").pack(side=tk.LEFT)
        self.seed_var = tk.IntVar(value=42)
        ttk.Entry(row, textvariable=self.seed_var, width=8).pack(side=tk.RIGHT)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- Action Buttons ---
        act = ttk.Frame(sf)
        act.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(act, text="1. Compute All PSDs",
                   command=self._compute_all).pack(fill=tk.X, pady=2)
        ttk.Button(act, text="2. Fit PCA",
                   command=self._fit_pca).pack(fill=tk.X, pady=2)
        ttk.Button(act, text="3. Generate Samples",
                   command=self._generate).pack(fill=tk.X, pady=2)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=4)
        ttk.Button(act, text="Export Samples (Persson)",
                   command=lambda: self._export_samples('persson')).pack(fill=tk.X, pady=2)
        ttk.Button(act, text="Export Samples (CSV)",
                   command=lambda: self._export_samples('csv')).pack(fill=tk.X, pady=2)

        ttk.Separator(sf, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5, pady=8)

        # --- Info ---
        info_lf = ttk.LabelFrame(sf, text="Info / Status")
        info_lf.pack(fill=tk.X, padx=5, pady=5)
        self.info_text = tk.Text(info_lf, height=12, width=42, state=tk.DISABLED,
                                 font=('Consolas', 8), bg='#f5f5f5')
        self.info_text.pack(fill=tk.X, padx=3, pady=3)

    # ------------------------------------------------------------------
    # Plot sub-tabs
    # ------------------------------------------------------------------

    def _build_plot_tabs(self):
        """Create 4 sub-tabs: Profiles, PSDs, Ensemble, PCA."""
        tabs = {}
        for name in ['Profiles', 'PSDs', 'Ensemble', 'PCA']:
            f = ttk.Frame(self.plot_nb)
            self.plot_nb.add(f, text=name)
            tabs[name] = f

        # --- Profiles tab ---
        self.fig_prof = Figure(figsize=(10, 7), dpi=100)
        self.ax_prof = self.fig_prof.add_subplot(111)
        self.ax_prof.set_title('Road Surface Profiles')
        self.ax_prof.grid(True, alpha=0.3)
        self.fig_prof.tight_layout()
        tb = ttk.Frame(tabs['Profiles']); tb.pack(side=tk.TOP, fill=tk.X)
        self.cv_prof = FigureCanvasTkAgg(self.fig_prof, master=tabs['Profiles'])
        NavigationToolbar2Tk(self.cv_prof, tb)
        self.cv_prof.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # --- PSDs tab ---
        self.fig_psd = Figure(figsize=(10, 7), dpi=100)
        self.ax_psd = self.fig_psd.add_subplot(111)
        self.ax_psd.set_title('Power Spectral Densities')
        self.ax_psd.grid(True, alpha=0.3)
        self.fig_psd.tight_layout()
        tb2 = ttk.Frame(tabs['PSDs']); tb2.pack(side=tk.TOP, fill=tk.X)
        self.cv_psd = FigureCanvasTkAgg(self.fig_psd, master=tabs['PSDs'])
        NavigationToolbar2Tk(self.cv_psd, tb2)
        self.cv_psd.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # --- Ensemble tab ---
        self.fig_ens = Figure(figsize=(10, 7), dpi=100)
        self.ax_ens = self.fig_ens.add_subplot(111)
        self.ax_ens.set_title('PSD Ensemble')
        self.ax_ens.grid(True, alpha=0.3)
        self.fig_ens.tight_layout()
        tb3 = ttk.Frame(tabs['Ensemble']); tb3.pack(side=tk.TOP, fill=tk.X)
        self.cv_ens = FigureCanvasTkAgg(self.fig_ens, master=tabs['Ensemble'])
        NavigationToolbar2Tk(self.cv_ens, tb3)
        self.cv_ens.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # --- PCA tab ---
        self.fig_pca = Figure(figsize=(10, 7), dpi=100)
        self.ax_pca_vec = self.fig_pca.add_subplot(1, 2, 1)
        self.ax_pca_var = self.fig_pca.add_subplot(1, 2, 2)
        self.fig_pca.tight_layout(pad=3.0)
        tb4 = ttk.Frame(tabs['PCA']); tb4.pack(side=tk.TOP, fill=tk.X)
        self.cv_pca = FigureCanvasTkAgg(self.fig_pca, master=tabs['PCA'])
        NavigationToolbar2Tk(self.cv_pca, tb4)
        self.cv_pca.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ------------------------------------------------------------------
    # File management
    # ------------------------------------------------------------------

    def _add_files(self):
        paths = filedialog.askopenfilenames(
            title="Select Profile Files",
            filetypes=[("CSV/Text/DAT", "*.csv *.txt *.dat"), ("All", "*.*")])
        for p in paths:
            if p not in self.file_paths:
                self.file_paths.append(p)
                name = os.path.basename(p)
                self.file_names.append(name)
                self.file_listbox.insert(tk.END, name)
        self.file_count_var.set(f"{len(self.file_paths)} files loaded")

    def _remove_files(self):
        sel = list(self.file_listbox.curselection())
        for i in reversed(sel):
            self.file_listbox.delete(i)
            del self.file_paths[i]
            del self.file_names[i]
        self.file_count_var.set(f"{len(self.file_paths)} files loaded")

    def _clear_files(self):
        self.file_listbox.delete(0, tk.END)
        self.file_paths.clear()
        self.file_names.clear()
        self.computers.clear()
        self.psd_results.clear()
        self.q_grid = None; self.C_matrix = None
        self.samples = None; self.Y_mean = None
        self.file_count_var.set("0 files loaded")

    # ------------------------------------------------------------------
    # Info helper
    # ------------------------------------------------------------------

    def _set_info(self, text):
        self.info_text.config(state=tk.NORMAL)
        self.info_text.delete('1.0', tk.END)
        self.info_text.insert(tk.END, text)
        self.info_text.config(state=tk.DISABLED)

    def _append_info(self, text):
        self.info_text.config(state=tk.NORMAL)
        self.info_text.insert(tk.END, text)
        self.info_text.config(state=tk.DISABLED)

    # ------------------------------------------------------------------
    # Step 1: Compute all PSDs
    # ------------------------------------------------------------------

    def _compute_all(self):
        if len(self.file_paths) < 2:
            messagebox.showwarning("Warning", "Add at least 2 profile files!")
            return

        self._set_info("Computing PSDs...\n")
        self.root.update_idletasks()

        params = _get_psd_params(self.pv)
        self.computers = []
        self.psd_results = []

        for i, fpath in enumerate(self.file_paths):
            try:
                comp = PSDComputer()
                info = comp.load_profile(fpath,
                                         x_unit=self.x_unit_var.get(),
                                         h_unit=self.h_unit_var.get())
                q, C, *_ = comp.compute_psd(**params)
                self.computers.append(comp)
                self.psd_results.append((q, C))
                self._append_info(f"  [{i+1}] {self.file_names[i]}: {len(q)} bins\n")
                self.root.update_idletasks()
            except Exception as e:
                self._append_info(f"  [{i+1}] FAILED: {e}\n")
                self.computers.append(None)
                self.psd_results.append(None)

        # Build common q grid
        valid = [r for r in self.psd_results if r is not None]
        n_fail = sum(1 for r in self.psd_results if r is None)
        if len(valid) < 2:
            # Collect error details for user
            err_lines = [f"Valid: {len(valid)}, Failed: {n_fail}\n\n"]
            info_content = self.info_text.get('1.0', tk.END)
            for line in info_content.split('\n'):
                if 'FAILED' in line:
                    err_lines.append(line.strip() + '\n')
            messagebox.showerror("Error",
                                 "Need at least 2 valid PSDs\n\n" + "".join(err_lines[:10]))
            return

        q_min = max(r[0][0] for r in valid)
        q_max = min(r[0][-1] for r in valid)
        self.q_grid = np.logspace(np.log10(q_min), np.log10(q_max), 100)

        n_valid = len(valid)
        self.C_matrix = np.zeros((n_valid, 100))
        idx = 0
        for r in self.psd_results:
            if r is None:
                continue
            q_i, C_i = r
            log_C_interp = np.interp(np.log10(self.q_grid),
                                     np.log10(q_i), np.log10(C_i))
            self.C_matrix[idx] = 10.0 ** log_C_interp
            idx += 1

        self._append_info(f"\nCommon q grid: [{self.q_grid[0]:.1f}, {self.q_grid[-1]:.2e}]\n"
                          f"Valid profiles: {n_valid}\n")

        self._plot_profiles()
        self._plot_psds()

    # ------------------------------------------------------------------
    # Step 2: Fit PCA
    # ------------------------------------------------------------------

    def _fit_pca(self):
        if self.C_matrix is None:
            messagebox.showwarning("Warning", "Compute All PSDs first!")
            return

        n = self.C_matrix.shape[0]
        Y = np.log(self.C_matrix)
        self.Y_mean = np.mean(Y, axis=0)
        dY = Y - self.Y_mean

        U, S, Vt = np.linalg.svd(dY, full_matrices=False)
        self.eigenvalues = S ** 2 / (n - 1)
        self.eigenvectors = Vt

        total = np.sum(self.eigenvalues)
        cum = np.cumsum(self.eigenvalues) / total
        thresh = self.var_thresh_var.get()
        self.K = min(int(np.searchsorted(cum, thresh) + 1), len(self.eigenvalues))

        lines = [f"PCA (threshold={thresh*100:.0f}%):\n"]
        for k in range(min(self.K + 2, len(self.eigenvalues))):
            pct = self.eigenvalues[k] / total * 100
            marker = " <--" if k < self.K else ""
            lines.append(f"  PC{k+1}: {pct:5.1f}% (cum {cum[k]*100:5.1f}%){marker}\n")
        lines.append(f"\nSelected K = {self.K}\n")
        self._set_info("".join(lines))

        self._plot_pca()

    # ------------------------------------------------------------------
    # Step 3: Generate samples
    # ------------------------------------------------------------------

    def _generate(self):
        if self.eigenvalues is None:
            messagebox.showwarning("Warning", "Fit PCA first!")
            return

        self._set_info("Generating samples...\n")
        self.root.update_idletasks()

        n_samples = self.nsamp_var.get()
        seed = self.seed_var.get()

        rms_orig = np.array([
            np.sqrt(_trapz(self.C_matrix[i], self.q_grid))
            for i in range(self.C_matrix.shape[0])
        ])
        rms_mean = np.mean(rms_orig)
        rms_std = np.std(rms_orig, ddof=0)
        rms_lo = rms_mean - 3.0 * rms_std
        rms_hi = rms_mean + 3.0 * rms_std

        rng = np.random.default_rng(seed)
        std_k = np.sqrt(self.eigenvalues[:self.K])
        ev_k = self.eigenvectors[:self.K]

        samples = []
        n_rej = 0
        max_att = n_samples * 50

        for _ in range(max_att):
            if len(samples) >= n_samples:
                break
            z = rng.normal(0.0, std_k)
            Y_new = self.Y_mean + z @ ev_k
            C_new = np.exp(Y_new)

            rms_new = np.sqrt(_trapz(C_new, self.q_grid))
            if rms_new < rms_lo or rms_new > rms_hi:
                n_rej += 1; continue
            if np.mean(np.diff(C_new)) > 0:
                n_rej += 1; continue
            samples.append(C_new)

        self.samples = np.array(samples)
        self._set_info(
            f"Generated: {len(self.samples)} / {n_samples}\n"
            f"Rejected:  {n_rej}\n"
            f"RMS range: [{rms_lo:.4e}, {rms_hi:.4e}]\n"
            f"K = {self.K} principal components\n"
        )
        self._plot_ensemble()

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _plot_profiles(self):
        """Plot all loaded road profiles overlaid."""
        self.ax_prof.clear()
        for i, comp in enumerate(self.computers):
            if comp is None:
                continue
            h = comp.h_raw
            x = np.arange(len(h)) * comp.dx * 1e3
            hu = h * 1e6
            step = max(1, len(h) // 3000)
            c = _COLORS[i % len(_COLORS)]
            self.ax_prof.plot(x[::step], hu[::step], color=c, lw=0.4,
                              label=self.file_names[i], alpha=0.8)
        self.ax_prof.set_xlabel('Position (mm)')
        self.ax_prof.set_ylabel('Height (um)')
        self.ax_prof.set_title(f'Road Surface Profiles ({len(self.computers)} files)')
        if len(self.computers) <= 15:
            self.ax_prof.legend(fontsize=7, loc='upper right', ncol=2)
        self.ax_prof.grid(True, alpha=0.3)
        self.cv_prof.draw()
        self.plot_nb.select(0)  # switch to Profiles tab

    def _plot_psds(self):
        """Plot all computed C(q) curves overlaid."""
        self.ax_psd.clear()
        for i, r in enumerate(self.psd_results):
            if r is None:
                continue
            q, C = r
            c = _COLORS[i % len(_COLORS)]
            self.ax_psd.plot(np.log10(q), np.log10(C), color=c, lw=1.2,
                             label=self.file_names[i], alpha=0.8)
        # Plot mean on common grid
        if self.C_matrix is not None:
            C_mean = np.mean(self.C_matrix, axis=0)
            self.ax_psd.plot(np.log10(self.q_grid), np.log10(C_mean),
                             'k--', lw=2.5, label='Mean', zorder=10)
        self.ax_psd.set_xlabel('log10(q) [1/m]')
        self.ax_psd.set_ylabel('log10(C(q)) [m^4]')
        self.ax_psd.set_title(f'PSD — {sum(1 for r in self.psd_results if r)} curves')
        if len(self.psd_results) <= 15:
            self.ax_psd.legend(fontsize=7, loc='upper right', ncol=2)
        self.ax_psd.grid(True, alpha=0.3)
        self.cv_psd.draw()

    def _plot_ensemble(self):
        """Plot generated samples + originals + mean."""
        self.ax_ens.clear()
        lq = np.log10(self.q_grid)

        # Generated samples (blue translucent)
        n_show = len(self.samples)
        # Scale alpha so dense regions don't saturate
        alpha = max(0.03, min(0.15, 15.0 / n_show))
        for i in range(n_show):
            self.ax_ens.plot(lq, np.log10(self.samples[i]),
                             color='steelblue', alpha=alpha, lw=0.5)

        # Originals (grey)
        for i in range(self.C_matrix.shape[0]):
            lbl = 'Originals' if i == 0 else None
            self.ax_ens.plot(lq, np.log10(self.C_matrix[i]),
                             color='grey', alpha=0.7, lw=1.5, label=lbl)

        # Mean
        C_mean = np.exp(self.Y_mean)
        self.ax_ens.plot(lq, np.log10(C_mean), 'r-', lw=2.5, label='Mean')

        self.ax_ens.plot([], [], color='steelblue', alpha=0.5, lw=1,
                         label=f'Generated ({n_show} shown)')

        self.ax_ens.set_xlabel('log10(q) [1/m]')
        self.ax_ens.set_ylabel('log10(C(q)) [m^4]')
        self.ax_ens.set_title(f'Ensemble: {len(self.samples)} samples')
        self.ax_ens.legend(loc='upper right', fontsize=9)
        self.ax_ens.grid(True, alpha=0.3)
        self.cv_ens.draw()
        self.plot_nb.select(2)  # switch to Ensemble tab

    def _plot_pca(self):
        """Plot PC eigenvectors and variance explained."""
        self.ax_pca_vec.clear()
        self.ax_pca_var.clear()
        lq = np.log10(self.q_grid)

        total = np.sum(self.eigenvalues)
        n_show = min(3, len(self.eigenvalues))
        colors = ['tab:blue', 'tab:orange', 'tab:green']

        for k in range(n_show):
            pct = self.eigenvalues[k] / total * 100
            self.ax_pca_vec.plot(lq, self.eigenvectors[k], color=colors[k],
                                 lw=1.5, label=f'PC{k+1} ({pct:.1f}%)')
        self.ax_pca_vec.axhline(0, color='grey', lw=0.5, ls='--')
        self.ax_pca_vec.set_xlabel('log10(q)')
        self.ax_pca_vec.set_ylabel('Eigenvector')
        self.ax_pca_vec.set_title('Principal Components')
        self.ax_pca_vec.legend(fontsize=9)
        self.ax_pca_vec.grid(True, alpha=0.3)

        # Variance bar chart
        n_bar = min(10, len(self.eigenvalues))
        pcts = self.eigenvalues[:n_bar] / total * 100
        cum = np.cumsum(self.eigenvalues[:n_bar]) / total * 100
        x = np.arange(n_bar)
        bars = self.ax_pca_var.bar(x, pcts, color='steelblue', alpha=0.7, label='Individual')
        self.ax_pca_var.plot(x, cum, 'ro-', lw=1.5, label='Cumulative')
        self.ax_pca_var.axhline(self.var_thresh_var.get() * 100, color='red',
                                ls='--', lw=1, alpha=0.5, label=f'Threshold')
        self.ax_pca_var.set_xlabel('PC index')
        self.ax_pca_var.set_ylabel('Variance explained (%)')
        self.ax_pca_var.set_title('Explained Variance')
        self.ax_pca_var.set_xticks(x, [f'PC{i+1}' for i in x], fontsize=7)
        self.ax_pca_var.legend(fontsize=8)
        self.ax_pca_var.grid(True, alpha=0.3, axis='y')

        self.fig_pca.tight_layout()
        self.cv_pca.draw()
        self.plot_nb.select(3)  # switch to PCA tab

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export_samples(self, fmt):
        if self.samples is None:
            messagebox.showwarning("Warning", "Generate samples first!")
            return

        save_dir = filedialog.askdirectory(title="Select output directory")
        if not save_dir:
            return

        for idx in range(len(self.samples)):
            fname = os.path.join(save_dir, f"sample_{idx+1:04d}.csv")
            with open(fname, 'w', newline='') as f:
                w = csv.writer(f)
                if fmt == 'csv':
                    w.writerow(['q (1/m)', 'C_2D (m^4)'])
                    for qi, ci in zip(self.q_grid, self.samples[idx]):
                        w.writerow([f'{qi:.6E}', f'{ci:.6E}'])
                else:
                    for qi, ci in zip(self.q_grid, self.samples[idx]):
                        w.writerow([f'{np.log10(qi):.2E}', f'{np.log10(ci):.2E}'])

        self._append_info(f"\nExported {len(self.samples)} samples to\n{save_dir}/\n")


# ==============================================================================
# Main Application  —  Notebook with 2 tabs
# ==============================================================================

class PSDGeneratorApp:
    """Main application window with tab-based layout."""

    def __init__(self, root):
        self.root = root
        self.root.title("PSD Generator - Persson Friction Model")
        self.root.geometry("1500x950")

        nb = ttk.Notebook(root)
        nb.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Tab 1: Single PSD
        tab1 = ttk.Frame(nb)
        nb.add(tab1, text="  Single PSD  ")
        SinglePSDTab(tab1, root)

        # Tab 2: Ensemble
        tab2 = ttk.Frame(nb)
        nb.add(tab2, text="  Ensemble Analysis  ")
        EnsembleTab(tab2, root)


def main():
    root = tk.Tk()
    app = PSDGeneratorApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
