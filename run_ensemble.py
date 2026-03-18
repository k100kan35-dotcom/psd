#!/usr/bin/env python3
"""
Run PSD Ensemble Generation
============================
Example script demonstrating the full PSDEnsemble workflow:
  1. Load multiple road surface profiles
  2. Fit PCA in log-space
  3. Generate 1000 physically plausible C(q) samples
  4. Visualise and export results

Usage:
    python run_ensemble.py

Profile files are expected in the current directory or a 'profiles/' subdirectory.
Edit the `files` list below to match your actual profile paths.
"""

import glob
import os
import sys
from psd_ensemble import PSDEnsemble


def main():
    # ----- Locate profile files -----
    # Try common locations; edit as needed for your setup
    search_patterns = [
        "profiles/*.dat",
        "profiles/*.csv",
        "profiles/*.txt",
        "*.dat",
    ]
    files = []
    for pat in search_patterns:
        files = sorted(glob.glob(pat))
        if len(files) >= 2:
            break

    # Fallback: if only a single IDADA profile exists, duplicate it for demo
    if len(files) < 2:
        single = "IDADA_road_profile.csv"
        if os.path.exists(single):
            print(f"NOTE: Only 1 profile found ({single}).")
            print("      Using it multiple times for demonstration.\n")
            files = [single] * 5
        else:
            print("ERROR: No profile files found.")
            print("Place IDADA-format profiles in ./profiles/ or current directory.")
            sys.exit(1)

    print(f"=== PSD Ensemble Generator ===\n")
    print(f"Profiles: {len(files)} files")

    # ----- PSD computation parameters (optional override) -----
    psd_params = {
        'detrend': 'linear',
        'window': 'none',
        'use_top_psd': False,
        'conversion_method': 'standard',
        'correction_factor': 1.1615,
        'n_bins': 88,
    }

    # ----- Run ensemble pipeline -----
    ens = PSDEnsemble(psd_params=psd_params)

    # Step 1: Load profiles and build common-q PSD matrix
    ens.load_profiles(files)

    # Step 2: PCA decomposition
    ens.fit_pca(var_threshold=0.90)

    # Step 3: Generate 1000 C(q) samples
    samples = ens.generate_samples(n_samples=1000, random_seed=42)
    print(f"  Generated sample shape: {samples.shape}")

    # Step 4: Visualise
    ens.plot_ensemble(n_show=100, save_path="ensemble_result.png")

    # Step 5: Export
    ens.export_samples(save_dir="psd_samples/", format='persson')

    print("\n=== Done ===")


if __name__ == '__main__':
    main()
