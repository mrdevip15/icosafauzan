#!/usr/bin/env python3
"""
Finite-size scaling analysis for layer-separated quasi-3D simulation.
Analyzes critical behavior and scaling properties.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from scipy.optimize import curve_fit

def load_all_data():
    """Load all simulation data."""
    data_dir = "simulation_runs"
    layer_data = {}
    
    for size_dir in glob.glob(f"{data_dir}/nx*_ny*"):
        size_name = os.path.basename(size_dir)
        parts = size_name.split('_')
        nx = int(parts[0][2:])
        ny = int(parts[1][2:])
        
        if nx == ny:
            layer_data[nx] = {}
            for layer in [1, 2, 3]:
                layer_file = f"{size_dir}/layer_{layer}.txt"
                if os.path.exists(layer_file):
                    data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                     names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                    layer_data[nx][layer] = data
    
    return layer_data

def find_critical_temps():
    """Find critical temperatures using peak heat capacity method."""
    layer_data = load_all_data()
    lattice_sizes = sorted(layer_data.keys())
    
    tc_data = {}
    
    for size in lattice_sizes:
        tc_data[size] = {}
        
        for layer in [1, 2, 3]:
            if layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            
            # Focus on reasonable temperature range
            mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 2.0)
            temp_range = data['Temperature'][mask]
            cv_range = data['Cv'][mask]
            
            # Find peak (excluding anomalous high values)
            reasonable_cv = cv_range < 100  # Filter out computational artifacts
            if reasonable_cv.any():
                filtered_temp = temp_range[reasonable_cv]
                filtered_cv = cv_range[reasonable_cv]
                
                if len(filtered_cv) > 0:
                    peak_idx = filtered_cv.idxmax()
                    tc = filtered_temp.loc[peak_idx]
                    tc_data[size][layer] = tc
    
    return tc_data

def power_law(x, a, b):
    """Power law function for fitting."""
    return a * x**b

def plot_finite_size_scaling():
    """Plot finite-size scaling analysis."""
    tc_data = find_critical_temps()
    layer_data = load_all_data()
    lattice_sizes = sorted(layer_data.keys())
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    
    print("=== Finite-Size Scaling Analysis ===")
    
    # Top row: Critical temperature scaling
    for layer in [1, 2, 3]:
        ax = axes[0, layer-1]
        
        sizes = []
        tcs = []
        
        for size in lattice_sizes:
            if size in tc_data and layer in tc_data[size]:
                sizes.append(size)
                tcs.append(tc_data[size][layer])
        
        if len(sizes) > 2:
            sizes = np.array(sizes)
            tcs = np.array(tcs)
            
            # Plot Tc vs L
            ax.plot(sizes, tcs, 'o-', color=layer_colors[layer-1], 
                   linewidth=2, markersize=8, label='Data')
            
            # Fit to scaling form: Tc(L) = Tc_inf + a/L^(1/nu)
            # For 2D Ising: nu = 1, so 1/nu = 1
            try:
                # Fit inverse relationship
                inv_sizes = 1.0 / sizes
                popt, _ = curve_fit(lambda x, a, b: a + b*x, inv_sizes, tcs)
                tc_inf, coeff = popt
                
                # Plot fit
                x_fit = np.linspace(inv_sizes.max(), inv_sizes.min(), 100)
                y_fit = tc_inf + coeff * x_fit
                size_fit = 1.0 / x_fit
                ax.plot(size_fit, y_fit, '--', color=layer_colors[layer-1], 
                       alpha=0.7, label=f'Fit: $T_c(∞) = {tc_inf:.3f}$')
                
                print(f"Layer {layer}: Tc(∞) = {tc_inf:.3f} ± finite-size correction")
                
            except:
                print(f"Layer {layer}: Fitting failed")
        
        ax.set_xlabel('Lattice Size $L$')
        ax.set_ylabel('Critical Temperature $T_c$')
        ax.set_title(f'Layer {layer}')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    # Bottom row: Peak heat capacity scaling
    for layer in [1, 2, 3]:
        ax = axes[1, layer-1]
        
        sizes = []
        max_cvs = []
        
        for size in lattice_sizes:
            if size not in layer_data or layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 2.0)
            cv_range = data['Cv'][mask]
            
            # Filter reasonable values
            reasonable_mask = cv_range < 100
            if reasonable_mask.any():
                max_cv = cv_range[reasonable_mask].max()
                sizes.append(size)
                max_cvs.append(max_cv)
        
        if len(sizes) > 2:
            sizes = np.array(sizes)
            max_cvs = np.array(max_cvs)
            
            # Plot Cv_max vs L
            ax.loglog(sizes, max_cvs, 'o-', color=layer_colors[layer-1], 
                     linewidth=2, markersize=8, label='Data')
            
            # For 2D Ising: Cv_max ~ L^(alpha/nu) with alpha = 0, nu = 1
            # But for Potts models, there can be logarithmic corrections
            try:
                log_sizes = np.log(sizes)
                log_cvs = np.log(max_cvs)
                
                # Fit power law
                popt, _ = curve_fit(lambda x, a, b: a + b*x, log_sizes, log_cvs)
                intercept, slope = popt
                
                # Plot fit
                x_fit = np.logspace(np.log10(sizes.min()), np.log10(sizes.max()), 100)
                y_fit = np.exp(intercept) * x_fit**slope
                ax.loglog(x_fit, y_fit, '--', color=layer_colors[layer-1], 
                         alpha=0.7, label=f'$C_v ∝ L^{{{slope:.2f}}}$')
                
                print(f"Layer {layer}: Cv_max ∝ L^{slope:.3f}")
                
            except:
                print(f"Layer {layer}: Heat capacity scaling fit failed")
        
        ax.set_xlabel('Lattice Size $L$')
        ax.set_ylabel('Peak Heat Capacity $C_{v,max}$')
        ax.set_title(f'Layer {layer}')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig('finite_size_scaling_analysis.png', dpi=300, bbox_inches='tight')
    print("\nSaved: finite_size_scaling_analysis.png")
    plt.show()

def analyze_universality():
    """Analyze universal scaling behavior."""
    tc_data = find_critical_temps()
    layer_data = load_all_data()
    
    print("\n=== Universality Analysis ===")
    
    # Compare critical temperatures between layers
    lattice_sizes = sorted(layer_data.keys())
    
    print("Critical temperature comparison:")
    print("Size\tLayer 1\tLayer 2\tLayer 3\tΔT(1-2)\tΔT(1-3)\tΔT(2-3)")
    print("-" * 70)
    
    for size in lattice_sizes:
        if size in tc_data:
            tcs = []
            for layer in [1, 2, 3]:
                if layer in tc_data[size]:
                    tcs.append(tc_data[size][layer])
                else:
                    tcs.append(float('nan'))
            
            if len([t for t in tcs if not np.isnan(t)]) >= 2:
                t1, t2, t3 = tcs
                dt12 = t1 - t2 if not (np.isnan(t1) or np.isnan(t2)) else float('nan')
                dt13 = t1 - t3 if not (np.isnan(t1) or np.isnan(t3)) else float('nan')
                dt23 = t2 - t3 if not (np.isnan(t2) or np.isnan(t3)) else float('nan')
                
                print(f"{size}\t{t1:.3f}\t{t2:.3f}\t{t3:.3f}\t{dt12:.3f}\t{dt13:.3f}\t{dt23:.3f}")

def plot_binder_cumulant():
    """Plot Binder cumulant for each layer."""
    layer_data = load_all_data()
    lattice_sizes = sorted(layer_data.keys())
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    
    for layer in [1, 2, 3]:
        ax = axes[layer-1]
        
        for i, size in enumerate(lattice_sizes):
            if size not in layer_data or layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            
            # Calculate Binder cumulant: U = 1 - <M^4>/(3<M^2>^2)
            mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 1.8)
            mask &= (data['M2'] > 1e-10) & (data['M4'] > 1e-10)
            
            if mask.sum() > 10:
                temp = data['Temperature'][mask]
                m2 = data['M2'][mask]
                m4 = data['M4'][mask]
                
                binder = 1 - m4 / (3 * m2**2)
                
                # Filter reasonable values
                reasonable = (binder > 0) & (binder < 1) & np.isfinite(binder)
                if reasonable.sum() > 5:
                    ax.plot(temp[reasonable], binder[reasonable], 
                           linewidth=2, alpha=0.8, label=f'L = {size}')
        
        # Add universal value line for 2D Ising
        ax.axhline(y=0.61, color='black', linestyle='--', alpha=0.5, 
                  label='2D Ising Universal')
        
        ax.set_xlabel('Temperature $T$')
        ax.set_ylabel('Binder Cumulant $U$')
        ax.set_title(f'Layer {layer}')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig('binder_cumulant_analysis.png', dpi=300, bbox_inches='tight')
    print("Saved: binder_cumulant_analysis.png")
    plt.show()

def main():
    """Main analysis function."""
    print("=== Advanced Finite-Size Scaling Analysis ===")
    print("Quasi-3D 8-state cube spin model analysis")
    print("=" * 50)
    
    plot_finite_size_scaling()
    analyze_universality()
    plot_binder_cumulant()
    
    print("\n=== Analysis Complete ===")
    print("Generated plots:")
    print("- finite_size_scaling_analysis.png")
    print("- binder_cumulant_analysis.png")

if __name__ == "__main__":
    main()