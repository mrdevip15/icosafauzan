#!/usr/bin/env python3
"""
Focused analysis of specific heat behavior and critical temperature determination.
Detailed examination of heat capacity peaks and critical scaling.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import warnings
warnings.filterwarnings('ignore')

def load_simulation_data():
    """Load all layer simulation data."""
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

def analyze_specific_heat_peaks():
    """Detailed analysis of specific heat peaks for critical temperature determination."""
    layer_data = load_simulation_data()
    lattice_sizes = sorted(layer_data.keys())
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    
    print("=== SPECIFIC HEAT CRITICAL TEMPERATURE ANALYSIS ===")
    print("Analyzing heat capacity peaks for each layer and lattice size")
    print("="*60)
    
    # Create comprehensive specific heat plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Raw specific heat data
    ax1 = axes[0, 0]
    ax1.set_title('Specific Heat vs Temperature - All Layers & Sizes')
    ax1.set_xlabel('Temperature T')
    ax1.set_ylabel('Specific Heat $C_v$')
    
    # Plot 2: Zoomed critical region
    ax2 = axes[0, 1] 
    ax2.set_title('Critical Region (T = 0.8 - 2.0)')
    ax2.set_xlabel('Temperature T')
    ax2.set_ylabel('Specific Heat $C_v$')
    
    # Plot 3: Peak positions vs lattice size
    ax3 = axes[1, 0]
    ax3.set_title('Critical Temperature vs Lattice Size')
    ax3.set_xlabel('Lattice Size L')
    ax3.set_ylabel('Critical Temperature $T_c$')
    
    # Plot 4: Peak heights vs lattice size (log scale)
    ax4 = axes[1, 1]
    ax4.set_title('Peak Height vs Lattice Size')
    ax4.set_xlabel('Lattice Size L')
    ax4.set_ylabel('Peak $C_v$ Value')
    ax4.set_yscale('log')
    
    peak_analysis = {}
    
    for layer in [1, 2, 3]:
        peak_analysis[layer] = {}
        
        print(f"\n--- LAYER {layer} ANALYSIS ---")
        
        sizes_for_plot = []
        tcs_for_plot = []
        peak_cvs_for_plot = []
        
        for size in lattice_sizes:
            if size not in layer_data or layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            
            # Filter temperature range and reasonable Cv values
            temp_mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.5)
            cv_mask = (data['Cv'] > 0) & (data['Cv'] < 10000)  # Remove artifacts
            combined_mask = temp_mask & cv_mask
            
            if combined_mask.sum() < 10:
                continue
            
            temp = data['Temperature'][combined_mask].values
            cv = data['Cv'][combined_mask].values
            
            # Plot raw data
            alpha = 0.6 + 0.1 * (lattice_sizes.index(size))
            linestyle = '-' if size == lattice_sizes[-1] else '--'
            
            ax1.plot(temp, cv, color=layer_colors[layer-1], alpha=alpha, 
                    linewidth=1.5, linestyle=linestyle,
                    label=f'L={size}, Layer {layer}' if size == lattice_sizes[-1] else '')
            
            # Zoomed critical region
            critical_mask = (temp >= 0.8) & (temp <= 2.0)
            if critical_mask.sum() > 5:
                ax2.plot(temp[critical_mask], cv[critical_mask], 
                        color=layer_colors[layer-1], alpha=alpha,
                        linewidth=2, linestyle=linestyle,
                        label=f'L={size}, Layer {layer}' if size == lattice_sizes[-1] else '')
            
            # Find peak
            if len(cv) > 10:
                # Smooth data for peak finding
                if len(cv) > 20:
                    cv_smooth = savgol_filter(cv, window_length=min(11, len(cv)//3*2+1), polyorder=2)
                else:
                    cv_smooth = cv
                
                peak_idx = np.argmax(cv_smooth)
                tc_peak = temp[peak_idx]
                cv_peak = cv_smooth[peak_idx]
                
                peak_analysis[layer][size] = {
                    'tc': tc_peak,
                    'cv_max': cv_peak,
                    'temperature_range': (temp.min(), temp.max()),
                    'cv_range': (cv.min(), cv.max())
                }
                
                sizes_for_plot.append(size)
                tcs_for_plot.append(tc_peak)
                peak_cvs_for_plot.append(cv_peak)
                
                print(f"Size {size}x{size}: Tc = {tc_peak:.4f}, Peak Cv = {cv_peak:.3f}")
        
        # Plot critical temperature vs size
        if len(sizes_for_plot) > 1:
            ax3.plot(sizes_for_plot, tcs_for_plot, 'o-', 
                    color=layer_colors[layer-1], linewidth=2, markersize=8,
                    label=f'Layer {layer}')
            
            # Finite-size extrapolation
            if len(sizes_for_plot) > 2:
                try:
                    sizes_array = np.array(sizes_for_plot)
                    tcs_array = np.array(tcs_for_plot)
                    
                    # Fit Tc(L) = Tc_inf + a/L
                    inv_sizes = 1.0 / sizes_array
                    coeffs = np.polyfit(inv_sizes, tcs_array, 1)
                    tc_inf = coeffs[1]  # y-intercept
                    
                    # Plot extrapolation
                    x_extrap = np.linspace(0, inv_sizes.max(), 100)
                    y_extrap = coeffs[0] * x_extrap + tc_inf
                    size_extrap = 1.0 / (x_extrap + 1e-10)
                    
                    valid_extrap = x_extrap > 1e-4
                    ax3.plot(size_extrap[valid_extrap], y_extrap[valid_extrap], '--', 
                            color=layer_colors[layer-1], alpha=0.7,
                            label=f'Layer {layer}: Tc(∞) = {tc_inf:.3f}')
                    
                    print(f"  → Finite-size extrapolation: Tc(∞) = {tc_inf:.4f}")
                    
                except Exception as e:
                    print(f"  → Extrapolation failed: {e}")
        
        # Plot peak heights vs size
        if len(sizes_for_plot) > 1:
            ax4.plot(sizes_for_plot, peak_cvs_for_plot, 'o-', 
                    color=layer_colors[layer-1], linewidth=2, markersize=8,
                    label=f'Layer {layer}')
    
    # Finalize plots
    ax1.set_xlim(0.5, 2.5)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlim(0.8, 2.0)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('specific_heat_critical_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved: specific_heat_critical_analysis.png")
    plt.show()
    
    return peak_analysis

def analyze_critical_scaling():
    """Analyze critical scaling behavior of specific heat."""
    layer_data = load_simulation_data()
    lattice_sizes = sorted(layer_data.keys())
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    
    print(f"\n=== CRITICAL SCALING ANALYSIS ===")
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for layer in [1, 2, 3]:
        ax = axes[layer-1]
        
        print(f"\n--- Layer {layer} Scaling ---")
        
        # Collect peak data
        sizes = []
        peak_cvs = []
        
        for size in lattice_sizes:
            if size not in layer_data or layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            
            # Filter reasonable values
            mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 2.5)
            mask &= (data['Cv'] > 0) & (data['Cv'] < 10000)
            
            if mask.sum() > 10:
                cv_filtered = data['Cv'][mask]
                peak_cv = cv_filtered.max()
                
                sizes.append(size)
                peak_cvs.append(peak_cv)
        
        if len(sizes) > 2:
            sizes = np.array(sizes)
            peak_cvs = np.array(peak_cvs)
            
            # Plot log-log
            ax.loglog(sizes, peak_cvs, 'o-', color=layer_colors[layer-1], 
                     linewidth=2, markersize=8, label=f'Layer {layer} Data')
            
            # Fit power law: Cv_max ~ L^α/ν
            try:
                log_sizes = np.log(sizes)
                log_cvs = np.log(peak_cvs)
                
                coeffs = np.polyfit(log_sizes, log_cvs, 1)
                scaling_exp = coeffs[0]
                
                # Plot fit
                size_fit = np.logspace(np.log10(sizes.min()), np.log10(sizes.max()), 100)
                cv_fit = np.exp(coeffs[1]) * size_fit**scaling_exp
                
                ax.loglog(size_fit, cv_fit, '--', color=layer_colors[layer-1], 
                         alpha=0.7, label=f'Fit: $C_v ∝ L^{{{scaling_exp:.2f}}}$')
                
                print(f"  Peak Cv scaling: Cv_max ∝ L^{scaling_exp:.3f}")
                
                # Compare with theoretical expectations
                print(f"  Expected for 2D Ising: α/ν = 0 (logarithmic)")
                print(f"  Expected for 2D Potts: α/ν ~ 0.5-1.0")
                
            except Exception as e:
                print(f"  Scaling fit failed: {e}")
        
        ax.set_xlabel('Lattice Size L')
        ax.set_ylabel('Peak Specific Heat')
        ax.set_title(f'Layer {layer}: Cv Peak Scaling')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig('specific_heat_scaling_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved: specific_heat_scaling_analysis.png")
    plt.show()

def detailed_critical_behavior():
    """Examine detailed critical behavior near Tc."""
    layer_data = load_simulation_data()
    lattice_sizes = sorted(layer_data.keys())
    
    # Use largest system for detailed analysis
    if not lattice_sizes:
        return
        
    size = lattice_sizes[-1]
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    
    print(f"\n=== DETAILED CRITICAL BEHAVIOR (L={size}) ===")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    for layer in [1, 2, 3]:
        if size not in layer_data or layer not in layer_data[size]:
            continue
            
        data = layer_data[size][layer]
        
        # Find approximate Tc from peak
        mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 2.5)
        mask &= (data['Cv'] > 0) & (data['Cv'] < 10000)
        
        if mask.sum() < 10:
            continue
            
        temp_filt = data['Temperature'][mask]
        cv_filt = data['Cv'][mask]
        
        peak_idx = np.argmax(cv_filt)
        tc_approx = temp_filt.iloc[peak_idx]
        
        print(f"\nLayer {layer}: Approximate Tc = {tc_approx:.3f}")
        
        # Critical region analysis
        critical_window = 0.4
        crit_mask = (data['Temperature'] >= tc_approx - critical_window) & \
                   (data['Temperature'] <= tc_approx + critical_window)
        
        temp_crit = data['Temperature'][crit_mask]
        cv_crit = data['Cv'][crit_mask]
        m2_crit = data['M2'][crit_mask]
        energy_crit = data['Energy'][crit_mask]
        
        # Plot 1: Specific heat in critical region
        ax1 = axes[0, 0]
        reasonable_cv = cv_crit < 2000
        if reasonable_cv.sum() > 5:
            ax1.plot(temp_crit[reasonable_cv], cv_crit[reasonable_cv], 'o-', 
                    color=layer_colors[layer-1], alpha=0.8, linewidth=2,
                    label=f'Layer {layer}')
            ax1.axvline(x=tc_approx, color=layer_colors[layer-1], 
                       linestyle='--', alpha=0.5)
        
        # Plot 2: Magnetization
        ax2 = axes[0, 1]
        ax2.plot(temp_crit, m2_crit, 'o-', color=layer_colors[layer-1], 
                alpha=0.8, linewidth=2, label=f'Layer {layer}')
        ax2.axvline(x=tc_approx, color=layer_colors[layer-1], 
                   linestyle='--', alpha=0.5)
        
        # Plot 3: Energy
        ax3 = axes[1, 0]
        ax3.plot(temp_crit, energy_crit, 'o-', color=layer_colors[layer-1], 
                alpha=0.8, linewidth=2, label=f'Layer {layer}')
        ax3.axvline(x=tc_approx, color=layer_colors[layer-1], 
                   linestyle='--', alpha=0.5)
        
        # Plot 4: Energy derivative (dE/dT)
        ax4 = axes[1, 1]
        if len(temp_crit) > 5:
            de_dt = np.gradient(energy_crit, temp_crit)
            ax4.plot(temp_crit, de_dt, 'o-', color=layer_colors[layer-1], 
                    alpha=0.8, linewidth=2, label=f'Layer {layer}')
            ax4.axvline(x=tc_approx, color=layer_colors[layer-1], 
                       linestyle='--', alpha=0.5)
    
    # Set labels and formatting
    axes[0, 0].set_xlabel('Temperature T')
    axes[0, 0].set_ylabel('Specific Heat $C_v$')
    axes[0, 0].set_title('Critical Region: Specific Heat')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    
    axes[0, 1].set_xlabel('Temperature T')
    axes[0, 1].set_ylabel('Magnetization $\\langle M^2 \\rangle$')
    axes[0, 1].set_title('Critical Region: Magnetization')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    
    axes[1, 0].set_xlabel('Temperature T')
    axes[1, 0].set_ylabel('Energy per site')
    axes[1, 0].set_title('Critical Region: Energy')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()
    
    axes[1, 1].set_xlabel('Temperature T')
    axes[1, 1].set_ylabel('$dE/dT$')
    axes[1, 1].set_title('Critical Region: Energy Derivative')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig('detailed_critical_behavior.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved: detailed_critical_behavior.png")
    plt.show()

def main():
    """Main analysis function."""
    print("=== COMPREHENSIVE SPECIFIC HEAT CRITICAL ANALYSIS ===")
    print("Focus on heat capacity peaks and critical scaling")
    print("="*60)
    
    # Run all analyses
    peak_data = analyze_specific_heat_peaks()
    analyze_critical_scaling()
    detailed_critical_behavior()
    
    print(f"\n{'='*60}")
    print("SPECIFIC HEAT ANALYSIS COMPLETE")
    print("="*60)
    print("Generated files:")
    print("- specific_heat_critical_analysis.png")
    print("- specific_heat_scaling_analysis.png")
    print("- detailed_critical_behavior.png")
    
    return peak_data

if __name__ == "__main__":
    main()