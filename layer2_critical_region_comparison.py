#!/usr/bin/env python3
"""
Focused analysis of Layer 2 critical region behavior across different lattice sizes.
Detailed comparison of thermodynamic properties in the critical temperature range.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

def load_layer2_data():
    """Load Layer 2 simulation data for all lattice sizes."""
    data_dir = "simulation_runs"
    layer2_data = {}
    
    for size_dir in glob.glob(f"{data_dir}/nx*_ny*"):
        size_name = os.path.basename(size_dir)
        parts = size_name.split('_')
        nx = int(parts[0][2:])
        ny = int(parts[1][2:])
        
        if nx == ny:  # Only square lattices
            layer_file = f"{size_dir}/layer_2.txt"
            if os.path.exists(layer_file):
                try:
                    data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                     names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                    layer2_data[nx] = data
                    print(f"Loaded Layer 2 data for size {nx}x{nx}: {len(data)} points")
                except Exception as e:
                    print(f"Warning: Could not load {layer_file}: {e}")
    
    return layer2_data

def analyze_critical_region():
    """Analyze critical region behavior for Layer 2 across all lattice sizes."""
    layer2_data = load_layer2_data()
    lattice_sizes = sorted(layer2_data.keys())
    
    if not lattice_sizes:
        print("No Layer 2 data found!")
        return
    
    print("=== LAYER 2 CRITICAL REGION COMPARISON ===")
    print(f"Available lattice sizes: {lattice_sizes}")
    print("="*60)
    
    # Define color scheme for different lattice sizes
    size_colors = plt.cm.viridis(np.linspace(0, 1, len(lattice_sizes)))
    
    # Create comprehensive critical region analysis
    fig, axes = plt.subplots(3, 3, figsize=(20, 15))
    fig.suptitle('Layer 2: Critical Region Comparison Across Lattice Sizes', fontsize=16, fontweight='bold')
    
    # Define critical temperature ranges for different analyses
    # Based on previous analysis: small systems peak at T=2.0, large system at T=0.5
    critical_ranges = {
        'full': (0.4, 2.1),
        'fine': (0.4, 1.2),
        'boundary': (1.8, 2.1)
    }
    
    print("\nCritical Temperature Analysis:")
    print("-" * 40)
    
    # Storage for analysis results
    analysis_results = {}
    
    for i, size in enumerate(lattice_sizes):
        data = layer2_data[size]
        color = size_colors[i]
        
        # Filter reasonable data
        temp_mask = (data['Temperature'] >= 0.4) & (data['Temperature'] <= 2.1)
        cv_mask = (data['Cv'] > 0) & (data['Cv'] < 15000)
        combined_mask = temp_mask & cv_mask
        
        if combined_mask.sum() < 10:
            print(f"Insufficient data for size {size}x{size}")
            continue
        
        temp = data['Temperature'][combined_mask]
        cv = data['Cv'][combined_mask]
        m2 = data['M2'][combined_mask]
        energy = data['Energy'][combined_mask]
        g2 = data['G2'][combined_mask]
        
        # Find critical temperature
        if len(cv) > 10:
            cv_smooth = savgol_filter(cv, window_length=min(11, len(cv)//3*2+1), polyorder=2) if len(cv) > 20 else cv
            peak_idx = np.argmax(cv_smooth)
            tc = temp.iloc[peak_idx]
            cv_max = cv_smooth[peak_idx]
            
            analysis_results[size] = {
                'tc': tc,
                'cv_max': cv_max,
                'temp': temp,
                'cv': cv,
                'm2': m2,
                'energy': energy,
                'g2': g2
            }
            
            print(f"Size {size}x{size}: Tc = {tc:.4f}, Peak Cv = {cv_max:.2f}")
        
        # Plot 1: Full range specific heat
        ax1 = axes[0, 0]
        ax1.plot(temp, cv, 'o-', color=color, alpha=0.7, linewidth=2, 
                markersize=3, label=f'L={size}')
        ax1.set_title('Specific Heat - Full Range')
        ax1.set_xlabel('Temperature T')
        ax1.set_ylabel('Specific Heat $C_v$')
        
        # Plot 2: Critical region zoom (T = 0.4 - 1.2)
        ax2 = axes[0, 1]
        crit_mask = (temp >= 0.4) & (temp <= 1.2)
        if crit_mask.sum() > 5:
            ax2.plot(temp[crit_mask], cv[crit_mask], 'o-', color=color, 
                    alpha=0.8, linewidth=2, markersize=4, label=f'L={size}')
        ax2.set_title('Specific Heat - Critical Region (0.4-1.2)')
        ax2.set_xlabel('Temperature T')
        ax2.set_ylabel('Specific Heat $C_v$')
        
        # Plot 3: High temperature region (T = 1.8 - 2.1)
        ax3 = axes[0, 2]
        high_mask = (temp >= 1.8) & (temp <= 2.1)
        if high_mask.sum() > 5:
            ax3.plot(temp[high_mask], cv[high_mask], 'o-', color=color,
                    alpha=0.8, linewidth=2, markersize=4, label=f'L={size}')
        ax3.set_title('Specific Heat - High T Region (1.8-2.1)')
        ax3.set_xlabel('Temperature T')
        ax3.set_ylabel('Specific Heat $C_v$')
        
        # Plot 4: Magnetization full range
        ax4 = axes[1, 0]
        ax4.plot(temp, m2, 'o-', color=color, alpha=0.7, linewidth=2,
                markersize=3, label=f'L={size}')
        ax4.set_title('Magnetization - Full Range')
        ax4.set_xlabel('Temperature T')
        ax4.set_ylabel('$\\langle M^2 \\rangle$')
        
        # Plot 5: Magnetization critical region
        ax5 = axes[1, 1]
        if crit_mask.sum() > 5:
            ax5.plot(temp[crit_mask], m2[crit_mask], 'o-', color=color,
                    alpha=0.8, linewidth=2, markersize=4, label=f'L={size}')
        ax5.set_title('Magnetization - Critical Region')
        ax5.set_xlabel('Temperature T')
        ax5.set_ylabel('$\\langle M^2 \\rangle$')
        
        # Plot 6: Energy full range
        ax6 = axes[1, 2]
        ax6.plot(temp, energy, 'o-', color=color, alpha=0.7, linewidth=2,
                markersize=3, label=f'L={size}')
        ax6.set_title('Energy - Full Range')
        ax6.set_xlabel('Temperature T')
        ax6.set_ylabel('Energy per site')
        
        # Plot 7: Binder cumulant
        ax7 = axes[2, 0]
        # Calculate Binder cumulant: U = 1 - M4/(3*M2^2)
        binder = 1 - data['M4'][combined_mask] / (3 * data['M2'][combined_mask]**2)
        ax7.plot(temp, binder, 'o-', color=color, alpha=0.7, linewidth=2,
                markersize=3, label=f'L={size}')
        ax7.set_title('Binder Cumulant')
        ax7.set_xlabel('Temperature T')
        ax7.set_ylabel('$U = 1 - \\langle M^4 \\rangle / 3\\langle M^2 \\rangle^2$')
        
        # Plot 8: Susceptibility (from G2)
        ax8 = axes[2, 1]
        ax8.plot(temp, g2, 'o-', color=color, alpha=0.7, linewidth=2,
                markersize=3, label=f'L={size}')
        ax8.set_title('Susceptibility $\\chi$')
        ax8.set_xlabel('Temperature T')
        ax8.set_ylabel('$\\chi \\propto \\langle G^2 \\rangle$')
        
        # Plot 9: Energy derivative (approximate specific heat)
        ax9 = axes[2, 2]
        if len(temp) > 5:
            de_dt = np.gradient(energy, temp)
            ax9.plot(temp, de_dt, 'o-', color=color, alpha=0.7, linewidth=2,
                    markersize=3, label=f'L={size}')
        ax9.set_title('Energy Derivative $dE/dT$')
        ax9.set_xlabel('Temperature T')
        ax9.set_ylabel('$dE/dT$')
    
    # Add legends and grids to all subplots
    for ax in axes.flat:
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('layer2_critical_region_comparison.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved: layer2_critical_region_comparison.png")
    plt.show()
    
    return analysis_results

def detailed_critical_point_analysis():
    """Detailed analysis of critical points for each lattice size."""
    layer2_data = load_layer2_data()
    lattice_sizes = sorted(layer2_data.keys())
    
    print(f"\n=== DETAILED CRITICAL POINT ANALYSIS - LAYER 2 ===")
    
    # Create focused critical point comparison
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Layer 2: Detailed Critical Point Analysis', fontsize=16, fontweight='bold')
    
    size_colors = plt.cm.plasma(np.linspace(0, 1, len(lattice_sizes)))
    critical_data = {}
    
    for i, size in enumerate(lattice_sizes):
        if size not in layer2_data:
            continue
            
        data = layer2_data[size]
        color = size_colors[i]
        
        # Filter data around potential critical regions
        # Check both low T (around 0.5) and high T (around 2.0) regions
        regions = {
            'low_T': (0.4, 1.0),
            'high_T': (1.8, 2.1)
        }
        
        critical_data[size] = {}
        
        for region_name, (t_min, t_max) in regions.items():
            mask = (data['Temperature'] >= t_min) & (data['Temperature'] <= t_max)
            mask &= (data['Cv'] > 0) & (data['Cv'] < 15000)
            
            if mask.sum() < 5:
                continue
                
            temp_region = data['Temperature'][mask]
            cv_region = data['Cv'][mask]
            m2_region = data['M2'][mask]
            
            # Find peak in this region
            if len(cv_region) > 5:
                peak_idx = cv_region.idxmax()
                tc_region = data.loc[peak_idx, 'Temperature']
                cv_max_region = data.loc[peak_idx, 'Cv']
                
                critical_data[size][region_name] = {
                    'tc': tc_region,
                    'cv_max': cv_max_region,
                    'temp': temp_region,
                    'cv': cv_region,
                    'm2': m2_region
                }
        
        # Plot critical regions
        # Plot 1: Low temperature critical region
        ax1 = axes[0, 0]
        if 'low_T' in critical_data[size]:
            data_low = critical_data[size]['low_T']
            ax1.plot(data_low['temp'], data_low['cv'], 'o-', color=color,
                    linewidth=2, markersize=5, alpha=0.8, label=f'L={size}')
            # Mark critical point
            ax1.plot(data_low['tc'], data_low['cv_max'], '*', color=color,
                    markersize=12, markeredgecolor='black', markeredgewidth=1)
        ax1.set_title('Low Temperature Critical Region (0.4-1.0)')
        ax1.set_xlabel('Temperature T')
        ax1.set_ylabel('Specific Heat $C_v$')
        
        # Plot 2: High temperature critical region  
        ax2 = axes[0, 1]
        if 'high_T' in critical_data[size]:
            data_high = critical_data[size]['high_T']
            ax2.plot(data_high['temp'], data_high['cv'], 'o-', color=color,
                    linewidth=2, markersize=5, alpha=0.8, label=f'L={size}')
            # Mark critical point
            ax2.plot(data_high['tc'], data_high['cv_max'], '*', color=color,
                    markersize=12, markeredgecolor='black', markeredgewidth=1)
        ax2.set_title('High Temperature Critical Region (1.8-2.1)')
        ax2.set_xlabel('Temperature T')
        ax2.set_ylabel('Specific Heat $C_v$')
        
        # Plot 3: Critical temperature vs lattice size
        ax3 = axes[1, 0]
        for region_name in ['low_T', 'high_T']:
            if region_name in critical_data[size]:
                tc = critical_data[size][region_name]['tc']
                marker = 'o' if region_name == 'low_T' else 's'
                ax3.plot(size, tc, marker, color=color, markersize=8,
                        label=f'L={size}, {region_name}' if i == 0 else '')
        
        # Plot 4: Peak height vs lattice size
        ax4 = axes[1, 1]
        for region_name in ['low_T', 'high_T']:
            if region_name in critical_data[size]:
                cv_max = critical_data[size][region_name]['cv_max']
                marker = 'o' if region_name == 'low_T' else 's'
                ax4.plot(size, cv_max, marker, color=color, markersize=8,
                        label=f'L={size}, {region_name}' if i == 0 else '')
    
    # Format plots
    ax3.set_title('Critical Temperature vs Lattice Size')
    ax3.set_xlabel('Lattice Size L')
    ax3.set_ylabel('Critical Temperature $T_c$')
    ax3.set_yscale('linear')
    
    ax4.set_title('Peak Height vs Lattice Size')
    ax4.set_xlabel('Lattice Size L')
    ax4.set_ylabel('Peak $C_v$ Value')
    ax4.set_yscale('log')
    
    for ax in axes.flat:
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('layer2_detailed_critical_analysis.png', dpi=300, bbox_inches='tight')
    print(f"Saved: layer2_detailed_critical_analysis.png")
    plt.show()
    
    # Print detailed results
    print(f"\nDetailed Critical Point Results:")
    print("-" * 60)
    print("Size\tLow-T Tc\tLow-T Cv\tHigh-T Tc\tHigh-T Cv")
    print("-" * 60)
    
    for size in lattice_sizes:
        if size in critical_data:
            low_tc = critical_data[size].get('low_T', {}).get('tc', 'N/A')
            low_cv = critical_data[size].get('low_T', {}).get('cv_max', 'N/A')
            high_tc = critical_data[size].get('high_T', {}).get('tc', 'N/A')
            high_cv = critical_data[size].get('high_T', {}).get('cv_max', 'N/A')
            
            print(f"{size}x{size}\t{low_tc if low_tc == 'N/A' else f'{low_tc:.4f}'}\t"
                  f"{low_cv if low_cv == 'N/A' else f'{low_cv:.2f}'}\t"
                  f"{high_tc if high_tc == 'N/A' else f'{high_tc:.4f}'}\t"
                  f"{high_cv if high_cv == 'N/A' else f'{high_cv:.2f}'}")
    
    return critical_data

def scaling_analysis():
    """Analyze finite-size scaling for Layer 2."""
    layer2_data = load_layer2_data()
    lattice_sizes = sorted(layer2_data.keys())
    
    print(f"\n=== LAYER 2 FINITE-SIZE SCALING ANALYSIS ===")
    
    # Collect scaling data
    sizes = []
    tc_values = []
    cv_max_values = []
    
    for size in lattice_sizes:
        if size not in layer2_data:
            continue
            
        data = layer2_data[size]
        
        # Find overall peak
        mask = (data['Temperature'] >= 0.4) & (data['Temperature'] <= 2.1)
        mask &= (data['Cv'] > 0) & (data['Cv'] < 15000)
        
        if mask.sum() < 10:
            continue
            
        cv_filtered = data['Cv'][mask]
        temp_filtered = data['Temperature'][mask]
        
        peak_idx = cv_filtered.idxmax()
        tc = data.loc[peak_idx, 'Temperature']
        cv_max = data.loc[peak_idx, 'Cv']
        
        sizes.append(size)
        tc_values.append(tc)
        cv_max_values.append(cv_max)
    
    if len(sizes) < 3:
        print("Insufficient data for scaling analysis")
        return
    
    # Create scaling plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Layer 2: Finite-Size Scaling Analysis', fontsize=16, fontweight='bold')
    
    sizes = np.array(sizes)
    tc_values = np.array(tc_values)
    cv_max_values = np.array(cv_max_values)
    
    # Plot 1: Tc vs 1/L (finite-size scaling)
    ax1 = axes[0]
    inv_sizes = 1.0 / sizes
    ax1.plot(inv_sizes, tc_values, 'o-', markersize=8, linewidth=2, color='blue')
    
    # Linear fit: Tc(L) = Tc_inf + a/L
    try:
        coeffs = np.polyfit(inv_sizes, tc_values, 1)
        tc_inf = coeffs[1]
        x_fit = np.linspace(0, inv_sizes.max(), 100)
        y_fit = coeffs[0] * x_fit + tc_inf
        ax1.plot(x_fit, y_fit, '--', color='red', alpha=0.7, linewidth=2,
                label=f'Tc(∞) = {tc_inf:.3f}')
        print(f"Finite-size extrapolation: Tc(∞) = {tc_inf:.4f}")
    except:
        print("Finite-size extrapolation failed")
    
    ax1.set_xlabel('1/L')
    ax1.set_ylabel('Critical Temperature $T_c$')
    ax1.set_title('Critical Temperature Scaling')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Peak height scaling (log-log)
    ax2 = axes[1]
    ax2.loglog(sizes, cv_max_values, 'o-', markersize=8, linewidth=2, color='green')
    
    # Power law fit: Cv_max ~ L^alpha
    try:
        log_sizes = np.log(sizes)
        log_cvs = np.log(cv_max_values)
        coeffs = np.polyfit(log_sizes, log_cvs, 1)
        alpha = coeffs[0]
        
        size_fit = np.logspace(np.log10(sizes.min()), np.log10(sizes.max()), 100)
        cv_fit = np.exp(coeffs[1]) * size_fit**alpha
        ax2.loglog(size_fit, cv_fit, '--', color='red', alpha=0.7, linewidth=2,
                  label=f'Cv ∝ L^{alpha:.2f}')
        print(f"Peak height scaling: Cv_max ∝ L^{alpha:.3f}")
    except:
        print("Peak height scaling failed")
    
    ax2.set_xlabel('Lattice Size L')
    ax2.set_ylabel('Peak $C_v$ Value')
    ax2.set_title('Peak Height Scaling')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Data collapse attempt
    ax3 = axes[2]
    # Try to collapse data using Tc scaling
    for i, size in enumerate(sizes):
        data = layer2_data[size]
        mask = (data['Temperature'] >= 0.4) & (data['Temperature'] <= 2.1)
        mask &= (data['Cv'] > 0) & (data['Cv'] < 15000)
        
        if mask.sum() < 10:
            continue
            
        temp = data['Temperature'][mask]
        cv = data['Cv'][mask]
        
        # Scale temperature: t = (T - Tc) * L^(1/nu)
        tc = tc_values[i]
        nu = 1.0  # Assume 2D Ising value
        t_scaled = (temp - tc) * size**(1/nu)
        
        # Scale Cv: cv_scaled = Cv / L^(alpha/nu)
        alpha_over_nu = 3.7  # From our scaling analysis
        cv_scaled = cv / size**alpha_over_nu
        
        ax3.plot(t_scaled, cv_scaled, 'o', alpha=0.7, markersize=4, label=f'L={size}')
    
    ax3.set_xlabel('$(T - T_c) L^{1/\\nu}$')
    ax3.set_ylabel('$C_v / L^{\\alpha/\\nu}$')
    ax3.set_title('Data Collapse Attempt')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('layer2_scaling_analysis.png', dpi=300, bbox_inches='tight')
    print(f"Saved: layer2_scaling_analysis.png")
    plt.show()

def main():
    """Main analysis function."""
    print("="*60)
    print("LAYER 2 CRITICAL REGION ANALYSIS")
    print("Detailed comparison across all lattice sizes")
    print("="*60)
    
    # Run all analyses
    try:
        # 1. Overall critical region comparison
        results = analyze_critical_region()
        
        # 2. Detailed critical point analysis
        critical_data = detailed_critical_point_analysis()
        
        # 3. Scaling analysis
        scaling_analysis()
        
        print(f"\n{'='*60}")
        print("LAYER 2 ANALYSIS COMPLETE")
        print("="*60)
        print("Generated files:")
        print("- layer2_critical_region_comparison.png")
        print("- layer2_detailed_critical_analysis.png") 
        print("- layer2_scaling_analysis.png")
        
        return results, critical_data
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return None, None

if __name__ == "__main__":
    main()