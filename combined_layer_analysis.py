#!/usr/bin/env python3
"""
Comprehensive combined analysis of all three layers across all lattice sizes.
Compares critical temperatures, thermodynamic properties, and scaling behavior
between layers and across different system sizes.
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

def load_all_layer_data():
    """Load all simulation data for all layers and lattice sizes."""
    data_dir = "simulation_runs"
    all_data = {}
    
    for size_dir in glob.glob(f"{data_dir}/nx*_ny*"):
        size_name = os.path.basename(size_dir)
        parts = size_name.split('_')
        nx = int(parts[0][2:])
        ny = int(parts[1][2:])
        
        if nx == ny:  # Only square lattices
            all_data[nx] = {}
            for layer in [1, 2, 3]:
                layer_file = f"{size_dir}/layer_{layer}.txt"
                if os.path.exists(layer_file):
                    try:
                        data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                         names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                        all_data[nx][layer] = data
                    except Exception as e:
                        print(f"Warning: Could not load {layer_file}: {e}")
    
    return all_data

def analyze_critical_temperatures():
    """Analyze and compare critical temperatures across all layers and sizes."""
    all_data = load_all_layer_data()
    lattice_sizes = sorted(all_data.keys())
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    layer_names = ['Layer 1', 'Layer 2', 'Layer 3']
    
    print("=== COMBINED LAYER CRITICAL TEMPERATURE ANALYSIS ===")
    print("Comparing critical temperatures across all layers and lattice sizes")
    print("="*70)
    
    # Storage for critical temperature data
    critical_data = {}
    for layer in [1, 2, 3]:
        critical_data[layer] = {'sizes': [], 'tc_values': [], 'cv_peaks': []}
    
    # Create comprehensive comparison plot
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    
    # Plot 1: Specific heat comparison for all layers and sizes
    ax1 = axes[0, 0]
    ax1.set_title('Specific Heat: All Layers & Sizes Combined', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Temperature T')
    ax1.set_ylabel('Specific Heat $C_v$')
    
    # Plot 2: Critical temperature vs lattice size
    ax2 = axes[0, 1]
    ax2.set_title('Critical Temperature vs Lattice Size', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Lattice Size L')
    ax2.set_ylabel('Critical Temperature $T_c$')
    
    # Plot 3: Peak height comparison
    ax3 = axes[0, 2]
    ax3.set_title('Peak Cv Height vs Lattice Size', fontsize=14, fontweight='bold')
    ax3.set_xlabel('Lattice Size L')
    ax3.set_ylabel('Peak $C_v$ Value')
    ax3.set_yscale('log')
    
    # Plot 4: Layer comparison at each size
    ax4 = axes[1, 0]
    ax4.set_title('Layer Comparison: Critical Temperatures', fontsize=14, fontweight='bold')
    ax4.set_xlabel('Lattice Size L')
    ax4.set_ylabel('$T_c$ Difference from Layer 1')
    
    # Plot 5: Finite-size scaling comparison
    ax5 = axes[1, 1]
    ax5.set_title('Finite-Size Scaling: Tc Extrapolation', fontsize=14, fontweight='bold')
    ax5.set_xlabel('1/L')
    ax5.set_ylabel('Critical Temperature $T_c$')
    
    # Plot 6: Scaling exponents comparison
    ax6 = axes[1, 2]
    ax6.set_title('Scaling Exponents: Cv Peak Heights', fontsize=14, fontweight='bold')
    ax6.set_xlabel('Lattice Size L')
    ax6.set_ylabel('Peak $C_v$ Value')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    
    print("\nDetailed Critical Temperature Analysis:")
    print("-" * 50)
    
    # Process each layer
    for layer in [1, 2, 3]:
        print(f"\n=== LAYER {layer} ===")
        
        for size in lattice_sizes:
            if size not in all_data or layer not in all_data[size]:
                continue
                
            data = all_data[size][layer]
            
            # Filter reasonable data
            temp_mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.5)
            cv_mask = (data['Cv'] > 0) & (data['Cv'] < 15000)  # Remove artifacts
            combined_mask = temp_mask & cv_mask
            
            if combined_mask.sum() < 10:
                continue
            
            temp = data['Temperature'][combined_mask].values
            cv = data['Cv'][combined_mask].values
            
            # Find critical temperature from peak
            if len(cv) > 10:
                # Smooth for peak finding
                if len(cv) > 20:
                    cv_smooth = savgol_filter(cv, window_length=min(11, len(cv)//3*2+1), polyorder=2)
                else:
                    cv_smooth = cv
                
                peak_idx = np.argmax(cv_smooth)
                tc = temp[peak_idx]
                cv_peak = cv_smooth[peak_idx]
                
                # Store data
                critical_data[layer]['sizes'].append(size)
                critical_data[layer]['tc_values'].append(tc)
                critical_data[layer]['cv_peaks'].append(cv_peak)
                
                # Plot specific heat
                alpha = 0.5 + 0.1 * (lattice_sizes.index(size))
                linestyle = '-' if layer == 1 else ('--' if layer == 2 else ':')
                linewidth = 1.5 + 0.5 * (lattice_sizes.index(size))
                
                ax1.plot(temp, cv, color=layer_colors[layer-1], alpha=alpha,
                        linestyle=linestyle, linewidth=linewidth,
                        label=f'L={size}, {layer_names[layer-1]}' if size == lattice_sizes[-1] else '')
                
                print(f"  Size {size}x{size}: Tc = {tc:.4f}, Peak Cv = {cv_peak:.2f}")
        
        # Convert to arrays for analysis
        if critical_data[layer]['sizes']:
            sizes_array = np.array(critical_data[layer]['sizes'])
            tc_array = np.array(critical_data[layer]['tc_values'])
            cv_array = np.array(critical_data[layer]['cv_peaks'])
            
            # Plot critical temperatures vs size
            ax2.plot(sizes_array, tc_array, 'o-', color=layer_colors[layer-1],
                    linewidth=2, markersize=8, label=layer_names[layer-1])
            
            # Plot peak heights vs size
            ax3.plot(sizes_array, cv_array, 'o-', color=layer_colors[layer-1],
                    linewidth=2, markersize=8, label=layer_names[layer-1])
            
            # Finite-size scaling analysis
            if len(sizes_array) > 2:
                try:
                    inv_sizes = 1.0 / sizes_array
                    
                    # Linear fit: Tc(L) = Tc_inf + a/L
                    coeffs = np.polyfit(inv_sizes, tc_array, 1)
                    tc_inf = coeffs[1]
                    
                    # Plot extrapolation
                    inv_range = np.linspace(0, inv_sizes.max(), 100)
                    tc_extrap = coeffs[0] * inv_range + tc_inf
                    
                    ax5.plot(inv_range, tc_extrap, '--', color=layer_colors[layer-1],
                            alpha=0.7, linewidth=2)
                    ax5.plot(inv_sizes, tc_array, 'o', color=layer_colors[layer-1],
                            markersize=8, label=f'{layer_names[layer-1]}: Tc(∞)={tc_inf:.3f}')
                    
                    print(f"  Finite-size extrapolation: Tc(∞) = {tc_inf:.4f}")
                    
                except Exception as e:
                    print(f"  Extrapolation failed: {e}")
            
            # Scaling analysis for peak heights
            if len(sizes_array) > 2:
                try:
                    # Power law fit: Cv_max ~ L^alpha
                    log_sizes = np.log(sizes_array)
                    log_cvs = np.log(cv_array)
                    
                    coeffs = np.polyfit(log_sizes, log_cvs, 1)
                    scaling_exp = coeffs[0]
                    
                    # Plot scaling
                    size_fit = np.logspace(np.log10(sizes_array.min()), np.log10(sizes_array.max()), 100)
                    cv_fit = np.exp(coeffs[1]) * size_fit**scaling_exp
                    
                    ax6.plot(size_fit, cv_fit, '--', color=layer_colors[layer-1],
                            alpha=0.7, linewidth=2)
                    ax6.plot(sizes_array, cv_array, 'o', color=layer_colors[layer-1],
                            markersize=8, label=f'{layer_names[layer-1]}: α={scaling_exp:.2f}')
                    
                    print(f"  Peak scaling: Cv_max ∝ L^{scaling_exp:.3f}")
                    
                except Exception as e:
                    print(f"  Scaling analysis failed: {e}")
    
    # Layer comparison analysis
    print(f"\n=== LAYER COMPARISON ANALYSIS ===")
    
    # Compare critical temperatures between layers
    common_sizes = set(critical_data[1]['sizes']) & set(critical_data[2]['sizes']) & set(critical_data[3]['sizes'])
    common_sizes = sorted(list(common_sizes))
    
    if len(common_sizes) > 1:
        print(f"\nCritical Temperature Comparison (Common sizes: {common_sizes}):")
        print("Size\tLayer1\tLayer2\tLayer3\tΔ(2-1)\tΔ(3-1)")
        print("-" * 60)
        
        for size in common_sizes:
            # Find Tc values for this size
            tc_values = {}
            for layer in [1, 2, 3]:
                try:
                    idx = critical_data[layer]['sizes'].index(size)
                    tc_values[layer] = critical_data[layer]['tc_values'][idx]
                except ValueError:
                    tc_values[layer] = None
            
            if all(tc_values[layer] is not None for layer in [1, 2, 3]):
                diff_21 = tc_values[2] - tc_values[1]
                diff_31 = tc_values[3] - tc_values[1]
                
                print(f"{size}\t{tc_values[1]:.4f}\t{tc_values[2]:.4f}\t{tc_values[3]:.4f}\t{diff_21:+.4f}\t{diff_31:+.4f}")
                
                # Plot differences
                ax4.plot(size, diff_21, 'o', color=layer_colors[1], markersize=8, label='Layer 2 - Layer 1' if size == common_sizes[0] else '')
                ax4.plot(size, diff_31, 's', color=layer_colors[2], markersize=8, label='Layer 3 - Layer 1' if size == common_sizes[0] else '')
    
    # Finalize all plots
    for ax in axes.flat:
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    ax4.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    ax4.set_ylabel('$ΔT_c$ (Layer - Layer 1)')
    
    plt.tight_layout()
    plt.savefig('combined_layer_critical_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved: combined_layer_critical_analysis.png")
    plt.show()
    
    return critical_data

def analyze_thermodynamic_comparison():
    """Compare all thermodynamic properties across layers and sizes."""
    all_data = load_all_layer_data()
    lattice_sizes = sorted(all_data.keys())
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    layer_names = ['Layer 1', 'Layer 2', 'Layer 3']
    
    print(f"\n=== THERMODYNAMIC PROPERTIES COMPARISON ===")
    
    # Use the largest system for detailed comparison
    if not lattice_sizes:
        print("No data available!")
        return
    
    size = lattice_sizes[-1]
    print(f"Detailed analysis using largest system: L={size}")
    
    # Create comprehensive thermodynamic comparison
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    
    properties = ['Cv', 'M2', 'Energy']
    property_labels = ['Specific Heat $C_v$', 'Magnetization $\\langle M^2 \\rangle$', 'Energy per site']
    
    for prop_idx, (prop, prop_label) in enumerate(zip(properties, property_labels)):
        # Full temperature range
        ax_full = axes[prop_idx, 0]
        ax_full.set_title(f'{prop_label} - Full Range (L={size})', fontweight='bold')
        ax_full.set_xlabel('Temperature T')
        ax_full.set_ylabel(prop_label)
        
        # Critical region
        ax_crit = axes[prop_idx, 1]
        ax_crit.set_title(f'{prop_label} - Critical Region', fontweight='bold')
        ax_crit.set_xlabel('Temperature T')
        ax_crit.set_ylabel(prop_label)
        
        # Layer differences
        ax_diff = axes[prop_idx, 2]
        ax_diff.set_title(f'{prop_label} - Layer Differences', fontweight='bold')
        ax_diff.set_xlabel('Temperature T')
        ax_diff.set_ylabel(f'Δ{prop} (Layer - Layer 1)')
        
        # Reference data from Layer 1
        layer1_data = None
        
        for layer in [1, 2, 3]:
            if size not in all_data or layer not in all_data[size]:
                continue
                
            data = all_data[size][layer]
            
            # Filter data
            if prop == 'Cv':
                mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.5)
                mask &= (data[prop] > 0) & (data[prop] < 10000)
            else:
                mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.5)
                mask &= data[prop].notna()
            
            if mask.sum() < 10:
                continue
            
            temp = data['Temperature'][mask]
            values = data[prop][mask]
            
            # Store Layer 1 as reference
            if layer == 1:
                layer1_data = {'temp': temp, 'values': values}
            
            # Full range plot
            ax_full.plot(temp, values, 'o-', color=layer_colors[layer-1],
                        alpha=0.7, linewidth=2, markersize=4, label=layer_names[layer-1])
            
            # Critical region (around T=0.8-1.2 based on our analysis)
            crit_mask = (temp >= 0.6) & (temp <= 1.4)
            if crit_mask.sum() > 5:
                ax_crit.plot(temp[crit_mask], values[crit_mask], 'o-',
                           color=layer_colors[layer-1], alpha=0.8, linewidth=2,
                           markersize=6, label=layer_names[layer-1])
            
            # Calculate differences from Layer 1
            if layer > 1 and layer1_data is not None:
                # Interpolate to common temperature grid
                common_temps = np.linspace(max(temp.min(), layer1_data['temp'].min()),
                                         min(temp.max(), layer1_data['temp'].max()), 100)
                
                try:
                    interp1 = np.interp(common_temps, layer1_data['temp'], layer1_data['values'])
                    interp_current = np.interp(common_temps, temp, values)
                    
                    differences = interp_current - interp1
                    
                    ax_diff.plot(common_temps, differences, '-', color=layer_colors[layer-1],
                               linewidth=2, alpha=0.8, label=f'{layer_names[layer-1]} - {layer_names[0]}')
                    
                except Exception as e:
                    print(f"Warning: Could not calculate differences for {prop}, layer {layer}: {e}")
        
        # Format difference plots
        ax_diff.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        
        # Add legends and grids
        for ax in [ax_full, ax_crit, ax_diff]:
            ax.grid(True, alpha=0.3)
            ax.legend()
    
    plt.tight_layout()
    plt.savefig('combined_thermodynamic_properties.png', dpi=300, bbox_inches='tight')
    print(f"Saved: combined_thermodynamic_properties.png")
    plt.show()

def generate_summary_statistics():
    """Generate comprehensive summary statistics for all layers and sizes."""
    all_data = load_all_layer_data()
    lattice_sizes = sorted(all_data.keys())
    
    print(f"\n=== COMPREHENSIVE SUMMARY STATISTICS ===")
    print("="*70)
    
    # Create summary tables
    summary_data = []
    
    for size in lattice_sizes:
        if size not in all_data:
            continue
            
        for layer in [1, 2, 3]:
            if layer not in all_data[size]:
                continue
                
            data = all_data[size][layer]
            
            # Filter reasonable data
            mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.5)
            mask &= (data['Cv'] > 0) & (data['Cv'] < 10000)
            
            if mask.sum() < 10:
                continue
            
            filtered_data = data[mask]
            
            # Calculate key statistics
            try:
                tc_idx = filtered_data['Cv'].idxmax()
                tc = filtered_data.loc[tc_idx, 'Temperature']
                cv_max = filtered_data.loc[tc_idx, 'Cv']
                
                # Energy and magnetization at Tc
                energy_at_tc = filtered_data.loc[tc_idx, 'Energy']
                m2_at_tc = filtered_data.loc[tc_idx, 'M2']
                
                # Low and high temperature limits
                low_temp_data = filtered_data[filtered_data['Temperature'] < 0.8]
                high_temp_data = filtered_data[filtered_data['Temperature'] > 1.8]
                
                cv_low = low_temp_data['Cv'].mean() if len(low_temp_data) > 0 else np.nan
                cv_high = high_temp_data['Cv'].mean() if len(high_temp_data) > 0 else np.nan
                
                m2_low = low_temp_data['M2'].mean() if len(low_temp_data) > 0 else np.nan
                m2_high = high_temp_data['M2'].mean() if len(high_temp_data) > 0 else np.nan
                
                summary_data.append({
                    'Size': size,
                    'Layer': layer,
                    'Tc': tc,
                    'Cv_max': cv_max,
                    'Energy_at_Tc': energy_at_tc,
                    'M2_at_Tc': m2_at_tc,
                    'Cv_low_T': cv_low,
                    'Cv_high_T': cv_high,
                    'M2_low_T': m2_low,
                    'M2_high_T': m2_high,
                    'N_points': mask.sum()
                })
                
            except Exception as e:
                print(f"Warning: Statistics calculation failed for L={size}, Layer={layer}: {e}")
    
    # Convert to DataFrame and display
    if summary_data:
        df = pd.DataFrame(summary_data)
        
        print("\nSummary Table: Critical Temperatures and Peak Values")
        print("-" * 80)
        print(df[['Size', 'Layer', 'Tc', 'Cv_max', 'Energy_at_Tc', 'M2_at_Tc']].to_string(index=False, float_format='%.4f'))
        
        print(f"\nLayer Comparison for each Size:")
        print("-" * 50)
        for size in sorted(df['Size'].unique()):
            size_data = df[df['Size'] == size]
            if len(size_data) == 3:  # All three layers
                tc_values = size_data['Tc'].values
                cv_values = size_data['Cv_max'].values
                
                tc_std = np.std(tc_values)
                cv_std = np.std(cv_values)
                
                print(f"Size {size}x{size}:")
                print(f"  Tc range: {tc_values.min():.4f} - {tc_values.max():.4f} (std: {tc_std:.4f})")
                print(f"  Cv_max range: {cv_values.min():.1f} - {cv_values.max():.1f} (std: {cv_std:.1f})")
        
        # Save summary to file
        df.to_csv('layer_comparison_summary.csv', index=False, float_format='%.6f')
        print(f"\nSaved detailed summary: layer_comparison_summary.csv")
        
        return df
    
    return None

def main():
    """Main analysis function."""
    print("="*70)
    print("COMPREHENSIVE COMBINED LAYER ANALYSIS")
    print("Comparing all three layers across all lattice sizes")
    print("="*70)
    
    # Load and verify data
    all_data = load_all_layer_data()
    
    if not all_data:
        print("Error: No simulation data found!")
        return
    
    lattice_sizes = sorted(all_data.keys())
    print(f"\nAvailable lattice sizes: {lattice_sizes}")
    
    # Count data availability
    data_count = {}
    for size in lattice_sizes:
        data_count[size] = len(all_data[size])
    
    print(f"Data availability: {data_count}")
    
    # Run comprehensive analyses
    print(f"\n{'='*70}")
    print("STARTING COMPREHENSIVE ANALYSIS")
    print("="*70)
    
    try:
        # 1. Critical temperature analysis
        critical_data = analyze_critical_temperatures()
        
        # 2. Thermodynamic properties comparison
        analyze_thermodynamic_comparison()
        
        # 3. Summary statistics
        summary_df = generate_summary_statistics()
        
        print(f"\n{'='*70}")
        print("ANALYSIS COMPLETE")
        print("="*70)
        print("Generated files:")
        print("- combined_layer_critical_analysis.png")
        print("- combined_thermodynamic_properties.png")
        print("- layer_comparison_summary.csv")
        
        return critical_data, summary_df
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return None, None

if __name__ == "__main__":
    main()