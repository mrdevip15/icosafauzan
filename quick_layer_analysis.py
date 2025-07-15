#!/usr/bin/env python3
"""
Quick thermodynamic analysis for layer-separated data with focused plotting.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob

def load_layer_data():
    """Load simulation data from existing files."""
    data_dir = "simulation_runs"
    layer_data = {}
    
    for size_dir in glob.glob(f"{data_dir}/nx*_ny*"):
        size_name = os.path.basename(size_dir)
        parts = size_name.split('_')
        nx = int(parts[0][2:])
        ny = int(parts[1][2:])
        
        if nx == ny:  # Square lattices only
            layer_data[nx] = {}
            
            for layer in [1, 2, 3]:
                layer_file = f"{size_dir}/layer_{layer}.txt"
                if os.path.exists(layer_file):
                    try:
                        data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                         names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                        layer_data[nx][layer] = data
                        print(f"Loaded {size_name}/layer_{layer}.txt: {len(data)} points")
                    except Exception as e:
                        print(f"Error loading {layer_file}: {e}")
    
    return layer_data

def plot_key_properties():
    """Plot key thermodynamic properties."""
    layer_data = load_layer_data()
    lattice_sizes = sorted(layer_data.keys())
    
    print(f"Found lattice sizes: {lattice_sizes}")
    
    # Layer colors
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']  # Red, Blue, Green
    
    # Create subplots for key properties
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    properties = [
        ('M2', 'Magnetization $\\langle M^2 \\rangle$'),
        ('Energy', 'Energy per site'),
        ('Cv', 'Heat Capacity $C_v$'),
        ('G2', 'Correlation Function $G^2$')
    ]
    
    for i, (prop, label) in enumerate(properties):
        ax = axes[i//2, i%2]
        
        for size in lattice_sizes:
            for layer in [1, 2, 3]:
                if layer not in layer_data[size]:
                    continue
                    
                data = layer_data[size][layer]
                
                # Focus on critical region
                mask = (data['Temperature'] >= 0.7) & (data['Temperature'] <= 1.8)
                temp = data['Temperature'][mask]
                values = data[prop][mask]
                
                if len(temp) > 5:
                    linestyle = '-' if size == lattice_sizes[-1] else '--'
                    alpha = 0.8 if size == lattice_sizes[-1] else 0.6
                    
                    ax.plot(temp, values,
                           color=layer_colors[layer-1],
                           linestyle=linestyle,
                           linewidth=2,
                           alpha=alpha,
                           label=f'L={size}, Layer {layer}' if size == lattice_sizes[-1] else '')
        
        ax.set_xlabel('Temperature $T$')
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend()
    
    plt.tight_layout()
    plt.savefig('layer_thermodynamics_analysis.png', dpi=300, bbox_inches='tight')
    print("Saved: layer_thermodynamics_analysis.png")
    plt.show()

def analyze_critical_region():
    """Analyze critical region and find peak heat capacity."""
    layer_data = load_layer_data()
    
    print("\n=== Critical Region Analysis ===")
    
    for size in sorted(layer_data.keys()):
        print(f"\nLattice size: {size}x{size}")
        print("-" * 30)
        
        for layer in [1, 2, 3]:
            if layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            
            # Find peak heat capacity in critical region
            critical_mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 1.5)
            critical_data = data[critical_mask]
            
            if len(critical_data) > 0:
                max_cv_idx = critical_data['Cv'].idxmax()
                max_cv_temp = critical_data.loc[max_cv_idx, 'Temperature']
                max_cv_value = critical_data.loc[max_cv_idx, 'Cv']
                
                # Get corresponding values at peak
                peak_m2 = critical_data.loc[max_cv_idx, 'M2']
                peak_energy = critical_data.loc[max_cv_idx, 'Energy']
                
                print(f"  Layer {layer}:")
                print(f"    Peak Cv: {max_cv_value:.3f} at T = {max_cv_temp:.3f}")
                print(f"    M² at peak: {peak_m2:.3f}")
                print(f"    Energy at peak: {peak_energy:.3f}")

def compare_layers():
    """Compare properties between layers at same lattice size."""
    layer_data = load_layer_data()
    lattice_sizes = sorted(layer_data.keys())
    
    if not lattice_sizes:
        print("No data found!")
        return
    
    # Use largest lattice size for comparison
    size = lattice_sizes[-1]
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    
    properties = [('M2', 'Magnetization'), ('Cv', 'Heat Capacity'), ('Energy', 'Energy')]
    
    for i, (prop, title) in enumerate(properties):
        ax = axes[i]
        
        for layer in [1, 2, 3]:
            if layer not in layer_data[size]:
                continue
                
            data = layer_data[size][layer]
            mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 1.6)
            temp = data['Temperature'][mask]
            values = data[prop][mask]
            
            ax.plot(temp, values, 
                   color=layer_colors[layer-1],
                   linewidth=2.5,
                   label=f'Layer {layer}')
        
        ax.set_xlabel('Temperature')
        ax.set_ylabel(title)
        ax.set_title(f'{title} - L={size}x{size}')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig('layer_comparison_L{}.png'.format(size), dpi=300, bbox_inches='tight')
    print(f"Saved: layer_comparison_L{size}.png")
    plt.show()

def main():
    """Main analysis function."""
    print("=== Quick Layer-Separated Thermodynamic Analysis ===")
    
    try:
        plot_key_properties()
        analyze_critical_region()
        compare_layers()
        
        print("\n=== Analysis Complete ===")
        print("Check the generated PNG files for detailed plots.")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()