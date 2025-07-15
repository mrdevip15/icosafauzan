#!/usr/bin/env python3
"""
Comprehensive thermodynamic analysis for layer-separated quasi-3D simulation data.
Analyzes magnetization, energy, heat capacity, and other properties for each layer
and compares across different lattice sizes.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from pathlib import Path
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

class LayerThermodynamicAnalyzer:
    def __init__(self, data_dir="simulation_runs"):
        """Initialize analyzer with simulation data directory."""
        self.data_dir = data_dir
        self.lattice_sizes = []
        self.layer_data = {}
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
        self.layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']  # Red, Blue, Green for layers 1, 2, 3
        
    def load_data(self):
        """Load all simulation data from directory structure."""
        print("Loading simulation data...")
        
        # Find all lattice size directories
        for size_dir in glob.glob(f"{self.data_dir}/nx*_ny*"):
            size_name = os.path.basename(size_dir)
            # Extract nx, ny from directory name
            parts = size_name.split('_')
            nx = int(parts[0][2:])  # Remove 'nx' prefix
            ny = int(parts[1][2:])  # Remove 'ny' prefix
            
            if nx == ny:  # Only consider square lattices
                self.lattice_sizes.append(nx)
                self.layer_data[nx] = {}
                
                # Load data for each layer
                for layer in [1, 2, 3]:
                    layer_file = f"{size_dir}/layer_{layer}.txt"
                    if os.path.exists(layer_file):
                        try:
                            data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                             names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                            self.layer_data[nx][layer] = data
                            print(f"Loaded {size_name}/layer_{layer}.txt: {len(data)} temperature points")
                        except Exception as e:
                            print(f"Error loading {layer_file}: {e}")
        
        self.lattice_sizes.sort()
        print(f"Found lattice sizes: {self.lattice_sizes}")
        
    def critical_temp_analysis(self, size):
        """Estimate critical temperature using Binder cumulant method."""
        if size not in self.layer_data:
            return None, None
            
        results = {}
        for layer in [1, 2, 3]:
            if layer not in self.layer_data[size]:
                continue
                
            data = self.layer_data[size][layer]
            
            # Calculate Binder cumulant U = 1 - <M^4>/(3<M^2>^2)
            valid_mask = (data['M2'] > 1e-10) & (data['M4'] > 1e-10)
            if not valid_mask.any():
                continue
                
            valid_data = data[valid_mask]
            binder = 1 - valid_data['M4'] / (3 * valid_data['M2']**2)
            
            # Find temperature where Binder cumulant ≈ 0.61 (2D Ising universal value)
            target_binder = 0.61
            if len(binder) > 10:
                # Interpolate to find Tc
                temp_range = valid_data['Temperature']
                binder_interp = interp1d(temp_range, binder, kind='linear', fill_value='extrapolate')
                
                # Find crossing point
                temp_fine = np.linspace(temp_range.min(), temp_range.max(), 1000)
                binder_fine = binder_interp(temp_fine)
                crossing_idx = np.argmin(np.abs(binder_fine - target_binder))
                
                tc_estimate = temp_fine[crossing_idx]
                results[layer] = tc_estimate
                
        return results
    
    def plot_thermodynamic_properties(self):
        """Plot all thermodynamic properties for each layer and lattice size."""
        properties = [
            ('M2', 'Magnetization $\\langle M^2 \\rangle$'),
            ('Energy', 'Energy per site'),
            ('Cv', 'Heat Capacity $C_v$'),
            ('Corr', 'Correlation Length $\\xi$')
        ]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        for i, (prop, label) in enumerate(properties):
            ax = axes[i]
            
            for j, size in enumerate(self.lattice_sizes):
                if size not in self.layer_data:
                    continue
                    
                for layer in [1, 2, 3]:
                    if layer not in self.layer_data[size]:
                        continue
                        
                    data = self.layer_data[size][layer]
                    
                    # Filter reasonable temperature range
                    mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.0)
                    temp = data['Temperature'][mask]
                    values = data[prop][mask]
                    
                    # Special handling for correlation length (remove infinite/zero values)
                    if prop == 'Corr':
                        corr_mask = (values > 0) & (values < 1000) & np.isfinite(values)
                        temp = temp[corr_mask]
                        values = values[corr_mask]
                    
                    if len(temp) > 0:
                        linestyle = '-' if layer == 1 else '--' if layer == 2 else ':'
                        alpha = 0.7 + 0.1 * j  # Vary transparency by size
                        ax.plot(temp, values, 
                               color=self.layer_colors[layer-1], 
                               linestyle=linestyle,
                               linewidth=2, 
                               alpha=alpha,
                               label=f'L={size}, Layer {layer}' if j == 0 else '')
            
            ax.set_xlabel('Temperature $T$')
            ax.set_ylabel(label)
            ax.grid(True, alpha=0.3)
            
            if i == 0:  # Only show legend for first subplot
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
        plt.tight_layout()
        plt.savefig('thermodynamic_properties_all_layers.png', dpi=300, bbox_inches='tight')
        print("Saved: thermodynamic_properties_all_layers.png")
        plt.show()
    
    def plot_lattice_size_comparison(self):
        """Compare properties across different lattice sizes for each layer."""
        fig, axes = plt.subplots(3, 2, figsize=(16, 18))
        
        properties = [
            ('M2', 'Magnetization $\\langle M^2 \\rangle$'),
            ('Cv', 'Heat Capacity $C_v$')
        ]
        
        for layer in [1, 2, 3]:
            for j, (prop, label) in enumerate(properties):
                ax = axes[layer-1, j]
                
                for i, size in enumerate(self.lattice_sizes):
                    if size not in self.layer_data or layer not in self.layer_data[size]:
                        continue
                        
                    data = self.layer_data[size][layer]
                    
                    # Focus on critical region
                    mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 1.5)
                    temp = data['Temperature'][mask]
                    values = data[prop][mask]
                    
                    if len(temp) > 0:
                        ax.plot(temp, values, 
                               color=self.colors[i % len(self.colors)],
                               linewidth=2.5, 
                               marker='o', markersize=3,
                               label=f'L = {size}')
                
                ax.set_xlabel('Temperature $T$')
                ax.set_ylabel(label)
                ax.set_title(f'Layer {layer}')
                ax.grid(True, alpha=0.3)
                ax.legend()
                
        plt.tight_layout()
        plt.savefig('lattice_size_comparison_by_layer.png', dpi=300, bbox_inches='tight')
        print("Saved: lattice_size_comparison_by_layer.png")
        plt.show()
    
    def plot_layer_comparison(self):
        """Compare properties between layers for each lattice size."""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        properties = [
            ('M2', 'Magnetization $\\langle M^2 \\rangle$'),
            ('Energy', 'Energy per site'),
            ('Cv', 'Heat Capacity $C_v$'),
            ('G2', 'Correlation $G^2$')
        ]
        
        for i, (prop, label) in enumerate(properties):
            ax = axes[i//2, i%2]
            
            for size in self.lattice_sizes[-2:]:  # Show only largest two sizes for clarity
                if size not in self.layer_data:
                    continue
                    
                for layer in [1, 2, 3]:
                    if layer not in self.layer_data[size]:
                        continue
                        
                    data = self.layer_data[size][layer]
                    
                    # Focus on critical region
                    mask = (data['Temperature'] >= 0.8) & (data['Temperature'] <= 1.8)
                    temp = data['Temperature'][mask]
                    values = data[prop][mask]
                    
                    if len(temp) > 0:
                        linestyle = '-' if size == self.lattice_sizes[-1] else '--'
                        ax.plot(temp, values, 
                               color=self.layer_colors[layer-1],
                               linestyle=linestyle,
                               linewidth=2.5, 
                               alpha=0.8,
                               label=f'L={size}, Layer {layer}')
            
            ax.set_xlabel('Temperature $T$')
            ax.set_ylabel(label)
            ax.grid(True, alpha=0.3)
            ax.legend()
                
        plt.tight_layout()
        plt.savefig('layer_comparison_properties.png', dpi=300, bbox_inches='tight')
        print("Saved: layer_comparison_properties.png")
        plt.show()
    
    def plot_critical_temperature_analysis(self):
        """Analyze and plot critical temperature estimates."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        tc_data = {}
        for size in self.lattice_sizes:
            tc_results = self.critical_temp_analysis(size)
            if tc_results:
                tc_data[size] = tc_results
        
        # Plot 1: Critical temperature vs lattice size for each layer
        for layer in [1, 2, 3]:
            sizes = []
            tcs = []
            for size, tc_dict in tc_data.items():
                if layer in tc_dict:
                    sizes.append(size)
                    tcs.append(tc_dict[layer])
            
            if len(sizes) > 1:
                ax1.plot(sizes, tcs, 'o-', 
                        color=self.layer_colors[layer-1],
                        linewidth=2, markersize=8,
                        label=f'Layer {layer}')
        
        ax1.set_xlabel('Lattice Size $L$')
        ax1.set_ylabel('Critical Temperature $T_c$')
        ax1.set_title('Critical Temperature vs Lattice Size')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Plot 2: Finite-size scaling
        if len(self.lattice_sizes) > 2:
            for layer in [1, 2, 3]:
                sizes = []
                tcs = []
                for size, tc_dict in tc_data.items():
                    if layer in tc_dict:
                        sizes.append(1.0/size)  # 1/L for finite-size scaling
                        tcs.append(tc_dict[layer])
                
                if len(sizes) > 1:
                    ax2.plot(sizes, tcs, 'o-', 
                            color=self.layer_colors[layer-1],
                            linewidth=2, markersize=8,
                            label=f'Layer {layer}')
        
        ax2.set_xlabel('$1/L$')
        ax2.set_ylabel('Critical Temperature $T_c$')
        ax2.set_title('Finite-Size Scaling of $T_c$')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig('critical_temperature_analysis.png', dpi=300, bbox_inches='tight')
        print("Saved: critical_temperature_analysis.png")
        plt.show()
        
        # Print critical temperature estimates
        print("\n=== Critical Temperature Estimates ===")
        for size, tc_dict in tc_data.items():
            print(f"Lattice size {size}x{size}:")
            for layer, tc in tc_dict.items():
                print(f"  Layer {layer}: Tc ≈ {tc:.3f}")
    
    def generate_summary_statistics(self):
        """Generate summary statistics for all layers and sizes."""
        print("\n=== Summary Statistics ===")
        
        for size in self.lattice_sizes:
            if size not in self.layer_data:
                continue
                
            print(f"\nLattice size: {size}x{size}")
            print("-" * 40)
            
            for layer in [1, 2, 3]:
                if layer not in self.layer_data[size]:
                    continue
                    
                data = self.layer_data[size][layer]
                
                # Critical region analysis (T = 1.0 to 1.5)
                critical_mask = (data['Temperature'] >= 1.0) & (data['Temperature'] <= 1.5)
                critical_data = data[critical_mask]
                
                if len(critical_data) > 0:
                    max_cv_idx = critical_data['Cv'].idxmax()
                    max_cv_temp = critical_data.loc[max_cv_idx, 'Temperature']
                    max_cv_value = critical_data.loc[max_cv_idx, 'Cv']
                    
                    print(f"  Layer {layer}:")
                    print(f"    Max Cv: {max_cv_value:.3f} at T = {max_cv_temp:.3f}")
                    print(f"    Energy range: {critical_data['Energy'].min():.3f} to {critical_data['Energy'].max():.3f}")
                    print(f"    M² range: {critical_data['M2'].min():.3f} to {critical_data['M2'].max():.3f}")
    
    def run_full_analysis(self):
        """Run complete thermodynamic analysis."""
        print("=== Layer-Separated Quasi-3D Thermodynamic Analysis ===")
        print("Monte Carlo simulation of 8-state cube spins on 3-layer lattice")
        print("=" * 60)
        
        self.load_data()
        
        if not self.layer_data:
            print("No simulation data found!")
            return
        
        print(f"\nAnalyzing {len(self.lattice_sizes)} lattice sizes: {self.lattice_sizes}")
        print("Generating comprehensive plots...")
        
        # Generate all plots
        self.plot_thermodynamic_properties()
        self.plot_lattice_size_comparison() 
        self.plot_layer_comparison()
        self.plot_critical_temperature_analysis()
        
        # Generate statistics
        self.generate_summary_statistics()
        
        print("\n=== Analysis Complete ===")
        print("Generated plots:")
        print("- thermodynamic_properties_all_layers.png")
        print("- lattice_size_comparison_by_layer.png") 
        print("- layer_comparison_properties.png")
        print("- critical_temperature_analysis.png")

def main():
    """Main analysis function."""
    analyzer = LayerThermodynamicAnalyzer()
    analyzer.run_full_analysis()

if __name__ == "__main__":
    main()