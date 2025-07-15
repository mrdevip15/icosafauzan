#!/usr/bin/env python3
"""
Python runner for the enhanced cubic3lay Monte Carlo simulation.
Features:
- Enhanced 8-state Potts theory with optimized temperature range
- Icosahedral spin model (13 states) with void parameter
- Hybrid Monte Carlo algorithm (Metropolis + Wolff clusters)
- First-order transition detection
- Timestamped hierarchical directory structure
"""

import subprocess
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def run_simulation(nx, ny, nmcs1, nmcs2, seed):
    """Run the enhanced Monte Carlo simulation."""
    try:
        # Prepare input
        input_data = f"{nx} {ny} {nmcs1} {nmcs2} {seed}\n"
        
        print(f"Running simulation: {nx}x{ny}, eq={nmcs1}, meas={nmcs2}, seed={seed}")
        
        # Run simulation
        result = subprocess.run(
            ['./cubic3lay_enhanced'], 
            input=input_data, 
            text=True, 
            capture_output=True, 
            timeout=300
        )
        
        if result.returncode == 0:
            # Extract output directory from stdout
            output_lines = result.stdout.strip().split('\n')
            for line in output_lines:
                if "Layer data saved to" in line:
                    output_dir = line.split("to ")[-1].rstrip('/')
                    print(f"✓ Success: Data saved to {output_dir}")
                    return output_dir
            
            print("✓ Simulation completed but no output directory found")
            return None
        else:
            print(f"✗ Simulation failed: {result.stderr}")
            return None
            
    except subprocess.TimeoutExpired:
        print("✗ Simulation timed out (>5 minutes)")
        return None
    except FileNotFoundError:
        print("✗ Executable './cubic3lay_enhanced' not found")
        print("   Compile first: gcc -o cubic3lay_enhanced cubic3lay.c rn32.c -lm")
        return None
    except Exception as e:
        print(f"✗ Error: {e}")
        return None

def analyze_simulation_data(output_dir):
    """Analyze the simulation results."""
    if not output_dir or not os.path.exists(output_dir):
        print("No valid output directory to analyze")
        return None
    
    print(f"\n=== ANALYZING RESULTS FROM {output_dir} ===")
    
    results = {}
    
    # Analyze each layer
    for layer in [1, 2, 3]:
        layer_file = f"{output_dir}/layer_{layer}.txt"
        
        if os.path.exists(layer_file):
            print(f"\nLayer {layer} analysis:")
            
            # Read data (skip comment lines)
            try:
                data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                 names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                
                # Basic statistics
                n_points = len(data)
                temp_range = (data['Temperature'].min(), data['Temperature'].max())
                
                print(f"  Data points: {n_points}")
                print(f"  Temperature range: {temp_range[0]:.3f} - {temp_range[1]:.3f}")
                
                # Critical region analysis (around theoretical Tc ≈ 0.751)
                critical_mask = (data['Temperature'] >= 0.6) & (data['Temperature'] <= 0.9)
                critical_data = data[critical_mask]
                
                if len(critical_data) > 0:
                    max_cv_idx = critical_data['Cv'].idxmax()
                    tc_apparent = critical_data.loc[max_cv_idx, 'Temperature']
                    cv_max = critical_data.loc[max_cv_idx, 'Cv']
                    
                    print(f"  Critical region points: {len(critical_data)}")
                    print(f"  Apparent Tc (Cv peak): {tc_apparent:.4f}")
                    print(f"  Peak Cv value: {cv_max:.2f}")
                    
                    # Compare with theory
                    theoretical_tc = 1.0 / np.log(1 + np.sqrt(8))  # 8-state Potts
                    print(f"  Theoretical Tc (8-state): {theoretical_tc:.4f}")
                    print(f"  Difference: {tc_apparent - theoretical_tc:+.4f}")
                else:
                    tc_apparent = np.nan
                    cv_max = np.nan
                    print("  No critical region data found")
                
                # Energy analysis
                energy_mean = data['Energy'].mean()
                energy_std = data['Energy'].std()
                print(f"  Average energy: {energy_mean:.4f} ± {energy_std:.4f}")
                
                # Store results
                results[layer] = {
                    'n_points': n_points,
                    'temp_range': temp_range,
                    'tc_apparent': tc_apparent,
                    'cv_max': cv_max,
                    'energy_mean': energy_mean,
                    'energy_std': energy_std,
                    'data': data
                }
                
            except Exception as e:
                print(f"  Error reading data: {e}")
                results[layer] = None
        else:
            print(f"Layer {layer} file not found: {layer_file}")
            results[layer] = None
    
    return results

def plot_results(results, output_dir):
    """Create plots of the simulation results."""
    if not results or not any(results.values()):
        print("No data to plot")
        return
    
    print(f"\n=== CREATING PLOTS ===")
    
    # Create plot directory
    plot_dir = f"{output_dir}/plots"
    os.makedirs(plot_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Enhanced Quasi-3D Monte Carlo Simulation Results\n'
                 'Icosahedral Spins (13-state) with Hybrid Algorithm', fontsize=16, fontweight='bold')
    
    layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']
    layer_names = ['Layer 1', 'Layer 2', 'Layer 3']
    
    # Theoretical critical temperature
    theoretical_tc = 1.0 / np.log(1 + np.sqrt(8))
    
    for layer in [1, 2, 3]:
        if results[layer] is None:
            continue
            
        data = results[layer]['data']
        color = layer_colors[layer-1]
        name = layer_names[layer-1]
        
        # Plot 1: Specific Heat vs Temperature
        ax1 = axes[0, 0]
        ax1.plot(data['Temperature'], data['Cv'], 'o-', color=color, alpha=0.7, 
                markersize=3, linewidth=1.5, label=name)
        ax1.axvline(theoretical_tc, color='black', linestyle='--', alpha=0.5, 
                   label='Theoretical Tc' if layer == 1 else '')
        ax1.set_xlabel('Temperature T')
        ax1.set_ylabel('Specific Heat $C_v$')
        ax1.set_title('Specific Heat vs Temperature')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Magnetization vs Temperature
        ax2 = axes[0, 1]
        ax2.plot(data['Temperature'], data['M2'], 'o-', color=color, alpha=0.7,
                markersize=3, linewidth=1.5, label=name)
        ax2.axvline(theoretical_tc, color='black', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Temperature T')
        ax2.set_ylabel('Magnetization $\\langle M^2 \\rangle$')
        ax2.set_title('Magnetization vs Temperature')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Energy vs Temperature
        ax3 = axes[0, 2]
        ax3.plot(data['Temperature'], data['Energy'], 'o-', color=color, alpha=0.7,
                markersize=3, linewidth=1.5, label=name)
        ax3.axvline(theoretical_tc, color='black', linestyle='--', alpha=0.5)
        ax3.set_xlabel('Temperature T')
        ax3.set_ylabel('Energy per site')
        ax3.set_title('Energy vs Temperature')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Critical region zoom (Specific Heat)
        ax4 = axes[1, 0]
        critical_mask = (data['Temperature'] >= 0.6) & (data['Temperature'] <= 0.9)
        if critical_mask.sum() > 0:
            ax4.plot(data['Temperature'][critical_mask], data['Cv'][critical_mask], 
                    'o-', color=color, alpha=0.8, markersize=4, linewidth=2, label=name)
        ax4.axvline(theoretical_tc, color='black', linestyle='--', alpha=0.5)
        ax4.set_xlabel('Temperature T')
        ax4.set_ylabel('Specific Heat $C_v$')
        ax4.set_title('Critical Region (0.6 - 0.9)')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        # Plot 5: Susceptibility
        ax5 = axes[1, 1]
        ax5.plot(data['Temperature'], data['G2'], 'o-', color=color, alpha=0.7,
                markersize=3, linewidth=1.5, label=name)
        ax5.axvline(theoretical_tc, color='black', linestyle='--', alpha=0.5)
        ax5.set_xlabel('Temperature T')
        ax5.set_ylabel('Susceptibility $\\chi$')
        ax5.set_title('Susceptibility vs Temperature')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
        
        # Plot 6: Correlation Length
        ax6 = axes[1, 2]
        # Only plot finite correlation lengths
        finite_corr = data['Corr'] > 0
        if finite_corr.sum() > 0:
            ax6.plot(data['Temperature'][finite_corr], data['Corr'][finite_corr], 
                    'o-', color=color, alpha=0.7, markersize=3, linewidth=1.5, label=name)
        ax6.axvline(theoretical_tc, color='black', linestyle='--', alpha=0.5)
        ax6.set_xlabel('Temperature T')
        ax6.set_ylabel('Correlation Length $\\xi$')
        ax6.set_title('Correlation Length vs Temperature')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plots
    plot_file = f"{plot_dir}/thermodynamic_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plots saved to: {plot_file}")
    
    plt.show()

def print_summary(results, output_dir):
    """Print a summary of the simulation results."""
    if not results:
        return
    
    print(f"\n{'='*60}")
    print("ENHANCED SIMULATION SUMMARY")
    print(f"{'='*60}")
    
    # Extract metadata from first layer file
    layer1_file = f"{output_dir}/layer_1.txt"
    if os.path.exists(layer1_file):
        print("Model Configuration:")
        with open(layer1_file, 'r') as f:
            for line in f:
                if line.startswith('#') and any(key in line for key in 
                    ['Spin model', 'Algorithm', 'Theoretical Tc', 'Temperature range', 'System size', 'Void parameter']):
                    print(f"  {line[1:].strip()}")
                elif not line.startswith('#'):
                    break
    
    print(f"\nResults Summary:")
    theoretical_tc = 1.0 / np.log(1 + np.sqrt(8))
    
    for layer in [1, 2, 3]:
        if results[layer] is not None:
            r = results[layer]
            print(f"\nLayer {layer}:")
            print(f"  Data points: {r['n_points']}")
            print(f"  Temperature range: {r['temp_range'][0]:.3f} - {r['temp_range'][1]:.3f}")
            if not np.isnan(r['tc_apparent']):
                print(f"  Apparent Tc: {r['tc_apparent']:.4f}")
                print(f"  Peak Cv: {r['cv_max']:.2f}")
                print(f"  Tc deviation from theory: {r['tc_apparent'] - theoretical_tc:+.4f}")
            print(f"  Average energy: {r['energy_mean']:.4f} ± {r['energy_std']:.4f}")
    
    # Layer consistency check
    tc_values = [results[layer]['tc_apparent'] for layer in [1, 2, 3] 
                if results[layer] is not None and not np.isnan(results[layer]['tc_apparent'])]
    
    if len(tc_values) > 1:
        tc_std = np.std(tc_values)
        tc_mean = np.mean(tc_values)
        print(f"\nLayer Consistency:")
        print(f"  Mean Tc across layers: {tc_mean:.4f}")
        print(f"  Standard deviation: {tc_std:.4f}")
        print(f"  Consistency: {'Excellent' if tc_std < 0.01 else 'Good' if tc_std < 0.05 else 'Poor'}")
    
    print(f"\nKey Enhancements:")
    print(f"  ✓ Optimized temperature range for 8-state Potts critical region")
    print(f"  ✓ Icosahedral spin model (13 states vs 8 cube states)")
    print(f"  ✓ Hybrid Monte Carlo algorithm (Metropolis + Wolff)")
    print(f"  ✓ First-order transition detection capabilities")
    print(f"  ✓ Timestamped hierarchical output structure")

def main():
    """Main function to run enhanced simulations."""
    print("="*60)
    print("ENHANCED QUASI-3D MONTE CARLO SIMULATION RUNNER")
    print("="*60)
    print("Features:")
    print("- 8-state Potts theory with optimized temperature range")
    print("- Icosahedral spin model (13 states with void parameter)")
    print("- Hybrid Monte Carlo algorithm (70% Metropolis + 30% Wolff)")
    print("- First-order transition detection and histogram analysis")
    print("- Timestamped hierarchical directory structure")
    print("-"*60)
    
    # Simulation parameters
    lattice_sizes = [8, 16]  # Start with smaller sizes
    nmcs1 = 1000   # Equilibration steps  
    nmcs2 = 2000   # Measurement steps
    base_seed = 12345
    
    print(f"Simulation Parameters:")
    print(f"- Lattice sizes: {lattice_sizes}")
    print(f"- Equilibration steps: {nmcs1}")
    print(f"- Measurement steps: {nmcs2}")
    print(f"- Theoretical Tc (8-state Potts): {1.0/np.log(1+np.sqrt(8)):.4f}")
    print(f"- Temperature range: 0.1 - 1.5 (fine resolution near Tc)")
    print()
    
    # Check executable
    if not os.path.exists('./cubic3lay_enhanced'):
        print("✗ Executable './cubic3lay_enhanced' not found")
        print("  Compile first: gcc -o cubic3lay_enhanced cubic3lay.c rn32.c -lm")
        return False
    
    # Run simulations
    simulation_results = []
    
    for i, nx in enumerate(lattice_sizes):
        ny = nx  # Square lattices
        seed = base_seed + i
        
        print(f"{'='*40}")
        print(f"SIMULATION {i+1}/{len(lattice_sizes)}: {nx}x{nx}")
        print(f"{'='*40}")
        
        output_dir = run_simulation(nx, ny, nmcs1, nmcs2, seed)
        
        if output_dir:
            # Analyze results
            results = analyze_simulation_data(output_dir)
            
            if results:
                # Create plots
                plot_results(results, output_dir)
                
                # Print summary
                print_summary(results, output_dir)
                
                simulation_results.append({
                    'size': nx,
                    'output_dir': output_dir,
                    'results': results
                })
        
        print()  # Spacing between simulations
    
    # Final summary
    print(f"{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    
    if simulation_results:
        print(f"✓ Completed {len(simulation_results)}/{len(lattice_sizes)} simulations")
        
        print(f"\nGenerated directories:")
        for sim in simulation_results:
            print(f"  L={sim['size']}: {sim['output_dir']}")
        
        print(f"\nDirectory structure example:")
        if simulation_results:
            example_dir = simulation_results[0]['output_dir']
            parts = example_dir.split('/')
            if len(parts) >= 3:
                print(f"  simulation_runs/")
                print(f"  ├── {parts[1]}/  (timestamp_model_algorithm)")
                print(f"  │   ├── {parts[2]}/  (lattice_size)")
                print(f"  │   │   ├── layer_1.txt")
                print(f"  │   │   ├── layer_2.txt")
                print(f"  │   │   └── layer_3.txt")
                print(f"  │   └── plots/")
    else:
        print(f"✗ No simulations completed successfully")
    
    return len(simulation_results) > 0

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)