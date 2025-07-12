#!/usr/bin/env python3
"""
Automated Multiple Monte Carlo Simulations for Smooth Plots
Runs multiple simulations with different seeds and averages the results
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
from datetime import datetime

def run_single_simulation(nmcs1, nmcs2, seed, output_file):
    """Run a single Monte Carlo simulation"""
    print(f"Running simulation with seed {seed}...")
    
    # Prepare input for the simulation
    input_data = f"{nmcs1} {nmcs2} {seed}\n"
    
    try:
        # Run the simulation
        process = subprocess.Popen(
            ['./icosa8.exe'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate(input=input_data)
        
        if process.returncode != 0:
            print(f"Simulation failed with error: {stderr}")
            return None
        
        # Parse the output
        lines = stdout.strip().split('\n')
        data_lines = []
        
        for line in lines:
            if line.startswith('#'):
                print(f"Header: {line}")
            elif line.strip() and not line.startswith('#'):
                data_lines.append(line)
        
        # Save data to file
        with open(output_file, 'w') as f:
            f.write("# Temperature Magnetization2 Magnetization4 Correlation2 Correlation4 Energy Energy2 HeatCapacity CorrelationLength\n")
            for line in data_lines:
                f.write(line + '\n')
        
        print(f"Simulation completed successfully!")
        print(f"Data saved to {output_file}")
        print(f"Number of temperature points: {len(data_lines)}")
        
        return output_file
        
    except Exception as e:
        print(f"Error running simulation: {e}")
        return None

def load_simulation_data(filename):
    """Load simulation data from file"""
    try:
        data = pd.read_csv(filename, delim_whitespace=True, 
                          names=['temp', 'm2', 'm4', 'g2', 'g4', 'e1', 'cv', 'corr'],
                          comment='#')
        return data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def average_simulation_data(data_files, output_file='averaged_simulation_data.txt'):
    """Average multiple simulation data files"""
    print(f"\nAveraging {len(data_files)} simulation runs...")
    
    all_data = []
    for file in data_files:
        data = load_simulation_data(file)
        if data is not None:
            all_data.append(data)
    
    if not all_data:
        print("No valid data files found!")
        return None
    
    # Calculate average and standard deviation
    avg_data = pd.concat(all_data).groupby(level=0).mean()
    std_data = pd.concat(all_data).groupby(level=0).std()
    
    # Save averaged data
    avg_data.to_csv(output_file, sep='\t', float_format='%.6e')
    std_data.to_csv('error_simulation_data.txt', sep='\t', float_format='%.6e')
    
    print(f"Averaged data saved to {output_file}")
    print(f"Error data saved to error_simulation_data.txt")
    
    return output_file

def create_smooth_plots(data_file):
    """Create smooth plots from averaged data"""
    print("\nCreating smooth plots from averaged data...")
    
    # Import the analyzer
    from analyze_thermodynamics import ThermodynamicAnalyzer
    
    # Load the analyzer with averaged data
    analyzer = ThermodynamicAnalyzer(data_file)
    
    # Create comprehensive plots
    analyzer.plot_thermodynamic_properties('smooth_thermodynamic_properties.png')
    tc_binder, tc_sus, tc_cv, exponents = analyzer.plot_critical_analysis('smooth_critical_analysis.png')
    
    # Generate detailed report
    analyzer.generate_report('smooth_thermodynamic_report.txt')
    
    return tc_binder, tc_sus, tc_cv, exponents

def run_multiple_simulations_automated(nmcs1=1000, nmcs2=100000, num_seeds=10, 
                                     start_seed=12345, output_dir='simulation_runs'):
    """Automate running multiple simulations and averaging results"""
    
    print("="*60)
    print("AUTOMATED MULTIPLE MONTE CARLO SIMULATIONS")
    print("="*60)
    print(f"Parameters:")
    print(f"  nmcs1 (equilibration steps): {nmcs1}")
    print(f"  nmcs2 (measurement steps): {nmcs2}")
    print(f"  Number of seeds: {num_seeds}")
    print(f"  Start seed: {start_seed}")
    print(f"  Output directory: {output_dir}")
    print("="*60)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate seeds
    seeds = [start_seed + i for i in range(num_seeds)]
    
    # Run simulations
    data_files = []
    start_time = time.time()
    
    for i, seed in enumerate(seeds):
        print(f"\n[{i+1}/{num_seeds}] Running simulation with seed {seed}")
        output_file = os.path.join(output_dir, f'simulation_seed_{seed}.txt')
        
        result = run_single_simulation(nmcs1, nmcs2, seed, output_file)
        if result:
            data_files.append(result)
        
        # Estimate remaining time
        if i > 0:
            elapsed = time.time() - start_time
            avg_time = elapsed / i
            remaining = avg_time * (num_seeds - i - 1)
            print(f"Estimated time remaining: {remaining/60:.1f} minutes")
    
    print(f"\nCompleted {len(data_files)} simulations out of {num_seeds}")
    
    if data_files:
        # Average the results
        avg_file = average_simulation_data(data_files, 'smooth_averaged_data.txt')
        
        # Create smooth plots
        print("\nCreating smooth plots...")
        tc_binder, tc_sus, tc_cv, exponents = create_smooth_plots(avg_file)
        
        # Print summary
        print("\n" + "="*60)
        print("SIMULATION SUMMARY")
        print("="*60)
        print(f"Total simulations completed: {len(data_files)}")
        print(f"Total measurement steps per temperature: {nmcs2 * len(data_files):,}")
        print(f"Critical temperatures:")
        print(f"  Binder cumulant method: {tc_binder:.6f}")
        print(f"  Susceptibility peak:    {tc_sus:.6f}")
        print(f"  Heat capacity peak:     {tc_cv:.6f}")
        
        if exponents:
            print(f"Critical exponents:")
            print(f"  β = {exponents['beta']:.4f}")
            print(f"  γ = {exponents['gamma']:.4f}")
        
        print(f"\nFiles generated:")
        print(f"  - smooth_averaged_data.txt (averaged data)")
        print(f"  - error_simulation_data.txt (error bars)")
        print(f"  - smooth_thermodynamic_properties.png")
        print(f"  - smooth_critical_analysis.png")
        print(f"  - smooth_thermodynamic_report.txt")
        print(f"  - Individual runs in {output_dir}/")
        
        total_time = time.time() - start_time
        print(f"\nTotal computation time: {total_time/60:.1f} minutes")
        print("="*60)
        
        return avg_file
    else:
        print("No simulations completed successfully!")
        return None

def main():
    """Main function with user-configurable parameters"""
    print("Automated Multiple Monte Carlo Simulations for Smooth Plots")
    print("="*60)
    
    # Configuration - modify these parameters as needed
    config = {
        'nmcs1': 1000,        # Equilibration steps
        'nmcs2': 100000,      # Measurement steps (increase for smoother plots)
        'num_seeds': 10,       # Number of different random seeds
        'start_seed': 12345,   # Starting seed number
        'output_dir': 'simulation_runs'  # Directory for individual runs
    }
    
    print("Configuration:")
    for key, value in config.items():
        print(f"  {key}: {value}")
    
    print(f"\nThis will run {config['num_seeds']} simulations with {config['nmcs2']:,} measurement steps each.")
    print(f"Total measurement steps per temperature: {config['nmcs2'] * config['num_seeds']:,}")
    
    # Estimate computation time (rough estimate)
    estimated_time = config['num_seeds'] * config['nmcs2'] / 1000000  # Rough estimate in minutes
    print(f"Estimated computation time: {estimated_time:.1f} minutes")
    
    # Ask for confirmation
    response = input("\nProceed with these parameters? (y/n): ").lower().strip()
    if response != 'y':
        print("Simulation cancelled.")
        return
    
    # Run the simulations
    result = run_multiple_simulations_automated(**config)
    
    if result:
        print(f"\n🎉 Analysis complete! Smooth plots and data saved.")
        print(f"Ready for Q1 journal submission with high-quality results!")
    else:
        print(f"\n❌ Analysis failed. Please check the error messages above.")

if __name__ == "__main__":
    main() 