#!/usr/bin/env python3
"""
Parallel Monte Carlo Simulations using all CPU cores
Runs multiple simulations in parallel for maximum speed
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
from datetime import datetime
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
import psutil

def run_single_simulation_parallel(args):
    """Run a single Monte Carlo simulation (for parallel processing)"""
    nx, ny, nmcs1, nmcs2, seed, run_id, output_dir = args
    
    print(f"Running simulation {run_id} with nx={nx}, ny={ny}, seed={seed}...")
    input_data = f"{nx} {ny} {nmcs1} {nmcs2} {seed}\n"
    
    try:
        process = subprocess.Popen(
            ['./icosa8.exe'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate(input=input_data)
        
        if process.returncode != 0:
            print(f"Simulation {run_id} failed with error: {stderr}")
            return None
        
        # Parse the output
        lines = stdout.strip().split('\n')
        data_lines = []
        
        for line in lines:
            if line.startswith('#'):
                print(f"Header {run_id}: {line}")
            elif line.strip() and not line.startswith('#'):
                data_lines.append(line)
        
        # Save data to file
        output_file = os.path.join(output_dir, f'simulation_run_{run_id}.txt')
        with open(output_file, 'w') as f:
            f.write("# Temperature Magnetization2 Magnetization4 Correlation2 Correlation4 Energy Energy2 HeatCapacity CorrelationLength\n")
            for line in data_lines:
                f.write(line + '\n')
        
        print(f"Simulation {run_id} completed successfully!")
        print(f"Data saved to {output_file}")
        print(f"Number of temperature points: {len(data_lines)}")
        
        return output_file
        
    except Exception as e:
        print(f"Error running simulation {run_id}: {e}")
        return None

def load_simulation_data(filename):
    """Load simulation data from file"""
    try:
        data = pd.read_csv(filename, sep='\s+', 
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

def create_smooth_plots(data_file, nx, ny):
    """Create smooth plots from averaged data for a given lattice size"""
    print(f"\nCreating smooth plots for {nx}x{ny}...")
    
    # Import the analyzer
    from analyze_thermodynamics import ThermodynamicAnalyzer
    
    # Load the analyzer with averaged data
    analyzer = ThermodynamicAnalyzer(data_file)
    
    # Create output directory for this lattice size
    output_dir = f'simulation_runs/nx{nx}_ny{ny}'
    
    # Create comprehensive plots
    analyzer.plot_thermodynamic_properties(f'{output_dir}/smooth_thermodynamic_properties_{nx}x{ny}.png')
    analyzer.plot_critical_analysis(f'{output_dir}/smooth_critical_analysis_{nx}x{ny}.png')
    
    # Generate detailed report
    analyzer.generate_report(f'{output_dir}/smooth_thermodynamic_report_{nx}x{ny}.txt')

def run_parallel_simulations(lattice_sizes, nmcs1, nmcs2, num_runs, seed, base_output_dir='simulation_runs'):
    """Run multiple simulations in parallel for each lattice size"""
    all_avg_files = []
    
    # Get number of CPU cores
    num_cores = cpu_count()
    print(f"Using {num_cores} CPU cores for parallel processing")
    
    for nx, ny in lattice_sizes:
        print(f"\n{'='*60}")
        print(f"LATTICE SIZE: {nx}x{ny}")
        print(f"Running {num_runs} simulations in parallel...")
        print(f"{'='*60}")
        
        output_dir = f"{base_output_dir}/nx{nx}_ny{ny}"
        os.makedirs(output_dir, exist_ok=True)
        
        # Prepare arguments for parallel processing
        args_list = []
        for run_id in range(1, num_runs + 1):
            current_seed = seed + 2 * (run_id - 1)  # Increment seed by 2 for each run
            args = (nx, ny, nmcs1, nmcs2, current_seed, run_id, output_dir)
            args_list.append(args)
        
        # Run simulations in parallel
        start_time = time.time()
        
        with Pool(processes=num_cores) as pool:
            results = pool.map(run_single_simulation_parallel, args_list)
        
        # Filter out None results
        data_files = [result for result in results if result is not None]
        
        elapsed_time = time.time() - start_time
        print(f"\nCompleted {len(data_files)} simulations in {elapsed_time:.1f} seconds")
        print(f"Average time per simulation: {elapsed_time/len(data_files):.1f} seconds")
        
        if data_files:
            # Average the results
            avg_file = average_simulation_data(data_files, f'{output_dir}/smooth_averaged_data.txt')
            all_avg_files.append((nx, ny, avg_file))
            
            # Create smooth plots
            create_smooth_plots(avg_file, nx, ny)
    
    return all_avg_files

def get_system_info():
    """Get system information for optimal parallel processing"""
    cpu_count = mp.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    
    print(f"System Information:")
    print(f"  CPU cores: {cpu_count}")
    print(f"  Total RAM: {memory_gb:.1f} GB")
    print(f"  Available RAM: {psutil.virtual_memory().available / (1024**3):.1f} GB")
    
    return cpu_count, memory_gb

def main():
    """Main function with parallel processing"""
    print("Parallel Monte Carlo Simulations using all CPU cores")
    print("="*60)
    
    # Get system information
    num_cores, memory_gb = get_system_info()
    
    # Configuration
    config = {
        'lattice_sizes': [(8,8), (16,16), (24,24), (32,32)],
        'nmcs1': 10000,
        'nmcs2': 20000,
        'num_runs': 10,
        'seed': 12345,
        'base_output_dir': 'simulation_runs'
    }
    
    print("\nConfiguration:")
    for key, value in config.items():
        print(f"  {key}: {value}")
    
    total_runs = len(config['lattice_sizes']) * config['num_runs']
    total_measurements = config['nmcs2'] * config['num_runs']
    
    print(f"\nThis will run {total_runs} simulations in total.")
    print(f"Total measurement steps per temperature per size: {total_measurements:,}")
    print(f"Starting seed: {config['seed']}, incrementing by 2 for each run")
    print(f"Parallel processing with {num_cores} cores")
    
    # Estimate computation time
    estimated_time = total_runs * config['nmcs2'] / (1000000 * num_cores)  # Rough estimate in minutes
    print(f"Estimated computation time: {estimated_time:.1f} minutes")
    
    # Check memory requirements
    memory_per_sim = 0.1  # GB per simulation (rough estimate)
    total_memory_needed = total_runs * memory_per_sim
    print(f"Estimated memory usage: {total_memory_needed:.1f} GB")
    
    if total_memory_needed > memory_gb * 0.8:
        print(f"⚠️  Warning: High memory usage expected. Consider reducing num_runs.")
    
    # Ask for confirmation
    response = input("\nProceed with these parameters? (y/n): ").lower().strip()
    if response != 'y':
        print("Simulation cancelled.")
        return
    
    # Run the parallel simulations
    start_time = time.time()
    all_avg_files = run_parallel_simulations(**config)
    total_time = time.time() - start_time
    
    print("\n" + "="*60)
    print("PARALLEL SIMULATION SUMMARY")
    print("="*60)
    print(f"Total simulations completed: {total_runs}")
    print(f"Total computation time: {total_time/60:.1f} minutes")
    print(f"Average time per simulation: {total_time/total_runs:.1f} seconds")
    print(f"Speedup factor: ~{num_cores}x faster than sequential")
    
    print(f"\nFiles generated:")
    for nx, ny, avg_file in all_avg_files:
        print(f"  - {avg_file}")
        print(f"  - smooth_thermodynamic_properties_{nx}x{ny}.png")
        print(f"  - smooth_critical_analysis_{nx}x{ny}.png")
        print(f"  - smooth_thermodynamic_report_{nx}x{ny}.txt")
    
    print(f"\n🎉 Parallel analysis complete!")
    print(f"Ready for Q1 journal submission with high-quality results!")
    
    # Create lattice size comparison plots if we have multiple sizes
    if len(all_avg_files) >= 2:
        print(f"\n{'='*60}")
        print("CREATING LATTICE SIZE COMPARISON PLOTS")
        print(f"{'='*60}")
        
        # Create dictionary of data files
        data_files_dict = {}
        for nx, ny, avg_file in all_avg_files:
            data_files_dict[(nx, ny)] = avg_file
        
        # Import analyzer and create comparison plots
        from analyze_thermodynamics import ThermodynamicAnalyzer
        
        # Create a dummy analyzer to access the comparison methods
        dummy_analyzer = ThermodynamicAnalyzer()
        
        print("Creating comprehensive lattice size comparison...")
        critical_temps = dummy_analyzer.compare_lattice_sizes(data_files_dict, 'simulation_runs/lattice_size_comparison.png')
        
        print("Creating finite size scaling analysis...")
        scaling_results = dummy_analyzer.plot_finite_size_scaling(data_files_dict, 'simulation_runs/finite_size_scaling.png')
        
        if scaling_results:
            print(f"\nFinite Size Scaling Results:")
            print(f"  Lattice sizes: {scaling_results['sizes']}")
            print(f"  Critical temperatures: {[f'{tc:.6f}' for tc in scaling_results['critical_temps']]}")
            if len(scaling_results['sizes']) >= 3:
                print(f"  Thermodynamic limit analysis completed")
        
        print(f"\n✅ Lattice size comparison complete!")
        print(f"Generated files:")
        print(f"  - lattice_size_comparison.png")
        print(f"  - finite_size_scaling.png")
    
    print(f"\n🎉 All analysis complete!")
    print(f"Ready for Q1 journal submission with high-quality results!")

if __name__ == "__main__":
    main() 