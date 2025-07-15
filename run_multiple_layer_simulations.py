#!/usr/bin/env python3
"""
Run multiple layer-separated simulations with different lattice sizes.
"""

import subprocess
import multiprocessing
import os
import time
import sys

def run_single_simulation(params):
    """Run a single simulation with given parameters."""
    nx, ny, nmcs1, nmcs2, seed = params
    
    executable = "./cubic3lay"
    if not os.path.exists(executable):
        return f"Error: {executable} not found for {nx}x{ny}"
    
    input_str = f"{nx} {ny} {nmcs1} {nmcs2} {seed}\n"
    
    try:
        process = subprocess.Popen(
            [executable],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate(input=input_str)
        
        if process.returncode == 0:
            return f"✓ Completed {nx}x{ny}x3 simulation (seed={seed})"
        else:
            return f"✗ Failed {nx}x{ny}x3 simulation: {stderr}"
            
    except Exception as e:
        return f"✗ Error in {nx}x{ny}x3 simulation: {e}"

def main():
    # Simulation parameters
    lattice_sizes = [
        (8, 8),
        (16, 16),
        (24, 24),
        (32, 32)
    ]
    
    nmcs1 = 1000  # Equilibration steps
    nmcs2 = 2000 # Measurement steps
    
    # Multiple seeds for statistical analysis
    seeds = [12345, 12346, 12347, 12348, 12349]
    
    print("=== Multiple Layer-separated Quasi-3D Simulations ===")
    print(f"Lattice sizes: {lattice_sizes}")
    print(f"Number of seeds per size: {len(seeds)}")
    print(f"Equilibration steps: {nmcs1}")
    print(f"Measurement steps: {nmcs2}")
    print(f"Total simulations: {len(lattice_sizes) * len(seeds)}")
    print()
    
    # Check if executable exists
    if not os.path.exists("./cubic3lay"):
        print("Error: cubic3lay executable not found!")
        print("Please compile first: gcc -o cubic3lay cubic3lay.c rn32.c -lm")
        sys.exit(1)
    
    # Prepare all simulation parameters
    all_params = []
    for nx, ny in lattice_sizes:
        for seed in seeds:
            all_params.append((nx, ny, nmcs1, nmcs2, seed))
    
    print(f"Starting {len(all_params)} simulations...")
    start_time = time.time()
    
    # Run simulations in parallel
    num_cores = min(multiprocessing.cpu_count(), len(all_params))
    print(f"Using {num_cores} CPU cores")
    
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(run_single_simulation, all_params)
    
    # Print results
    print("\n=== Simulation Results ===")
    for result in results:
        print(result)
    
    # Summary
    successful = sum(1 for r in results if r.startswith("✓"))
    failed = len(results) - successful
    
    end_time = time.time()
    duration = end_time - start_time
    
    print(f"\n=== Summary ===")
    print(f"Total simulations: {len(results)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total time: {duration:.1f} seconds")
    print(f"Average time per simulation: {duration/len(results):.1f} seconds")
    
    # List created directories
    print(f"\n=== Created Output Directories ===")
    for nx, ny in lattice_sizes:
        output_dir = f"simulation_runs/nx{nx}_ny{ny}"
        if os.path.exists(output_dir):
            files = os.listdir(output_dir)
            layer_files = [f for f in files if f.startswith("layer_")]
            print(f"{output_dir}: {len(layer_files)} layer files")
    
    if failed > 0:
        print(f"\nWarning: {failed} simulations failed!")
        sys.exit(1)
    else:
        print("\nAll simulations completed successfully!")

if __name__ == "__main__":
    main()