#!/usr/bin/env python3
"""
Run layer-separated quasi-3D Monte Carlo simulation for 8-state cube spins.
This script runs the modified cubic3lay simulation and organizes output by layers.
"""

import subprocess
import os
import sys
import argparse
from pathlib import Path

def run_simulation(nx, ny, nmcs1, nmcs2, seed):
    """
    Run the cubic3lay simulation with given parameters.
    
    Args:
        nx, ny: Lattice dimensions (always 3 layers in z-direction)
        nmcs1: Equilibration steps
        nmcs2: Measurement steps  
        seed: Random number seed
    """
    
    # Check if executable exists
    executable = "./cubic3lay"
    if not os.path.exists(executable):
        print(f"Error: {executable} not found. Please compile the simulation first:")
        print("gcc -o cubic3lay cubic3lay.c rn32.c -lm")
        return False
    
    # Prepare input string
    input_str = f"{nx} {ny} {nmcs1} {nmcs2} {seed}\n"
    
    print(f"Running simulation: nx={nx}, ny={ny}, nz=3, nmcs1={nmcs1}, nmcs2={nmcs2}, seed={seed}")
    
    try:
        # Run the simulation
        process = subprocess.Popen(
            [executable],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate(input=input_str)
        
        if process.returncode != 0:
            print(f"Simulation failed with return code {process.returncode}")
            print(f"STDERR: {stderr}")
            return False
        
        print("Simulation completed successfully")
        print(f"Output directory: simulation_runs/nx{nx}_ny{ny}/")
        print("Layer files created:")
        
        # List the created files
        output_dir = f"simulation_runs/nx{nx}_ny{ny}"
        if os.path.exists(output_dir):
            for i in range(1, 4):
                layer_file = f"{output_dir}/layer_{i}.txt"
                if os.path.exists(layer_file):
                    print(f"  - {layer_file}")
                    # Print first few lines to verify
                    with open(layer_file, 'r') as f:
                        lines = f.readlines()[:5]
                        print(f"    First few lines of layer_{i}.txt:")
                        for line in lines:
                            print(f"    {line.strip()}")
        
        return True
        
    except Exception as e:
        print(f"Error running simulation: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Run layer-separated quasi-3D Monte Carlo simulation')
    parser.add_argument('--nx', type=int, default=8, help='Lattice width (default: 8)')
    parser.add_argument('--ny', type=int, default=8, help='Lattice height (default: 8)')
    parser.add_argument('--nmcs1', type=int, default=1000, help='Equilibration steps (default: 1000)')
    parser.add_argument('--nmcs2', type=int, default=10000, help='Measurement steps (default: 10000)')
    parser.add_argument('--seed', type=int, default=12345, help='Random seed (default: 12345)')
    
    args = parser.parse_args()
    
    print("=== Layer-separated Quasi-3D Monte Carlo Simulation ===")
    print(f"Parameters: {args.nx}x{args.ny}x3 lattice, {args.nmcs1} equilibration + {args.nmcs2} measurement steps")
    print(f"Random seed: {args.seed}")
    print()
    
    success = run_simulation(args.nx, args.ny, args.nmcs1, args.nmcs2, args.seed)
    
    if success:
        print("\nSimulation completed successfully!")
        print("Each layer's thermodynamic data has been saved to separate files.")
    else:
        print("\nSimulation failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()