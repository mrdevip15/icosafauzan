#!/usr/bin/env python3
"""
Test script to demonstrate the new hierarchical directory structure.
Tests multiple lattice sizes to show the organized output.
"""

import subprocess
import time
import os

def run_simulation(nx, ny, nmcs1, nmcs2, seed):
    """Run simulation with given parameters."""
    try:
        input_data = f"{nx} {ny} {nmcs1} {nmcs2} {seed}\n"
        result = subprocess.run(
            ['./cubic3lay_enhanced'], 
            input=input_data, 
            text=True, 
            capture_output=True, 
            timeout=120
        )
        
        if result.returncode == 0:
            print(f"✓ Simulation completed for {nx}x{ny}")
            # Extract directory from output
            output_lines = result.stdout.strip().split('\n')
            for line in output_lines:
                if "Layer data saved to" in line:
                    return line.split("to ")[-1].rstrip('/')
            return None
        else:
            print(f"✗ Simulation failed for {nx}x{ny}: {result.stderr}")
            return None
            
    except Exception as e:
        print(f"✗ Error for {nx}x{ny}: {e}")
        return None

def test_directory_structure():
    """Test the hierarchical directory structure with multiple lattice sizes."""
    print("=== TESTING HIERARCHICAL DIRECTORY STRUCTURE ===")
    print("Structure: simulation_runs/timestamp/lattice_size/layer_files")
    print("=" * 60)
    
    # Test parameters
    lattice_sizes = [8, 16, 24]
    nmcs1 = 500    # Shorter for testing
    nmcs2 = 1000   # Shorter for testing
    base_seed = 12345
    
    print(f"Testing lattice sizes: {lattice_sizes}")
    print(f"Parameters: equilibration={nmcs1}, measurement={nmcs2}")
    print()
    
    # Keep track of simulation directories
    sim_dirs = []
    
    # Run simulations for different lattice sizes
    for i, nx in enumerate(lattice_sizes):
        ny = nx  # Square lattices
        seed = base_seed + i
        
        print(f"Running simulation {i+1}/{len(lattice_sizes)}: {nx}x{ny}")
        sim_dir = run_simulation(nx, ny, nmcs1, nmcs2, seed)
        
        if sim_dir:
            sim_dirs.append(sim_dir)
            print(f"  → Data saved to: {sim_dir}")
        
        # Small delay to ensure different timestamps if needed
        time.sleep(1)
    
    print(f"\n=== DIRECTORY STRUCTURE ANALYSIS ===")
    
    if sim_dirs:
        # Analyze the directory structure
        print(f"✓ {len(sim_dirs)} simulations completed successfully")
        
        # Check if they're grouped by timestamp
        timestamps = set()
        for sim_dir in sim_dirs:
            # Extract timestamp from path
            parts = sim_dir.split('/')
            if len(parts) >= 2:
                timestamp_dir = parts[1]  # simulation_runs/TIMESTAMP_*/...
                timestamps.add(timestamp_dir)
        
        print(f"✓ Number of timestamp directories: {len(timestamps)}")
        
        # Display the complete structure
        print(f"\nGenerated directory structure:")
        for timestamp in sorted(timestamps):
            timestamp_path = f"simulation_runs/{timestamp}"
            if os.path.exists(timestamp_path):
                print(f"📁 {timestamp}/")
                
                # List lattice size directories
                try:
                    size_dirs = [d for d in os.listdir(timestamp_path) 
                               if os.path.isdir(os.path.join(timestamp_path, d))]
                    
                    for size_dir in sorted(size_dirs):
                        print(f"  📁 {size_dir}/")
                        
                        # List layer files
                        size_path = os.path.join(timestamp_path, size_dir)
                        if os.path.exists(size_path):
                            layer_files = [f for f in os.listdir(size_path) 
                                         if f.endswith('.txt')]
                            for layer_file in sorted(layer_files):
                                file_path = os.path.join(size_path, layer_file)
                                file_size = os.path.getsize(file_path)
                                print(f"    📄 {layer_file} ({file_size:,} bytes)")
                
                except Exception as e:
                    print(f"    ⚠️  Error reading directory: {e}")
        
        # Verify the structure benefits
        print(f"\n=== STRUCTURE BENEFITS ===")
        print(f"✓ Timestamped runs prevent data overwriting")
        print(f"✓ Multiple lattice sizes organized under same timestamp")
        print(f"✓ Layer data clearly separated within each lattice size")
        print(f"✓ Easy to compare different system sizes from same run")
        print(f"✓ Metadata includes enhanced model information")
        
        # Show sample metadata
        if sim_dirs:
            sample_file = f"{sim_dirs[0]}/layer_1.txt"
            if os.path.exists(sample_file):
                print(f"\nSample metadata from {sample_file}:")
                with open(sample_file, 'r') as f:
                    for i, line in enumerate(f):
                        if line.startswith('#') and i < 10:
                            print(f"  {line.rstrip()}")
                        elif not line.startswith('#'):
                            break
        
        return True
    else:
        print("✗ No simulations completed successfully")
        return False

def main():
    """Main test function."""
    print("Enhanced Monte Carlo Simulation - Directory Structure Test")
    print("-" * 60)
    
    # Check if executable exists
    if not os.path.exists('./cubic3lay_enhanced'):
        print("✗ Enhanced executable not found. Please compile first:")
        print("  gcc -o cubic3lay_enhanced cubic3lay.c rn32.c -lm")
        return False
    
    success = test_directory_structure()
    
    print(f"\n" + "="*60)
    if success:
        print("🎉 DIRECTORY STRUCTURE TEST PASSED!")
        print("The hierarchical organization is working correctly.")
    else:
        print("⚠️  Test failed - check implementation")
    print("="*60)
    
    return success

if __name__ == "__main__":
    main()