import os
import subprocess
import sys
import time

def run_simulation(nx, ny):
    """Run the cubic3lay simulation for a given lattice size"""
    print(f"Running simulation for lattice size {nx}x{ny}...")
    print(f"Starting time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create input for the simulation
    input_str = f"{nx} {ny} 10000 20000 12345"
    
    # Run the simulation
    process = subprocess.Popen("cubic3lay.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, text=True)
    
    # Send input to the process
    print(f"Sending input: {input_str}")
    stdout_data, stderr_data = process.communicate(input=input_str)
    
    # Print first few lines of output to show progress
    print("\nOutput preview:")
    output_lines = stdout_data.splitlines()
    for i, line in enumerate(output_lines[:10]):  # Show first 10 lines
        print(line)
    
    if len(output_lines) > 10:
        print(f"... ({len(output_lines) - 10} more lines)")
    
    # Print the last few lines to confirm completion
    if len(output_lines) > 5:
        print("\nLast 5 lines of output:")
        for line in output_lines[-5:]:
            print(line)
    
    if process.returncode != 0:
        print(f"Error running simulation: {stderr_data}")
        return False
    
    print(f"\nSimulation for {nx}x{ny} completed successfully")
    print(f"Ending time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Check if output directories were created
    base_dir = "simulation_runs"
    if os.path.exists(base_dir):
        dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
        if dirs:
            latest_dir = max(dirs, key=lambda d: os.path.getmtime(os.path.join(base_dir, d)))
            size_dir = os.path.join(base_dir, latest_dir, f"nx{nx}_ny{ny}")
            if os.path.exists(size_dir):
                print(f"Output directory created: {size_dir}")
                layer_files = [f for f in os.listdir(size_dir) if f.startswith("layer_")]
                print(f"Layer files created: {', '.join(layer_files)}")
    
    return True

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_single_simulation.py <lattice_size1> [<lattice_size2> ...]")
        sys.exit(1)
    
    for arg in sys.argv[1:]:
        size = int(arg)
        success = run_simulation(size, size)
        if not success:
            print(f"Failed to run simulation for size {size}x{size}")
        print("\n" + "="*50 + "\n") 