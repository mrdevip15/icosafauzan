import os
import subprocess
import multiprocessing
import time
from datetime import datetime

# Parameters
LATTICE_SIZES = [8, 16, 24, 32]
MONTE_CARLO_STEPS_EQUILIBRATION = 10000  # nmcs1
MONTE_CARLO_STEPS_MEASUREMENT = 20000    # nmcs2
RANDOM_SEED = 12345                      # iri

def run_simulation(nx):
    """Run the cubic3lay simulation for a given lattice size"""
    print(f"Starting simulation for lattice size {nx}x{nx}...")
    start_time = time.time()
    
    # Create input for the simulation
    input_str = f"{nx} {nx} {MONTE_CARLO_STEPS_EQUILIBRATION} {MONTE_CARLO_STEPS_MEASUREMENT} {RANDOM_SEED}"
    
    # Run the simulation
    process = subprocess.Popen("cubic3lay.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, text=True)
    
    # Send input to the process
    stdout_data, stderr_data = process.communicate(input=input_str)
    
    # Check if simulation was successful
    if process.returncode != 0:
        print(f"Error running simulation for {nx}x{nx}: {stderr_data}")
        return False
    
    elapsed_time = time.time() - start_time
    print(f"Simulation for {nx}x{nx} completed in {elapsed_time:.2f} seconds")
    
    # Check if output directories were created
    base_dir = "simulation_runs"
    if os.path.exists(base_dir):
        dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
        if dirs:
            latest_dir = max(dirs, key=lambda d: os.path.getmtime(os.path.join(base_dir, d)))
            size_dir = os.path.join(base_dir, latest_dir, f"nx{nx}_ny{nx}")
            if os.path.exists(size_dir):
                print(f"Output directory created: {size_dir}")
                layer_files = [f for f in os.listdir(size_dir) if f.startswith("layer_")]
                print(f"Layer files created: {', '.join(layer_files)}")
    
    return True

def run_all_simulations():
    """Run simulations for all lattice sizes in parallel"""
    print(f"Starting parallel simulations for lattice sizes: {LATTICE_SIZES}")
    print(f"Current time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create a pool of workers
    num_cores = min(multiprocessing.cpu_count(), len(LATTICE_SIZES))
    print(f"Using {num_cores} CPU cores")
    
    # Run simulations in parallel
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(run_simulation, LATTICE_SIZES)
    
    # Check results
    successful = [size for size, success in zip(LATTICE_SIZES, results) if success]
    failed = [size for size, success in zip(LATTICE_SIZES, results) if not success]
    
    print("\nSimulation Summary:")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if not failed:
        print("\nAll simulations completed successfully. You can now run analyze_simulation_results.py to analyze the results.")
    else:
        print("\nSome simulations failed. Please check the output for errors.")

if __name__ == "__main__":
    run_all_simulations() 