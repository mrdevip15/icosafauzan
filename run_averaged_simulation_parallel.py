import os
import subprocess
import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import datetime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_simulation(params):
    """Run the cubic3lay simulation for a given lattice size"""
    nx, ny, run_index = params
    print(f"Starting simulation {run_index+1}/3 for lattice size {nx}x{ny}...")
    start_time = time.time()
    
    # Create input for the simulation
    input_str = f"{nx} {ny} 10000 20000 {12345 + run_index}"  # Different random seeds
    
    # Run the simulation
    process = subprocess.Popen("cubic3lay.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, text=True)
    
    # Send input to the process
    stdout_data, stderr_data = process.communicate(input=input_str)
    
    # Check if simulation was successful
    if process.returncode != 0:
        print(f"Error running simulation {run_index+1}: {stderr_data}")
        return None
    
    elapsed_time = time.time() - start_time
    print(f"Simulation {run_index+1}/3 completed in {elapsed_time:.2f} seconds")
    
    # Find the most recent simulation directory
    sim_dirs = []
    for dirname in os.listdir("simulation_runs"):
        full_path = os.path.join("simulation_runs", dirname)
        if os.path.isdir(full_path):
            sim_dirs.append(full_path)
    
    if not sim_dirs:
        print("No simulation directories found.")
        return None
    
    # Sort by modification time (most recent first)
    sim_dirs.sort(key=lambda d: os.path.getmtime(d), reverse=True)
    latest_dir = sim_dirs[0]
    
    # Load data from each layer
    size_dir = os.path.join(latest_dir, f"nx{nx}_ny{ny}")
    if not os.path.exists(size_dir):
        print(f"Directory not found: {size_dir}")
        return None
    
    data = []
    for layer in range(1, 4):
        layer_file = os.path.join(size_dir, f"layer_{layer}.txt")
        if os.path.exists(layer_file):
            layer_data = np.loadtxt(layer_file, comments='#')
            data.append(layer_data)
        else:
            print(f"Layer {layer} data not found")
            return None
    
    return data

def run_parallel_simulations(nx, ny, num_runs=3):
    """Run multiple simulations in parallel"""
    params_list = [(nx, ny, i) for i in range(num_runs)]
    all_runs_data = []
    
    # Determine number of cores to use (use at most num_runs cores)
    num_cores = min(multiprocessing.cpu_count(), num_runs)
    print(f"Running {num_runs} simulations in parallel using {num_cores} CPU cores")
    
    # Run simulations in parallel
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        # Submit all tasks
        future_to_params = {executor.submit(run_simulation, params): params for params in params_list}
        
        # Process results as they complete
        for future in as_completed(future_to_params):
            params = future_to_params[future]
            run_index = params[2]
            try:
                data = future.result()
                if data is not None:
                    all_runs_data.append((run_index, data))
                else:
                    print(f"Run {run_index+1} failed to return data")
                    return None
            except Exception as exc:
                print(f"Run {run_index+1} generated an exception: {exc}")
                return None
    
    # Sort results by run index to ensure consistent ordering
    all_runs_data.sort(key=lambda x: x[0])
    return [data for _, data in all_runs_data]

def average_simulations(all_runs_data):
    """Average the results from multiple simulation runs"""
    if not all_runs_data:
        return None
    
    # Average the data across runs
    avg_data = []
    for layer in range(3):  # 3 layers
        # Initialize with the first run's temperature values (column 0)
        avg_layer_data = np.copy(all_runs_data[0][layer])
        
        # Average columns 1-7 (thermodynamic properties) across runs
        for col in range(1, 8):
            for run in range(1, len(all_runs_data)):
                avg_layer_data[:, col] += all_runs_data[run][layer][:, col]
            avg_layer_data[:, col] /= len(all_runs_data)
        
        avg_data.append(avg_layer_data)
    
    return avg_data

def save_averaged_data(avg_data, nx, ny):
    """Save the averaged data to files"""
    # Create output directory
    output_dir = f"averaged_results_nx{nx}_ny{ny}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save data for each layer
    for layer in range(3):
        output_file = os.path.join(output_dir, f"avg_layer_{layer+1}.txt")
        np.savetxt(output_file, avg_data[layer], 
                  header="Temperature M^2 M^4 G^2 G^4 Energy Cv Corr", 
                  fmt='%13.6e')
        print(f"Saved averaged data for layer {layer+1} to {output_file}")
    
    return output_dir

def plot_averaged_results(avg_data, nx, ny, output_dir):
    """Plot the averaged thermodynamic properties"""
    # Define properties and their corresponding column indices
    properties = [
        ("Magnetization²", 1),
        ("Specific Heat", 6),
        ("Energy", 5),
        ("Correlation Length", 7),
        ("Binder Cumulant", None)  # Special case, calculated from M² and M⁴
    ]
    
    # Theoretical critical temperature for 8-state Potts model
    tc_theory = 1.0 / np.log(1 + np.sqrt(8))  # ~0.751
    
    # Get temperature values (same for all layers)
    temps = avg_data[0][:, 0]
    
    # Create a figure with subplots
    plt.figure(figsize=(18, 12))
    
    for i, (prop, col) in enumerate(properties):
        plt.subplot(2, 3, i+1)
        
        if prop == "Binder Cumulant":
            # Calculate Binder cumulant for each layer
            for layer in range(3):
                m2 = avg_data[layer][:, 1]  # M²
                m4 = avg_data[layer][:, 2]  # M⁴
                binder = 1 - m4 / (3 * m2 * m2)
                plt.plot(temps, binder, 'o-', label=f"Layer {layer+1}")
                
            # Also calculate average over all layers
            m2_avg = np.mean([avg_data[layer][:, 1] for layer in range(3)], axis=0)
            m4_avg = np.mean([avg_data[layer][:, 2] for layer in range(3)], axis=0)
            binder_avg = 1 - m4_avg / (3 * m2_avg * m2_avg)
            plt.plot(temps, binder_avg, 'k-', linewidth=2, label="All Layers")
        else:
            # Plot property for each layer
            for layer in range(3):
                plt.plot(temps, avg_data[layer][:, col], 'o-', label=f"Layer {layer+1}")
                
            # Also plot average over all layers
            avg_property = np.mean([avg_data[layer][:, col] for layer in range(3)], axis=0)
            plt.plot(temps, avg_property, 'k-', linewidth=2, label="All Layers")
        
        # Mark theoretical critical temperature
        plt.axvline(x=tc_theory, color='r', linestyle='--', alpha=0.7, 
                   label=f"Tc (theory) = {tc_theory:.3f}")
        
        plt.xlabel("Temperature")
        plt.ylabel(prop)
        plt.title(f"{prop} vs Temperature ({nx}x{ny} lattice)")
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"averaged_properties_nx{nx}_ny{ny}.png"), dpi=300)
    print(f"Saved plot to {os.path.join(output_dir, f'averaged_properties_nx{nx}_ny{ny}.png')}")
    
    # Also create individual plots for each property
    for prop, col in properties:
        plt.figure(figsize=(10, 6))
        
        if prop == "Binder Cumulant":
            # Calculate Binder cumulant
            for layer in range(3):
                m2 = avg_data[layer][:, 1]
                m4 = avg_data[layer][:, 2]
                binder = 1 - m4 / (3 * m2 * m2)
                plt.plot(temps, binder, 'o-', label=f"Layer {layer+1}")
                
            # Also calculate average over all layers
            m2_avg = np.mean([avg_data[layer][:, 1] for layer in range(3)], axis=0)
            m4_avg = np.mean([avg_data[layer][:, 2] for layer in range(3)], axis=0)
            binder_avg = 1 - m4_avg / (3 * m2_avg * m2_avg)
            plt.plot(temps, binder_avg, 'k-', linewidth=2, label="All Layers")
        else:
            # Plot property for each layer
            for layer in range(3):
                plt.plot(temps, avg_data[layer][:, col], 'o-', label=f"Layer {layer+1}")
                
            # Also plot average over all layers
            avg_property = np.mean([avg_data[layer][:, col] for layer in range(3)], axis=0)
            plt.plot(temps, avg_property, 'k-', linewidth=2, label="All Layers")
        
        # Mark theoretical critical temperature
        plt.axvline(x=tc_theory, color='r', linestyle='--', alpha=0.7, 
                   label=f"Tc (theory) = {tc_theory:.3f}")
        
        plt.xlabel("Temperature")
        plt.ylabel(prop)
        plt.title(f"{prop} vs Temperature ({nx}x{ny} lattice)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        filename = f"{prop.lower().replace('²', '2').replace(' ', '_')}_nx{nx}_ny{ny}.png"
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        print(f"Saved {prop} plot to {os.path.join(output_dir, filename)}")
        plt.close()

def main():
    # Lattice size
    nx = 8
    ny = 8
    num_runs = 3
    
    print(f"Running {nx}x{ny} simulation with {num_runs} parallel runs for averaging")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Run simulations in parallel
    all_runs_data = run_parallel_simulations(nx, ny, num_runs)
    if all_runs_data is None:
        print("Failed to complete all simulation runs")
        return
    
    # Average the results
    avg_data = average_simulations(all_runs_data)
    if avg_data is None:
        print("Failed to average simulation results")
        return
    
    # Save averaged data
    output_dir = save_averaged_data(avg_data, nx, ny)
    
    # Plot results
    plot_averaged_results(avg_data, nx, ny, output_dir)
    
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("All tasks completed successfully")

if __name__ == "__main__":
    # This is required for Windows to use multiprocessing
    multiprocessing.freeze_support()
    main() 