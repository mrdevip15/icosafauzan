import os
import subprocess
import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import datetime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import queue
import threading

# Configuration
LATTICE_SIZES = [8, 16, 24, 32]
RUNS_PER_LATTICE = 3
MONTE_CARLO_STEPS_EQUILIBRATION = 10000
MONTE_CARLO_STEPS_MEASUREMENT = 20000
BASE_RANDOM_SEED = 12345

def run_simulation(params):
    """Run the cubic3lay simulation for a given lattice size"""
    nx, ny, run_index = params
    print(f"Starting simulation for lattice size {nx}x{ny}, run {run_index+1}/{RUNS_PER_LATTICE}")
    start_time = time.time()
    
    # Create input for the simulation
    input_str = f"{nx} {ny} {MONTE_CARLO_STEPS_EQUILIBRATION} {MONTE_CARLO_STEPS_MEASUREMENT} {BASE_RANDOM_SEED + run_index}"
    
    # Run the simulation
    process = subprocess.Popen("cubic3lay.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, text=True)
    
    # Send input to the process
    stdout_data, stderr_data = process.communicate(input=input_str)
    
    # Check if simulation was successful
    if process.returncode != 0:
        print(f"Error running simulation for {nx}x{ny}, run {run_index+1}: {stderr_data}")
        return None
    
    elapsed_time = time.time() - start_time
    print(f"Completed simulation for {nx}x{ny}, run {run_index+1} in {elapsed_time:.2f} seconds")
    
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
            print(f"Layer {layer} data not found for {nx}x{ny}, run {run_index+1}")
            return None
    
    return (nx, ny, run_index, data)

def process_results_thread(result_queue, completed_runs, lock):
    """Process completed simulation results"""
    # Dictionary to store results for each lattice size
    lattice_results = {size: [] for size in LATTICE_SIZES}
    
    while True:
        try:
            # Get a result from the queue (timeout to check if we should exit)
            result = result_queue.get(timeout=1)
            
            if result is None:
                # None is our signal to exit
                break
                
            nx, ny, run_index, data = result
            
            if data is not None:
                # Add the result to the appropriate list
                lattice_results[(nx, ny)].append((run_index, data))
                
                # Update the completed runs count
                with lock:
                    completed_runs[(nx, ny)] += 1
                    print(f"Progress: {nx}x{ny} - {completed_runs[(nx, ny)]}/{RUNS_PER_LATTICE} runs completed")
                
                # Check if all runs for this lattice size are complete
                if completed_runs[(nx, ny)] == RUNS_PER_LATTICE:
                    print(f"All runs for lattice size {nx}x{ny} completed. Processing results...")
                    process_lattice_results(nx, ny, lattice_results[(nx, ny)])
            
            # Mark the task as done
            result_queue.task_done()
            
        except queue.Empty:
            # Queue is empty, check if we should exit
            with lock:
                total_completed = sum(completed_runs.values())
                total_needed = len(LATTICE_SIZES) * RUNS_PER_LATTICE
                if total_completed == total_needed:
                    break
                    
    print("Result processing thread exiting")

def process_lattice_results(nx, ny, results):
    """Process and save results for a specific lattice size"""
    # Sort results by run index
    results.sort(key=lambda x: x[0])
    
    # Extract just the data
    all_runs_data = [data for _, data in results]
    
    # Average the data across runs
    avg_data = average_simulations(all_runs_data)
    
    # Save the averaged data
    output_dir = save_averaged_data(avg_data, nx, ny)
    
    # Plot the results
    plot_averaged_results(avg_data, nx, ny, output_dir)

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
    print(f"Starting parallel simulations for all lattice sizes")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create a queue for tasks and results
    task_queue = queue.Queue()
    result_queue = queue.Queue()
    
    # Create a dictionary to track completed runs
    completed_runs = {(nx, ny): 0 for nx, ny in [(size, size) for size in LATTICE_SIZES]}
    lock = threading.Lock()
    
    # Add all tasks to the queue
    for size in LATTICE_SIZES:
        for run in range(RUNS_PER_LATTICE):
            task_queue.put((size, size, run))
    
    # Start a thread to process results
    result_thread = threading.Thread(target=process_results_thread, 
                                    args=(result_queue, completed_runs, lock))
    result_thread.start()
    
    # Determine number of workers (CPU cores)
    num_cores = multiprocessing.cpu_count()
    print(f"Using {num_cores} CPU cores for parallel processing")
    
    # Create a pool of workers
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        # Submit initial batch of tasks (one per core)
        futures = {}
        for _ in range(min(num_cores, task_queue.qsize())):
            if not task_queue.empty():
                params = task_queue.get()
                future = executor.submit(run_simulation, params)
                futures[future] = params
        
        # Process tasks as they complete and add new ones from the queue
        while futures:
            # Wait for the next task to complete
            done, not_done = as_completed(futures, timeout=None, return_when='FIRST_COMPLETED')
            
            for future in done:
                params = futures[future]
                try:
                    # Get the result and put it in the result queue
                    result = future.result()
                    if result is not None:
                        result_queue.put(result)
                    
                    # Remove the completed future
                    del futures[future]
                    
                    # Add a new task if there are any left
                    if not task_queue.empty():
                        new_params = task_queue.get()
                        new_future = executor.submit(run_simulation, new_params)
                        futures[new_future] = new_params
                        
                except Exception as exc:
                    nx, ny, run_index = params
                    print(f"Simulation for {nx}x{ny}, run {run_index+1} generated an exception: {exc}")
                    
                    # Remove the failed future
                    del futures[future]
                    
                    # Add a new task if there are any left
                    if not task_queue.empty():
                        new_params = task_queue.get()
                        new_future = executor.submit(run_simulation, new_params)
                        futures[new_future] = new_params
    
    # Signal the result thread to exit
    result_queue.put(None)
    
    # Wait for the result thread to finish
    result_thread.join()
    
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("All simulations and analysis completed")
    
    # Create a final summary plot comparing all lattice sizes
    create_combined_plots()

def create_combined_plots():
    """Create plots comparing results from all lattice sizes"""
    # Check which lattice sizes have completed results
    completed_sizes = []
    for size in LATTICE_SIZES:
        output_dir = f"averaged_results_nx{size}_ny{size}"
        if os.path.exists(output_dir):
            completed_sizes.append(size)
    
    if not completed_sizes:
        print("No completed results found for any lattice size")
        return
    
    print(f"Creating combined plots for lattice sizes: {completed_sizes}")
    
    # Load data for each completed lattice size
    all_data = {}
    for size in completed_sizes:
        output_dir = f"averaged_results_nx{size}_ny{size}"
        
        # Load data for all three layers
        size_data = []
        for layer in range(1, 4):
            file_path = os.path.join(output_dir, f"avg_layer_{layer}.txt")
            if os.path.exists(file_path):
                layer_data = np.loadtxt(file_path)
                size_data.append(layer_data)
            else:
                print(f"Missing data file for lattice size {size}, layer {layer}")
                break
        
        if len(size_data) == 3:
            all_data[size] = size_data
    
    if not all_data:
        print("No valid data found for any lattice size")
        return
    
    # Create plots comparing different lattice sizes
    properties = [
        ("Magnetization²", 1),
        ("Specific Heat", 6),
        ("Energy", 5),
        ("Correlation Length", 7),
        ("Binder Cumulant", None)
    ]
    
    # Theoretical critical temperature
    tc_theory = 1.0 / np.log(1 + np.sqrt(8))  # ~0.751
    
    # Create a figure for each property
    for prop_name, col_idx in properties:
        plt.figure(figsize=(12, 8))
        
        for size in all_data.keys():
            # Average over all three layers
            if prop_name == "Binder Cumulant":
                # Calculate Binder cumulant from M² and M⁴
                m2_avg = np.mean([all_data[size][layer][:, 1] for layer in range(3)], axis=0)
                m4_avg = np.mean([all_data[size][layer][:, 2] for layer in range(3)], axis=0)
                values = 1 - m4_avg / (3 * m2_avg * m2_avg)
            else:
                # Average the property over all three layers
                values = np.mean([all_data[size][layer][:, col_idx] for layer in range(3)], axis=0)
            
            # Get temperature values
            temps = all_data[size][0][:, 0]
            
            # Plot the data
            plt.plot(temps, values, 'o-', label=f"L = {size}")
        
        # Mark theoretical critical temperature
        plt.axvline(x=tc_theory, color='r', linestyle='--', alpha=0.7, 
                   label=f"Tc (theory) = {tc_theory:.3f}")
        
        plt.xlabel("Temperature", fontsize=14)
        plt.ylabel(prop_name, fontsize=14)
        plt.title(f"{prop_name} vs Temperature for Different Lattice Sizes", fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # Save the figure
        filename = f"combined_{prop_name.lower().replace('²', '2').replace(' ', '_')}.png"
        plt.savefig(filename, dpi=300)
        print(f"Saved combined plot to {filename}")
        plt.close()
    
    # Create finite-size scaling plot for critical temperature
    create_finite_size_scaling_plot(all_data)

def create_finite_size_scaling_plot(all_data):
    """Create finite-size scaling plot for critical temperature"""
    # Methods to estimate critical temperature
    methods = {
        "Specific Heat Peak": 6,
        "Binder Cumulant Crossing": None
    }
    
    # Theoretical critical temperature
    tc_theory = 1.0 / np.log(1 + np.sqrt(8))  # ~0.751
    
    # Find critical temperatures for each lattice size using different methods
    tc_values = {method: {} for method in methods}
    
    for size in all_data.keys():
        # Specific Heat Peak method
        cv_avg = np.mean([all_data[size][layer][:, 6] for layer in range(3)], axis=0)
        temps = all_data[size][0][:, 0]
        peak_idx = np.argmax(cv_avg)
        tc_values["Specific Heat Peak"][size] = temps[peak_idx]
        
        # Binder Cumulant method (if we have at least two lattice sizes)
        if len(all_data) >= 2 and size == max(all_data.keys()):
            # We'll compute this after collecting data for all sizes
            pass
    
    # Compute Binder cumulant crossing points if we have multiple lattice sizes
    if len(all_data) >= 2:
        # Calculate Binder cumulant for each lattice size
        binder_data = {}
        for size in all_data.keys():
            m2_avg = np.mean([all_data[size][layer][:, 1] for layer in range(3)], axis=0)
            m4_avg = np.mean([all_data[size][layer][:, 2] for layer in range(3)], axis=0)
            binder = 1 - m4_avg / (3 * m2_avg * m2_avg)
            binder_data[size] = (all_data[size][0][:, 0], binder)  # (temps, binder)
        
        # Find crossing points between all pairs of lattice sizes
        crossings = []
        sizes = list(binder_data.keys())
        for i in range(len(sizes)):
            for j in range(i+1, len(sizes)):
                size1, size2 = sizes[i], sizes[j]
                temps1, binder1 = binder_data[size1]
                temps2, binder2 = binder_data[size2]
                
                # Interpolate to a common temperature grid
                common_temps = np.linspace(max(temps1.min(), temps2.min()),
                                          min(temps1.max(), temps2.max()),
                                          1000)
                binder1_interp = np.interp(common_temps, temps1, binder1)
                binder2_interp = np.interp(common_temps, temps2, binder2)
                
                # Find where curves cross by finding where their difference changes sign
                diff = binder1_interp - binder2_interp
                for k in range(len(diff)-1):
                    if diff[k] * diff[k+1] <= 0:  # Sign change
                        # Linear interpolation to find crossing point
                        t1, t2 = common_temps[k], common_temps[k+1]
                        d1, d2 = diff[k], diff[k+1]
                        tc = t1 - d1 * (t2 - t1) / (d2 - d1)
                        crossings.append(tc)
                        print(f"Found crossing between L={size1} and L={size2} at T={tc:.4f}")
        
        # Average crossing points
        if crossings:
            avg_crossing = np.mean(crossings)
            for size in all_data.keys():
                tc_values["Binder Cumulant Crossing"][size] = avg_crossing
    
    # Create finite-size scaling plot
    plt.figure(figsize=(10, 8))
    
    for method, tc_dict in tc_values.items():
        if tc_dict:  # Check if we have data for this method
            sizes = np.array(list(tc_dict.keys()))
            tcs = np.array(list(tc_dict.values()))
            
            # Plot data points
            plt.scatter(1/sizes, tcs, label=f"{method} data")
            
            # Fit to Tc(L) = Tc(∞) + a/L^(1/ν)
            # For 8-state Potts model, ν = 1
            if len(sizes) >= 2:  # Need at least 2 points for fitting
                try:
                    from scipy.optimize import curve_fit
                    
                    def scaling_func(x, tc_inf, a):
                        return tc_inf + a / x
                    
                    params, _ = curve_fit(scaling_func, sizes, tcs)
                    tc_inf, a = params
                    
                    # Plot fit
                    x_fit = np.linspace(1/max(sizes), 1/min(sizes), 100)
                    plt.plot(x_fit, scaling_func(1/x_fit, tc_inf, a), '--', 
                            label=f"{method}: Tc(∞) = {tc_inf:.4f}")
                    
                    print(f"{method}: Extrapolated Tc(∞) = {tc_inf:.4f}, "
                         f"difference from theory: {abs(tc_inf - tc_theory):.4f}")
                except Exception as e:
                    print(f"Fitting failed for {method}: {e}")
    
    # Add theoretical value
    plt.axhline(y=tc_theory, color='r', linestyle='-', 
               label=f"Theoretical Tc = {tc_theory:.4f}")
    
    plt.xlabel("1/L", fontsize=14)
    plt.ylabel("Critical Temperature Tc(L)", fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.title("Finite Size Scaling of Critical Temperature", fontsize=16)
    
    # Save the figure
    plt.savefig("finite_size_scaling.png", dpi=300)
    print("Saved finite-size scaling plot to finite_size_scaling.png")
    plt.close()

if __name__ == "__main__":
    # This is required for Windows to use multiprocessing
    multiprocessing.freeze_support()
    main() 