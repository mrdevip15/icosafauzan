import os
import time

def check_simulation_progress():
    """Check the progress of the simulations"""
    print("Checking simulation progress...")
    
    base_dir = "simulation_runs"
    if not os.path.exists(base_dir):
        print("No simulation_runs directory found. Simulations may not have started yet.")
        return
    
    # Get all simulation directories
    sim_dirs = []
    for dirname in os.listdir(base_dir):
        full_path = os.path.join(base_dir, dirname)
        if os.path.isdir(full_path):
            sim_dirs.append(full_path)
    
    if not sim_dirs:
        print("No simulation directories found.")
        return
    
    # Sort by modification time (most recent first)
    sim_dirs.sort(key=lambda d: os.path.getmtime(d), reverse=True)
    
    # Check the most recent simulation directory
    latest_dir = sim_dirs[0]
    print(f"Most recent simulation directory: {latest_dir}")
    print(f"Created: {time.ctime(os.path.getctime(latest_dir))}")
    print(f"Last modified: {time.ctime(os.path.getmtime(latest_dir))}")
    
    # Check for lattice size directories
    lattice_sizes = []
    for size_dir in os.listdir(latest_dir):
        if size_dir.startswith("nx"):
            lattice_sizes.append(size_dir)
    
    print(f"Found lattice sizes: {', '.join(lattice_sizes) if lattice_sizes else 'None'}")
    
    # Check each lattice size directory
    for size_dir in lattice_sizes:
        full_size_path = os.path.join(latest_dir, size_dir)
        print(f"\nChecking {size_dir}:")
        
        # Check for layer files
        layer_files = [f for f in os.listdir(full_size_path) if f.startswith("layer_")]
        print(f"  Layer files: {', '.join(layer_files) if layer_files else 'None'}")
        
        # Check file sizes and last modification times
        for layer_file in layer_files:
            file_path = os.path.join(full_size_path, layer_file)
            file_size = os.path.getsize(file_path) / 1024  # KB
            mod_time = time.ctime(os.path.getmtime(file_path))
            
            # Check if file is still being written to
            is_active = (time.time() - os.path.getmtime(file_path)) < 60  # Modified in last minute
            status = "ACTIVE" if is_active else "Completed"
            
            print(f"  - {layer_file}: {file_size:.2f} KB, Last modified: {mod_time}, Status: {status}")
            
            # Check content of the file
            if os.path.exists(file_path):
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    data_lines = [l for l in lines if not l.startswith('#')]
                    print(f"    Data lines: {len(data_lines)}")
                    
                    # Show the last few data points
                    if data_lines:
                        print(f"    Last data point: {data_lines[-1].strip()}")
    
    # Check for other lattice sizes that might be in other simulation directories
    all_sizes = set()
    for sim_dir in sim_dirs:
        for item in os.listdir(sim_dir):
            if item.startswith("nx"):
                all_sizes.add(item)
    
    print("\nAll lattice sizes found across all simulations:")
    for size in sorted(all_sizes):
        print(f"- {size}")

if __name__ == "__main__":
    check_simulation_progress() 