import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def find_latest_simulation():
    """Find the most recent simulation directory"""
    base_dir = "simulation_runs"
    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} does not exist.")
        return None
    
    # Get all simulation directories
    sim_dirs = []
    for dirname in os.listdir(base_dir):
        full_path = os.path.join(base_dir, dirname)
        if os.path.isdir(full_path):
            sim_dirs.append(full_path)
    
    if not sim_dirs:
        print("No simulation directories found.")
        return None
    
    # Sort by modification time (most recent first)
    sim_dirs.sort(key=lambda d: os.path.getmtime(d), reverse=True)
    return sim_dirs[0]

def load_8x8_data(latest_dir):
    """Load data for 8x8 lattice from the specified directory"""
    size_dir = os.path.join(latest_dir, "nx8_ny8")
    if not os.path.exists(size_dir):
        print(f"No 8x8 lattice data found in {latest_dir}")
        return None
    
    data = []
    for layer in range(1, 4):
        layer_file = os.path.join(size_dir, f"layer_{layer}.txt")
        if os.path.exists(layer_file):
            layer_data = np.loadtxt(layer_file, comments='#')
            data.append(layer_data)
        else:
            print(f"Layer {layer} data not found for 8x8 lattice")
            return None
    
    print(f"Successfully loaded 8x8 data from {size_dir}")
    return data

def plot_thermodynamic_properties(data):
    """Plot thermodynamic properties for 8x8 lattice"""
    # Define properties and their corresponding column indices
    properties = [
        ("Magnetization²", 1),
        ("Magnetization⁴", 2),
        ("Energy", 5),
        ("Specific Heat", 6),
        ("Correlation Length", 7),
        ("Binder Cumulant", None)  # Special case, calculated from M² and M⁴
    ]
    
    # Create a figure with subplots
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3)
    
    # Theoretical critical temperature for 8-state Potts model
    tc_theory = 1.0 / np.log(1 + np.sqrt(8))  # ~0.751
    
    # Get temperature values (same for all layers)
    temps = data[0][:, 0]
    
    for i, (prop, col) in enumerate(properties):
        ax = fig.add_subplot(gs[i//3, i%3])
        
        if prop == "Binder Cumulant":
            # Calculate Binder cumulant: U₄ = 1 - ⟨m⁴⟩/(3⟨m²⟩²)
            for layer in range(3):
                m2 = data[layer][:, 1]  # M²
                m4 = data[layer][:, 2]  # M⁴
                binder = 1 - m4 / (3 * m2 * m2)
                ax.plot(temps, binder, 'o-', label=f"Layer {layer+1}")
        else:
            # Plot regular property for each layer
            for layer in range(3):
                ax.plot(temps, data[layer][:, col], 'o-', label=f"Layer {layer+1}")
        
        # Mark theoretical critical temperature
        ax.axvline(x=tc_theory, color='r', linestyle='--', alpha=0.7, 
                  label=f"Tc (theory) = {tc_theory:.3f}")
        
        # Add labels and legend
        ax.set_xlabel("Temperature", fontsize=12)
        ax.set_ylabel(prop, fontsize=12)
        ax.set_title(f"{prop} vs Temperature (8x8 lattice)", fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Add a title for the entire figure
    plt.suptitle("Thermodynamic Properties for 8x8 Lattice", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust for the suptitle
    
    # Save the figure
    plt.savefig("8x8_thermodynamic_properties.png", dpi=300)
    print("Plot saved as '8x8_thermodynamic_properties.png'")
    
    # Also create individual plots for each property
    for prop, col in properties:
        plt.figure(figsize=(10, 6))
        
        if prop == "Binder Cumulant":
            # Calculate Binder cumulant
            for layer in range(3):
                m2 = data[layer][:, 1]
                m4 = data[layer][:, 2]
                binder = 1 - m4 / (3 * m2 * m2)
                plt.plot(temps, binder, 'o-', label=f"Layer {layer+1}")
                
            # Also plot average over all layers
            m2_avg = np.mean([data[layer][:, 1] for layer in range(3)], axis=0)
            m4_avg = np.mean([data[layer][:, 2] for layer in range(3)], axis=0)
            binder_avg = 1 - m4_avg / (3 * m2_avg * m2_avg)
            plt.plot(temps, binder_avg, 'k-', linewidth=2, label="Average")
        else:
            # Plot regular property for each layer
            for layer in range(3):
                plt.plot(temps, data[layer][:, col], 'o-', label=f"Layer {layer+1}")
                
            # Also plot average over all layers
            avg_data = np.mean([data[layer][:, col] for layer in range(3)], axis=0)
            plt.plot(temps, avg_data, 'k-', linewidth=2, label="Average")
        
        # Mark theoretical critical temperature
        plt.axvline(x=tc_theory, color='r', linestyle='--', alpha=0.7, 
                   label=f"Tc (theory) = {tc_theory:.3f}")
        
        # Add labels and legend
        plt.xlabel("Temperature", fontsize=12)
        plt.ylabel(prop, fontsize=12)
        plt.title(f"{prop} vs Temperature (8x8 lattice)", fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save the individual plot
        filename = f"8x8_{prop.lower().replace('²', '2').replace('⁴', '4').replace(' ', '_')}.png"
        plt.savefig(filename, dpi=300)
        print(f"Individual plot saved as '{filename}'")
        plt.close()

def main():
    # Find the latest simulation directory
    latest_dir = find_latest_simulation()
    if not latest_dir:
        return
    
    print(f"Found latest simulation directory: {latest_dir}")
    
    # Load 8x8 data
    data = load_8x8_data(latest_dir)
    if not data:
        return
    
    # Plot thermodynamic properties
    plot_thermodynamic_properties(data)

if __name__ == "__main__":
    main() 