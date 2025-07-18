import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.gridspec as gridspec

# Parameters
LATTICE_SIZES = [8, 16, 24, 32]
THEORETICAL_TC = 1.0 / np.log(1 + np.sqrt(8))  # ~0.751 for 8-state Potts model
THEORETICAL_BETA = 1/8   # Critical exponent for magnetization
THEORETICAL_NU = 1       # Critical exponent for correlation length
THEORETICAL_ALPHA = 0    # Critical exponent for specific heat
THEORETICAL_GAMMA = 7/4  # Critical exponent for susceptibility

def load_simulation_data(nx, ny):
    """Load simulation data for all three layers"""
    data = []
    
    # Find the most recent simulation directory for this lattice size
    base_dir = "simulation_runs"
    target_dir = None
    
    # Check if simulation_runs directory exists
    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} does not exist yet. No simulation data found.")
        return None
    
    for dirname in sorted(os.listdir(base_dir), reverse=True):
        size_dir = os.path.join(base_dir, dirname, f"nx{nx}_ny{ny}")
        if os.path.exists(size_dir):
            target_dir = size_dir
            break
    
    if not target_dir:
        print(f"No simulation data found for lattice size {nx}x{ny}")
        return None
    
    # Load data from each layer
    for layer in range(1, 4):
        layer_file = os.path.join(target_dir, f"layer_{layer}.txt")
        if os.path.exists(layer_file):
            layer_data = np.loadtxt(layer_file, comments='#')
            data.append(layer_data)
        else:
            print(f"Layer {layer} data not found for lattice size {nx}x{ny}")
            return None
    
    print(f"Successfully loaded data for lattice size {nx}x{ny} from {target_dir}")
    return data

def plot_thermodynamic_properties(all_data):
    """Plot thermodynamic properties for all lattice sizes"""
    properties = ["Magnetization²", "Specific Heat", "Energy", "Correlation Length"]
    columns = [1, 6, 5, 7]  # Columns in data files
    
    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(2, 2)
    
    for i, (prop, col) in enumerate(zip(properties, columns)):
        ax = fig.add_subplot(gs[i//2, i%2])
        
        for size, data in all_data.items():
            # Average over all three layers
            avg_data = np.mean([layer[:, col] for layer in data], axis=0)
            temps = data[0][:, 0]  # Temperature column
            
            ax.plot(temps, avg_data, 'o-', label=f"L = {size}")
            
        ax.set_xlabel("Temperature", fontsize=14)
        ax.set_ylabel(prop, fontsize=14)
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Mark theoretical critical temperature
        ax.axvline(x=THEORETICAL_TC, color='r', linestyle='--', alpha=0.5, 
                  label=f"Theoretical Tc = {THEORETICAL_TC:.3f}")
    
    plt.tight_layout()
    plt.savefig("thermodynamic_properties_all_sizes.png", dpi=300)
    plt.close()
    
    print("Thermodynamic properties plot saved as 'thermodynamic_properties_all_sizes.png'")

def find_critical_temperature(all_data):
    """Find critical temperature using different methods"""
    methods = {
        "Specific Heat Peak": (6, find_peak),
        "Susceptibility Peak": (1, find_susceptibility_peak),
        "Binder Cumulant Crossing": (1, find_binder_crossing),
    }
    
    results = {}
    
    for method_name, (col_idx, method_func) in methods.items():
        tc_values = method_func(all_data, col_idx)
        results[method_name] = tc_values
    
    return results

def find_peak(all_data, col_idx):
    """Find temperature at which a property peaks for each lattice size"""
    tc_values = {}
    
    for size, data in all_data.items():
        # Average over all three layers
        avg_data = np.mean([layer[:, col_idx] for layer in data], axis=0)
        temps = data[0][:, 0]  # Temperature column
        
        # Find peak
        peak_idx = np.argmax(avg_data)
        tc_values[size] = temps[peak_idx]
    
    return tc_values

def find_susceptibility_peak(all_data, col_idx):
    """Find temperature at which susceptibility peaks for each lattice size"""
    tc_values = {}
    
    for size, data in all_data.items():
        # Average over all three layers
        avg_m2 = np.mean([layer[:, col_idx] for layer in data], axis=0)
        temps = data[0][:, 0]  # Temperature column
        
        # Calculate susceptibility: χ = N⋅(⟨m²⟩)/T
        susceptibility = size*size * avg_m2 / temps
        
        # Find peak
        peak_idx = np.argmax(susceptibility)
        tc_values[size] = temps[peak_idx]
    
    return tc_values

def find_binder_crossing(all_data, col_idx):
    """Find temperature at which Binder cumulant curves cross"""
    binder_data = {}
    temps = None
    
    for size, data in all_data.items():
        # Calculate Binder cumulant: U₄ = 1 - ⟨m⁴⟩/(3⟨m²⟩²)
        avg_m2 = np.mean([layer[:, 1] for layer in data], axis=0)  # m²
        avg_m4 = np.mean([layer[:, 2] for layer in data], axis=0)  # m⁴
        temps = data[0][:, 0]  # Temperature column
        
        binder = 1 - avg_m4 / (3 * avg_m2 * avg_m2)
        binder_data[size] = binder
    
    # Find crossing points between all pairs of lattice sizes
    crossings = []
    sizes = list(binder_data.keys())
    
    for i in range(len(sizes)):
        for j in range(i+1, len(sizes)):
            size1, size2 = sizes[i], sizes[j]
            b1, b2 = binder_data[size1], binder_data[size2]
            
            # Find where curves cross by finding where their difference changes sign
            diff = b1 - b2
            for k in range(len(diff)-1):
                if diff[k] * diff[k+1] <= 0:  # Sign change
                    # Linear interpolation to find crossing point
                    t1, t2 = temps[k], temps[k+1]
                    d1, d2 = diff[k], diff[k+1]
                    tc = t1 - d1 * (t2 - t1) / (d2 - d1)
                    crossings.append(tc)
    
    # Return average of crossing points
    if crossings:
        avg_tc = np.mean(crossings)
        return {size: avg_tc for size in sizes}
    else:
        return {size: np.nan for size in sizes}

def finite_size_scaling(critical_temps):
    """Perform finite size scaling analysis to extrapolate to infinite system"""
    plt.figure(figsize=(10, 8))
    
    for method, tc_values in critical_temps.items():
        sizes = np.array(list(tc_values.keys()))
        tcs = np.array(list(tc_values.values()))
        
        # Fit to Tc(L) = Tc(∞) + a/L^(1/ν)
        # For 8-state Potts model, ν = 1
        def scaling_func(x, tc_inf, a):
            return tc_inf + a / x
        
        try:
            params, _ = curve_fit(scaling_func, sizes, tcs)
            tc_inf, a = params
            
            # Plot data and fit
            plt.scatter(1/sizes, tcs, label=f"{method} data")
            
            x_fit = np.linspace(1/max(sizes), 1/min(sizes), 100)
            plt.plot(x_fit, scaling_func(1/x_fit, tc_inf, a), '--', 
                    label=f"{method}: Tc(∞) = {tc_inf:.4f}")
        except:
            print(f"Fitting failed for {method}")
    
    # Add theoretical value
    plt.axhline(y=THEORETICAL_TC, color='r', linestyle='-', 
               label=f"Theoretical Tc = {THEORETICAL_TC:.4f}")
    
    plt.xlabel("1/L", fontsize=14)
    plt.ylabel("Critical Temperature Tc(L)", fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.title("Finite Size Scaling of Critical Temperature", fontsize=16)
    plt.savefig("finite_size_scaling.png", dpi=300)
    plt.close()
    
    print("Finite size scaling plot saved as 'finite_size_scaling.png'")

def critical_exponents_analysis(all_data):
    """Analyze critical exponents using finite size scaling collapse"""
    # Extract magnetization and specific heat near critical point
    tc = THEORETICAL_TC
    window = 0.1  # Temperature window around Tc
    
    m_data = {}  # Magnetization data
    c_data = {}  # Specific heat data
    chi_data = {}  # Susceptibility data
    
    for size, data in all_data.items():
        temps = data[0][:, 0]
        
        # Get indices within temperature window
        idx = np.where((temps >= tc - window) & (temps <= tc + window))[0]
        
        # Extract data
        t_vals = temps[idx]
        t_scaled = (t_vals - tc) / tc  # (T-Tc)/Tc
        
        # Magnetization (column 1)
        m_vals = np.mean([layer[idx, 1] for layer in data], axis=0)
        m_data[size] = (t_scaled, np.sqrt(m_vals))  # m = sqrt(m²)
        
        # Specific heat (column 6)
        c_vals = np.mean([layer[idx, 6] for layer in data], axis=0)
        c_data[size] = (t_scaled, c_vals)
        
        # Susceptibility: χ = N⋅(⟨m²⟩)/T
        chi_vals = size*size * m_vals / t_vals
        chi_data[size] = (t_scaled, chi_vals)
    
    # Plot data collapses
    fig = plt.figure(figsize=(18, 6))
    
    # 1. Magnetization: M * L^(β/ν) vs (T-Tc)/Tc * L^(1/ν)
    ax1 = fig.add_subplot(131)
    for size, (t, m) in m_data.items():
        scaled_m = m * size**(THEORETICAL_BETA/THEORETICAL_NU)
        scaled_t = t * size**(1/THEORETICAL_NU)
        ax1.plot(scaled_t, scaled_m, 'o-', label=f"L = {size}")
    
    ax1.set_xlabel(r"$(T-T_c)/T_c \cdot L^{1/\nu}$", fontsize=12)
    ax1.set_ylabel(r"$M \cdot L^{\beta/\nu}$", fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f"Magnetization Scaling\nβ = {THEORETICAL_BETA}, ν = {THEORETICAL_NU}", fontsize=12)
    
    # 2. Specific Heat: C * L^(-α/ν) vs (T-Tc)/Tc * L^(1/ν)
    ax2 = fig.add_subplot(132)
    for size, (t, c) in c_data.items():
        scaled_c = c * size**(-THEORETICAL_ALPHA/THEORETICAL_NU)
        scaled_t = t * size**(1/THEORETICAL_NU)
        ax2.plot(scaled_t, scaled_c, 'o-', label=f"L = {size}")
    
    ax2.set_xlabel(r"$(T-T_c)/T_c \cdot L^{1/\nu}$", fontsize=12)
    ax2.set_ylabel(r"$C \cdot L^{-\alpha/\nu}$", fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_title(f"Specific Heat Scaling\nα = {THEORETICAL_ALPHA}, ν = {THEORETICAL_NU}", fontsize=12)
    
    # 3. Susceptibility: χ * L^(-γ/ν) vs (T-Tc)/Tc * L^(1/ν)
    ax3 = fig.add_subplot(133)
    for size, (t, chi) in chi_data.items():
        scaled_chi = chi * size**(-THEORETICAL_GAMMA/THEORETICAL_NU)
        scaled_t = t * size**(1/THEORETICAL_NU)
        ax3.plot(scaled_t, scaled_chi, 'o-', label=f"L = {size}")
    
    ax3.set_xlabel(r"$(T-T_c)/T_c \cdot L^{1/\nu}$", fontsize=12)
    ax3.set_ylabel(r"$\chi \cdot L^{-\gamma/\nu}$", fontsize=12)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_title(f"Susceptibility Scaling\nγ = {THEORETICAL_GAMMA}, ν = {THEORETICAL_NU}", fontsize=12)
    
    plt.tight_layout()
    plt.savefig("critical_exponents_analysis.png", dpi=300)
    plt.close()
    
    print("Critical exponents analysis plot saved as 'critical_exponents_analysis.png'")

def plot_binder_cumulant(all_data):
    """Plot Binder cumulant for different lattice sizes"""
    plt.figure(figsize=(10, 8))
    
    for size, data in all_data.items():
        # Calculate Binder cumulant: U₄ = 1 - ⟨m⁴⟩/(3⟨m²⟩²)
        avg_m2 = np.mean([layer[:, 1] for layer in data], axis=0)  # m²
        avg_m4 = np.mean([layer[:, 2] for layer in data], axis=0)  # m⁴
        temps = data[0][:, 0]  # Temperature column
        
        binder = 1 - avg_m4 / (3 * avg_m2 * avg_m2)
        plt.plot(temps, binder, 'o-', label=f"L = {size}")
    
    # Add theoretical critical temperature
    plt.axvline(x=THEORETICAL_TC, color='r', linestyle='--', alpha=0.5, 
               label=f"Theoretical Tc = {THEORETICAL_TC:.3f}")
    
    plt.xlabel("Temperature", fontsize=14)
    plt.ylabel("Binder Cumulant U₄", fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.title("Binder Cumulant Analysis", fontsize=16)
    plt.savefig("binder_cumulant_analysis.png", dpi=300)
    plt.close()
    
    print("Binder cumulant plot saved as 'binder_cumulant_analysis.png'")

def compare_with_theory(critical_temps):
    """Compare measured critical temperatures with theoretical prediction"""
    # Extract infinite-size extrapolations
    methods = list(critical_temps.keys())
    tc_values = {}
    
    for method in methods:
        sizes = np.array(list(critical_temps[method].keys()))
        tcs = np.array(list(critical_temps[method].values()))
        
        # Fit to Tc(L) = Tc(∞) + a/L^(1/ν)
        def scaling_func(x, tc_inf, a):
            return tc_inf + a / x
        
        try:
            params, _ = curve_fit(scaling_func, sizes, tcs)
            tc_inf, a = params
            tc_values[method] = tc_inf
        except:
            tc_values[method] = np.nan
    
    # Create comparison table
    comparison = pd.DataFrame({
        'Method': list(tc_values.keys()) + ['Theoretical'],
        'Critical Temperature': list(tc_values.values()) + [THEORETICAL_TC],
        'Difference from Theory': [abs(tc - THEORETICAL_TC) for tc in tc_values.values()] + [0],
        'Relative Error (%)': [abs(tc - THEORETICAL_TC) / THEORETICAL_TC * 100 for tc in tc_values.values()] + [0]
    })
    
    # Save comparison to CSV
    comparison.to_csv("critical_temperature_comparison.csv", index=False)
    
    # Create bar chart
    plt.figure(figsize=(12, 8))
    
    methods = list(tc_values.keys()) + ['Theoretical']
    temps = list(tc_values.values()) + [THEORETICAL_TC]
    
    bars = plt.bar(methods, temps, alpha=0.7)
    
    # Add horizontal line for theoretical value
    plt.axhline(y=THEORETICAL_TC, color='r', linestyle='--', alpha=0.5, 
               label=f"Theoretical Tc = {THEORETICAL_TC:.4f}")
    
    plt.xlabel("Method", fontsize=14)
    plt.ylabel("Critical Temperature", fontsize=14)
    plt.title("Critical Temperature Comparison", fontsize=16)
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.005,
                f'{height:.4f}', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig("critical_temperature_methods_comparison.png", dpi=300)
    plt.close()
    
    print("Critical temperature comparison saved as 'critical_temperature_methods_comparison.png'")
    print("Detailed comparison data saved as 'critical_temperature_comparison.csv'")
    
    # Generate summary report
    with open("CRITICAL_TEMPERATURE_SUMMARY.md", "w") as f:
        f.write("# Critical Temperature Analysis Summary\n\n")
        f.write(f"Theoretical critical temperature for 8-state Potts model: Tc = {THEORETICAL_TC:.6f}\n\n")
        f.write("## Measured Critical Temperatures\n\n")
        f.write("| Method | Critical Temperature | Difference from Theory | Relative Error (%) |\n")
        f.write("|--------|----------------------|------------------------|--------------------|\n")
        
        for method, tc in tc_values.items():
            diff = abs(tc - THEORETICAL_TC)
            rel_err = diff / THEORETICAL_TC * 100
            f.write(f"| {method} | {tc:.6f} | {diff:.6f} | {rel_err:.2f} |\n")
    
    print("Critical temperature summary report saved as 'CRITICAL_TEMPERATURE_SUMMARY.md'")

def main():
    """Main function to analyze simulation results"""
    print("Starting analysis of simulation results...")
    
    # Load simulation data for all lattice sizes
    all_data = {}
    available_sizes = []
    
    for size in LATTICE_SIZES:
        data = load_simulation_data(size, size)
        if data is not None:
            all_data[size] = data
            available_sizes.append(size)
    
    if not all_data:
        print("No simulation data available. Exiting.")
        return
    
    print(f"Found data for lattice sizes: {available_sizes}")
    
    # Plot thermodynamic properties
    plot_thermodynamic_properties(all_data)
    
    # Find critical temperatures using different methods
    critical_temps = find_critical_temperature(all_data)
    
    # Perform finite size scaling analysis
    finite_size_scaling(critical_temps)
    
    # Analyze critical exponents
    critical_exponents_analysis(all_data)
    
    # Plot Binder cumulant
    plot_binder_cumulant(all_data)
    
    # Compare with theory
    compare_with_theory(critical_temps)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main() 