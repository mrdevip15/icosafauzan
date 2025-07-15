#!/usr/bin/env python3
"""
Detailed critical temperature analysis using specific heat and other thermodynamic indicators.
Focuses on precise determination of critical points and critical exponents.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.signal import find_peaks, savgol_filter
import warnings
warnings.filterwarnings('ignore')

class CriticalTemperatureAnalyzer:
    def __init__(self, data_dir="simulation_runs"):
        self.data_dir = data_dir
        self.layer_data = {}
        self.lattice_sizes = []
        self.layer_colors = ['#e41a1c', '#377eb8', '#4daf4a']  # Red, Blue, Green
        
    def load_data(self):
        """Load all simulation data."""
        print("Loading simulation data for critical temperature analysis...")
        
        for size_dir in glob.glob(f"{self.data_dir}/nx*_ny*"):
            size_name = os.path.basename(size_dir)
            parts = size_name.split('_')
            nx = int(parts[0][2:])
            ny = int(parts[1][2:])
            
            if nx == ny:
                self.lattice_sizes.append(nx)
                self.layer_data[nx] = {}
                
                for layer in [1, 2, 3]:
                    layer_file = f"{size_dir}/layer_{layer}.txt"
                    if os.path.exists(layer_file):
                        data = pd.read_csv(layer_file, sep=r'\s+', comment='#',
                                         names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                        self.layer_data[nx][layer] = data
        
        self.lattice_sizes.sort()
        print(f"Loaded data for lattice sizes: {self.lattice_sizes}")
    
    def calculate_susceptibility(self, data, size):
        """Calculate magnetic susceptibility from magnetization fluctuations."""
        # χ = β * L^d * (<M²> - <|M|>²)
        # For our data, we have <M²> directly
        # Susceptibility ∝ <M²> * β * V where V = L²d (d=2 effective)
        beta = 1.0 / data['Temperature']
        volume = size * size  # Per layer volume
        susceptibility = beta * volume * data['M2']
        return susceptibility
    
    def find_critical_temperature_methods(self, size, layer):
        """Find critical temperature using multiple methods."""
        if size not in self.layer_data or layer not in self.layer_data[size]:
            return {}
        
        data = self.layer_data[size][layer]
        
        # Focus on reasonable temperature range
        mask = (data['Temperature'] >= 0.5) & (data['Temperature'] <= 2.5)
        temp = data['Temperature'][mask].values
        cv = data['Cv'][mask].values
        m2 = data['M2'][mask].values
        m4 = data['M4'][mask].values
        energy = data['Energy'][mask].values
        
        results = {}
        
        # Method 1: Peak of specific heat
        try:
            # Filter out computational artifacts (very high Cv values)
            reasonable_cv = cv < 1000  # Adjust threshold as needed
            if reasonable_cv.sum() > 10:
                temp_filt = temp[reasonable_cv]
                cv_filt = cv[reasonable_cv]
                
                # Smooth the data to find reliable peak
                if len(temp_filt) > 10:
                    cv_smooth = savgol_filter(cv_filt, window_length=min(11, len(cv_filt)//2*2+1), polyorder=3)
                    peaks, _ = find_peaks(cv_smooth, height=np.max(cv_smooth)*0.5)
                    
                    if len(peaks) > 0:
                        peak_idx = peaks[np.argmax(cv_smooth[peaks])]
                        tc_cv = temp_filt[peak_idx]
                        results['cv_peak'] = tc_cv
                        results['cv_max'] = cv_smooth[peak_idx]
        except:
            pass
        
        # Method 2: Inflection point of magnetization
        try:
            if len(temp) > 10:
                # Calculate derivative of M²
                dm2_dt = np.gradient(m2, temp)
                
                # Find minimum of derivative (steepest drop)
                min_idx = np.argmin(dm2_dt)
                if 0 < min_idx < len(temp) - 1:
                    tc_mag = temp[min_idx]
                    results['magnetization_inflection'] = tc_mag
        except:
            pass
        
        # Method 3: Peak of susceptibility
        try:
            susceptibility = self.calculate_susceptibility(data[mask], size)
            chi = susceptibility.values
            
            if len(chi) > 10 and np.max(chi) > 0:
                chi_smooth = savgol_filter(chi, window_length=min(11, len(chi)//2*2+1), polyorder=3)
                peaks, _ = find_peaks(chi_smooth, height=np.max(chi_smooth)*0.3)
                
                if len(peaks) > 0:
                    peak_idx = peaks[np.argmax(chi_smooth[peaks])]
                    tc_chi = temp[peak_idx]
                    results['susceptibility_peak'] = tc_chi
        except:
            pass
        
        # Method 4: Binder cumulant crossing (for comparison with universal value)
        try:
            valid_mask = (m2 > 1e-10) & (m4 > 1e-10)
            if valid_mask.sum() > 10:
                temp_valid = temp[valid_mask]
                m2_valid = m2[valid_mask]
                m4_valid = m4[valid_mask]
                
                binder = 1 - m4_valid / (3 * m2_valid**2)
                
                # Filter reasonable Binder values
                binder_mask = (binder > 0.3) & (binder < 0.9) & np.isfinite(binder)
                if binder_mask.sum() > 5:
                    temp_binder = temp_valid[binder_mask]
                    binder_filt = binder[binder_mask]
                    
                    # Find crossing with universal value 0.61
                    target_binder = 0.61
                    crossing_idx = np.argmin(np.abs(binder_filt - target_binder))
                    tc_binder = temp_binder[crossing_idx]
                    results['binder_crossing'] = tc_binder
        except:
            pass
        
        # Method 5: Peak of energy derivative (dE/dT ~ specific heat)
        try:
            if len(temp) > 10:
                de_dt = np.gradient(energy, temp)
                de_dt_smooth = savgol_filter(de_dt, window_length=min(11, len(de_dt)//2*2+1), polyorder=3)
                
                peaks, _ = find_peaks(-de_dt_smooth)  # Find minima (most negative derivative)
                if len(peaks) > 0:
                    peak_idx = peaks[np.argmin(de_dt_smooth[peaks])]
                    tc_energy = temp[peak_idx]
                    results['energy_derivative_peak'] = tc_energy
        except:
            pass
        
        return results
    
    def analyze_critical_exponents(self, size, layer):
        """Analyze critical exponents near the critical temperature."""
        if size not in self.layer_data or layer not in self.layer_data[size]:
            return {}
        
        data = self.layer_data[size][layer]
        
        # Get critical temperature estimate
        tc_methods = self.find_critical_temperature_methods(size, layer)
        if not tc_methods:
            return {}
        
        # Use specific heat peak as primary estimate
        tc = tc_methods.get('cv_peak', tc_methods.get('magnetization_inflection', 1.5))
        
        results = {'tc_used': tc}
        
        # Analyze specific heat exponent α
        try:
            # Cv ~ |t|^(-α) where t = (T-Tc)/Tc
            temp_range = (data['Temperature'] >= tc - 0.3) & (data['Temperature'] <= tc + 0.3)
            temp_fit = data['Temperature'][temp_range].values
            cv_fit = data['Cv'][temp_range].values
            
            # Filter reasonable values
            reasonable = (cv_fit > 0) & (cv_fit < 1000) & (temp_fit != tc)
            if reasonable.sum() > 10:
                temp_good = temp_fit[reasonable]
                cv_good = cv_fit[reasonable]
                
                # Reduced temperature
                t = np.abs((temp_good - tc) / tc)
                
                # Fit Cv ~ t^(-α)
                log_t = np.log(t + 1e-10)
                log_cv = np.log(cv_good + 1e-10)
                
                # Linear fit in log space
                valid_log = np.isfinite(log_t) & np.isfinite(log_cv)
                if valid_log.sum() > 5:
                    coeffs = np.polyfit(log_t[valid_log], log_cv[valid_log], 1)
                    alpha = -coeffs[0]  # Negative because Cv ~ t^(-α)
                    results['alpha_exponent'] = alpha
        except:
            pass
        
        # Analyze magnetization exponent β
        try:
            # Below Tc: M ~ |t|^β where t = (Tc-T)/Tc
            below_tc = (data['Temperature'] < tc) & (data['Temperature'] > tc - 0.3)
            temp_below = data['Temperature'][below_tc].values
            m2_below = data['M2'][below_tc].values
            
            if len(temp_below) > 5:
                t_below = (tc - temp_below) / tc
                m_below = np.sqrt(m2_below)  # |M| from <M²>
                
                # Fit M ~ t^β
                log_t = np.log(t_below + 1e-10)
                log_m = np.log(m_below + 1e-10)
                
                valid_log = np.isfinite(log_t) & np.isfinite(log_m) & (m_below > 0.1)
                if valid_log.sum() > 3:
                    coeffs = np.polyfit(log_t[valid_log], log_m[valid_log], 1)
                    beta = coeffs[0]
                    results['beta_exponent'] = beta
        except:
            pass
        
        # Analyze susceptibility exponent γ
        try:
            # χ ~ |t|^(-γ)
            susceptibility = self.calculate_susceptibility(data, size)
            temp_range = (data['Temperature'] >= tc - 0.3) & (data['Temperature'] <= tc + 0.3)
            temp_fit = data['Temperature'][temp_range].values
            chi_fit = susceptibility[temp_range].values
            
            reasonable = (chi_fit > 0) & np.isfinite(chi_fit) & (temp_fit != tc)
            if reasonable.sum() > 10:
                temp_good = temp_fit[reasonable]
                chi_good = chi_fit[reasonable]
                
                t = np.abs((temp_good - tc) / tc)
                
                log_t = np.log(t + 1e-10)
                log_chi = np.log(chi_good + 1e-10)
                
                valid_log = np.isfinite(log_t) & np.isfinite(log_chi)
                if valid_log.sum() > 5:
                    coeffs = np.polyfit(log_t[valid_log], log_chi[valid_log], 1)
                    gamma = -coeffs[0]
                    results['gamma_exponent'] = gamma
        except:
            pass
        
        return results
    
    def plot_critical_temperature_determination(self):
        """Plot detailed critical temperature determination for each method."""
        fig, axes = plt.subplots(3, 4, figsize=(20, 15))
        
        methods = ['cv_peak', 'magnetization_inflection', 'susceptibility_peak', 'binder_crossing']
        method_labels = ['Specific Heat Peak', 'Magnetization Inflection', 'Susceptibility Peak', 'Binder Crossing']
        
        all_tc_data = {}
        
        for layer in [1, 2, 3]:
            layer_tc_data = {}
            
            for size in self.lattice_sizes:
                tc_results = self.find_critical_temperature_methods(size, layer)
                if tc_results:
                    layer_tc_data[size] = tc_results
            
            all_tc_data[layer] = layer_tc_data
            
            # Plot each method
            for j, (method, label) in enumerate(zip(methods, method_labels)):
                ax = axes[layer-1, j]
                
                sizes = []
                tcs = []
                
                for size, tc_dict in layer_tc_data.items():
                    if method in tc_dict:
                        sizes.append(size)
                        tcs.append(tc_dict[method])
                
                if len(sizes) > 1:
                    sizes = np.array(sizes)
                    tcs = np.array(tcs)
                    
                    # Plot data points
                    ax.plot(sizes, tcs, 'o-', color=self.layer_colors[layer-1], 
                           linewidth=2, markersize=8, label=f'Layer {layer}')
                    
                    # Fit finite-size scaling: Tc(L) = Tc(∞) + a/L
                    if len(sizes) > 2:
                        try:
                            inv_sizes = 1.0 / sizes
                            coeffs = np.polyfit(inv_sizes, tcs, 1)
                            tc_inf = coeffs[1]  # Intercept
                            
                            # Plot extrapolation
                            x_fit = np.linspace(0, inv_sizes.max(), 100)
                            y_fit = coeffs[0] * x_fit + tc_inf
                            size_fit = 1.0 / (x_fit + 1e-10)
                            
                            # Plot finite-size extrapolation
                            ax.plot(size_fit[x_fit > 1e-3], y_fit[x_fit > 1e-3], '--', 
                                   color=self.layer_colors[layer-1], alpha=0.7)
                            
                            # Mark infinite size limit
                            ax.axhline(y=tc_inf, color=self.layer_colors[layer-1], 
                                     linestyle=':', alpha=0.5, 
                                     label=f'Tc(∞) = {tc_inf:.3f}')
                            
                        except:
                            pass
                
                ax.set_xlabel('Lattice Size L')
                ax.set_ylabel('Critical Temperature')
                ax.set_title(f'Layer {layer}: {label}')
                ax.grid(True, alpha=0.3)
                ax.legend()
        
        plt.tight_layout()
        plt.savefig('critical_temperature_methods_comparison.png', dpi=300, bbox_inches='tight')
        print("Saved: critical_temperature_methods_comparison.png")
        plt.show()
        
        return all_tc_data
    
    def plot_critical_exponents_analysis(self):
        """Plot critical exponents analysis."""
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        
        exponents = ['alpha_exponent', 'beta_exponent', 'gamma_exponent']
        exponent_labels = ['Specific Heat α', 'Magnetization β', 'Susceptibility γ']
        theoretical_values = [0.0, 0.125, 1.75]  # 2D Ising values
        
        print("\n=== Critical Exponents Analysis ===")
        
        for layer in [1, 2, 3]:
            print(f"\nLayer {layer}:")
            
            for j, (exp, label, theory) in enumerate(zip(exponents, exponent_labels, theoretical_values)):
                ax = axes[layer-1, j]
                
                sizes = []
                exps = []
                
                for size in self.lattice_sizes:
                    exp_results = self.analyze_critical_exponents(size, layer)
                    if exp in exp_results:
                        sizes.append(size)
                        exps.append(exp_results[exp])
                        print(f"  L={size}: {label} = {exp_results[exp]:.3f}")
                
                if len(sizes) > 0:
                    ax.plot(sizes, exps, 'o-', color=self.layer_colors[layer-1], 
                           linewidth=2, markersize=8, label=f'Layer {layer}')
                    
                    # Show theoretical value
                    ax.axhline(y=theory, color='black', linestyle='--', alpha=0.5, 
                             label=f'2D Ising: {theory}')
                
                ax.set_xlabel('Lattice Size L')
                ax.set_ylabel(f'Exponent {exp.split("_")[0]}')
                ax.set_title(f'Layer {layer}: {label}')
                ax.grid(True, alpha=0.3)
                ax.legend()
        
        plt.tight_layout()
        plt.savefig('critical_exponents_analysis.png', dpi=300, bbox_inches='tight')
        print("\nSaved: critical_exponents_analysis.png")
        plt.show()
    
    def plot_detailed_critical_region(self):
        """Plot detailed behavior in the critical region."""
        # Use largest lattice size for detailed analysis
        if not self.lattice_sizes:
            return
            
        size = self.lattice_sizes[-1]
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        properties = [
            ('Cv', 'Specific Heat $C_v$'),
            ('M2', 'Magnetization $\\langle M^2 \\rangle$'),
            ('Energy', 'Energy per site'),
        ]
        
        for layer in [1, 2, 3]:
            if size not in self.layer_data or layer not in self.layer_data[size]:
                continue
                
            data = self.layer_data[size][layer]
            
            # Get critical temperature
            tc_methods = self.find_critical_temperature_methods(size, layer)
            tc = tc_methods.get('cv_peak', tc_methods.get('magnetization_inflection', 1.5))
            
            # Focus on critical region
            critical_range = (data['Temperature'] >= tc - 0.5) & (data['Temperature'] <= tc + 0.5)
            temp_crit = data['Temperature'][critical_range]
            
            for j, (prop, label) in enumerate(properties):
                ax = axes[0, j] if layer == 1 else (axes[1, j] if layer == 2 else axes[0, j])
                
                values_crit = data[prop][critical_range]
                
                # Filter reasonable values for specific heat
                if prop == 'Cv':
                    reasonable = values_crit < 2000
                    temp_plot = temp_crit[reasonable]
                    values_plot = values_crit[reasonable]
                else:
                    temp_plot = temp_crit
                    values_plot = values_crit
                
                if len(temp_plot) > 5:
                    ax.plot(temp_plot, values_plot, 'o-', 
                           color=self.layer_colors[layer-1], 
                           linewidth=2, markersize=4, 
                           alpha=0.8, label=f'Layer {layer}')
                    
                    # Mark critical temperature
                    if prop == 'Cv' and 'cv_peak' in tc_methods:
                        ax.axvline(x=tc, color=self.layer_colors[layer-1], 
                                 linestyle='--', alpha=0.7, 
                                 label=f'Tc = {tc:.3f}')
            
            # Add susceptibility plot
            ax_sus = axes[1, 0] if layer != 1 else axes[1, 1]
            
            susceptibility = self.calculate_susceptibility(data[critical_range], size)
            if len(susceptibility) > 5:
                ax_sus.plot(temp_crit, susceptibility, 'o-', 
                           color=self.layer_colors[layer-1], 
                           linewidth=2, markersize=4, 
                           alpha=0.8, label=f'Layer {layer}')
                ax_sus.set_xlabel('Temperature T')
                ax_sus.set_ylabel('Susceptibility χ')
                ax_sus.set_title('Magnetic Susceptibility')
                ax_sus.grid(True, alpha=0.3)
                ax_sus.legend()
        
        # Set labels for main plots
        for j, (prop, label) in enumerate(properties):
            axes[0, j].set_xlabel('Temperature T')
            axes[0, j].set_ylabel(label)
            axes[0, j].set_title(f'{label} - Critical Region (L={size})')
            axes[0, j].grid(True, alpha=0.3)
            axes[0, j].legend()
        
        plt.tight_layout()
        plt.savefig('detailed_critical_region_analysis.png', dpi=300, bbox_inches='tight')
        print("Saved: detailed_critical_region_analysis.png")
        plt.show()
    
    def generate_critical_temperature_summary(self):
        """Generate comprehensive summary of critical temperature analysis."""
        print("\n" + "="*60)
        print("CRITICAL TEMPERATURE ANALYSIS SUMMARY")
        print("="*60)
        
        all_results = {}
        
        for layer in [1, 2, 3]:
            print(f"\n--- LAYER {layer} ---")
            layer_results = {}
            
            for size in self.lattice_sizes:
                tc_methods = self.find_critical_temperature_methods(size, layer)
                exp_results = self.analyze_critical_exponents(size, layer)
                
                if tc_methods or exp_results:
                    combined_results = {**tc_methods, **exp_results}
                    layer_results[size] = combined_results
                    
                    print(f"\nLattice Size {size}x{size}:")
                    
                    # Critical temperatures
                    if tc_methods:
                        print("  Critical Temperature Methods:")
                        for method, tc in tc_methods.items():
                            print(f"    {method:25}: T_c = {tc:.4f}")
                    
                    # Critical exponents
                    if exp_results:
                        print("  Critical Exponents:")
                        for exp, value in exp_results.items():
                            if 'exponent' in exp:
                                print(f"    {exp:25}: {value:.4f}")
            
            all_results[layer] = layer_results
        
        # Finite-size extrapolation summary
        print(f"\n--- FINITE-SIZE EXTRAPOLATION ---")
        for layer in [1, 2, 3]:
            print(f"\nLayer {layer}:")
            
            methods = ['cv_peak', 'magnetization_inflection', 'susceptibility_peak']
            
            for method in methods:
                sizes = []
                tcs = []
                
                for size, results in all_results[layer].items():
                    if method in results:
                        sizes.append(size)
                        tcs.append(results[method])
                
                if len(sizes) > 2:
                    try:
                        sizes = np.array(sizes)
                        tcs = np.array(tcs)
                        
                        # Fit Tc(L) = Tc(∞) + a/L
                        inv_sizes = 1.0 / sizes
                        coeffs = np.polyfit(inv_sizes, tcs, 1)
                        tc_inf = coeffs[1]
                        
                        print(f"  {method:25}: T_c(∞) = {tc_inf:.4f}")
                        
                    except:
                        print(f"  {method:25}: Extrapolation failed")
        
        return all_results
    
    def run_complete_analysis(self):
        """Run complete critical temperature analysis."""
        print("=== COMPREHENSIVE CRITICAL TEMPERATURE ANALYSIS ===")
        print("Quasi-3D 8-state cube spins - Layer-separated Monte Carlo")
        print("="*60)
        
        self.load_data()
        
        if not self.layer_data:
            print("No data found!")
            return
        
        print(f"\nAnalyzing critical behavior for lattice sizes: {self.lattice_sizes}")
        print("Using multiple thermodynamic indicators...")
        
        # Generate all analyses
        tc_data = self.plot_critical_temperature_determination()
        self.plot_critical_exponents_analysis()
        self.plot_detailed_critical_region()
        
        # Generate comprehensive summary
        results = self.generate_critical_temperature_summary()
        
        print(f"\n{'='*60}")
        print("ANALYSIS COMPLETE")
        print("="*60)
        print("Generated files:")
        print("- critical_temperature_methods_comparison.png")
        print("- critical_exponents_analysis.png") 
        print("- detailed_critical_region_analysis.png")
        
        return results

def main():
    """Main analysis function."""
    analyzer = CriticalTemperatureAnalyzer()
    results = analyzer.run_complete_analysis()
    return results

if __name__ == "__main__":
    main()