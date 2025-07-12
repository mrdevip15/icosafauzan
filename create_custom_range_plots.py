#!/usr/bin/env python3
"""
Create custom temperature range plots using all available simulation data
Focuses on temperature range 0.1 to 2.0 with all lattice sizes
"""

import os
import glob
from analyze_thermodynamics import ThermodynamicAnalyzer

def main():
    """Create custom temperature range plots using all simulation data"""
    print("Creating Custom Temperature Range Plots (0.1 to 2.0)")
    print("="*60)
    
    # Find all averaged data files
    data_files = {}
    for file in glob.glob("simulation_runs/*/smooth_averaged_data.txt"):
        # Extract lattice size from path (handle both Windows and Unix paths)
        path_parts = file.replace('\\', '/').split('/')
        if len(path_parts) >= 2:
            size_part = path_parts[-2]  # e.g., "nx8_ny8"
            if size_part.startswith('nx') and '_ny' in size_part:
                nx = int(size_part[2:size_part.find('_ny')])
                ny = int(size_part[size_part.find('ny')+2:])
                data_files[(nx, ny)] = file
    
    print(f"Found {len(data_files)} lattice sizes:")
    for (nx, ny), file in data_files.items():
        print(f"  {nx}×{ny}: {file}")
    
    if len(data_files) == 0:
        print("No averaged data files found!")
        return
    
    # Create analyzer instance
    analyzer = ThermodynamicAnalyzer()
    
    # Create custom temperature range plots
    print("\n" + "="*60)
    print("CREATING CUSTOM TEMPERATURE RANGE PLOTS")
    print("="*60)
    
    # 1. Custom temperature range (0.1 to 2.0)
    print("\n1. Creating custom temperature range plots (0.1 to 2.0)...")
    analyzer.plot_custom_temperature_range(
        data_files, 
        temp_min=0.1, 
        temp_max=2.0, 
        save_path='simulation_runs/custom_temperature_range_0.1_to_2.0.png'
    )
    
    # 2. Focused critical region (0.8 to 1.2)
    print("\n2. Creating focused critical region plots (0.8 to 1.2)...")
    analyzer.plot_focused_critical_region(
        data_files, 
        temp_min=0.8, 
        temp_max=1.2, 
        save_path='simulation_runs/focused_critical_region_0.8_to_1.2.png'
    )
    
    # 3. Extended critical region (0.5 to 1.5)
    print("\n3. Creating extended critical region plots (0.5 to 1.5)...")
    analyzer.plot_focused_critical_region(
        data_files, 
        temp_min=0.5, 
        temp_max=1.5, 
        save_path='simulation_runs/extended_critical_region_0.5_to_1.5.png'
    )
    
    # 4. High temperature region (1.0 to 2.0)
    print("\n4. Creating high temperature region plots (1.0 to 2.0)...")
    analyzer.plot_custom_temperature_range(
        data_files, 
        temp_min=1.0, 
        temp_max=2.0, 
        save_path='simulation_runs/high_temperature_region_1.0_to_2.0.png'
    )
    
    # 5. Low temperature region (0.1 to 1.0)
    print("\n5. Creating low temperature region plots (0.1 to 1.0)...")
    analyzer.plot_custom_temperature_range(
        data_files, 
        temp_min=0.1, 
        temp_max=1.0, 
        save_path='simulation_runs/low_temperature_region_0.1_to_1.0.png'
    )
    
    print("\n" + "="*60)
    print("CUSTOM TEMPERATURE RANGE PLOTS COMPLETE")
    print("="*60)
    print("Generated files:")
    print("  - simulation_runs/custom_temperature_range_0.1_to_2.0.png")
    print("  - simulation_runs/focused_critical_region_0.8_to_1.2.png")
    print("  - simulation_runs/extended_critical_region_0.5_to_1.5.png")
    print("  - simulation_runs/high_temperature_region_1.0_to_2.0.png")
    print("  - simulation_runs/low_temperature_region_0.1_to_1.0.png")
    print("\nAll plots use temperature range 0.1 to 2.0 with all available lattice sizes!")

if __name__ == "__main__":
    main() 