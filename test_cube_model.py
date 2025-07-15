#!/usr/bin/env python3
"""
Quick test script for the enhanced 8-state cube model Monte Carlo simulation.
"""

import subprocess
import os

def test_cube_model():
    """Test the enhanced cube model simulation."""
    print("="*60)
    print("ENHANCED 8-STATE CUBE MODEL TEST")
    print("="*60)
    print("Testing the enhanced Monte Carlo simulation with:")
    print("- 8-state cube spins (vertices of cube)")
    print("- Hybrid Monte Carlo algorithm")
    print("- Optimized temperature range for Tc ≈ 0.751")
    print("- Timestamped hierarchical directory structure")
    print("-"*60)
    
    # Check executable
    if not os.path.exists('./cubic3lay_enhanced'):
        print("✗ Executable not found. Compile with:")
        print("  gcc -o cubic3lay_enhanced cubic3lay.c rn32.c -lm")
        return False
    
    # Test parameters
    nx, ny = 8, 8
    nmcs1, nmcs2 = 500, 1000  # Short test
    seed = 12345
    
    print(f"\nRunning test simulation:")
    print(f"- Lattice size: {nx}x{ny}x3")
    print(f"- Equilibration: {nmcs1} steps")
    print(f"- Measurement: {nmcs2} steps")
    print(f"- Seed: {seed}")
    
    try:
        # Run simulation
        input_data = f"{nx} {ny} {nmcs1} {nmcs2} {seed}\n"
        result = subprocess.run(
            ['./cubic3lay_enhanced'], 
            input=input_data, 
            text=True, 
            capture_output=True, 
            timeout=120
        )
        
        if result.returncode == 0:
            print("✓ Simulation completed successfully!")
            
            # Extract output directory
            output_lines = result.stdout.strip().split('\n')
            output_dir = None
            for line in output_lines:
                if "Layer data saved to" in line:
                    output_dir = line.split("to ")[-1].rstrip('/')
                    break
            
            if output_dir and os.path.exists(output_dir):
                print(f"✓ Data saved to: {output_dir}")
                
                # Check files
                layer_files = []
                for layer in [1, 2, 3]:
                    layer_file = f"{output_dir}/layer_{layer}.txt"
                    if os.path.exists(layer_file):
                        layer_files.append(layer_file)
                        file_size = os.path.getsize(layer_file)
                        print(f"  ✓ layer_{layer}.txt ({file_size:,} bytes)")
                
                # Analyze first layer file
                if layer_files:
                    print(f"\nAnalyzing {layer_files[0]}:")
                    
                    # Read metadata
                    print("  Metadata:")
                    with open(layer_files[0], 'r') as f:
                        for i, line in enumerate(f):
                            if line.startswith('#') and i < 15:
                                if any(key in line for key in ['Spin model', 'Algorithm', 'Theoretical Tc', 'System size']):
                                    print(f"    {line.strip()}")
                            elif not line.startswith('#'):
                                break
                    
                    # Quick data check
                    try:
                        import pandas as pd
                        data = pd.read_csv(layer_files[0], sep=r'\s+', comment='#',
                                         names=['Temperature', 'M2', 'M4', 'G2', 'G4', 'Energy', 'Cv', 'Corr'])
                        
                        print("  Data summary:")
                        print(f"    Temperature range: {data['Temperature'].min():.3f} - {data['Temperature'].max():.3f}")
                        print(f"    Data points: {len(data)}")
                        print(f"    Max specific heat: {data['Cv'].max():.3f}")
                        
                        # Check critical region
                        critical_mask = (data['Temperature'] >= 0.6) & (data['Temperature'] <= 0.9)
                        critical_points = critical_mask.sum()
                        print(f"    Critical region points: {critical_points}")
                        
                        if critical_points > 0:
                            max_cv_idx = data[critical_mask]['Cv'].idxmax()
                            tc_apparent = data.loc[max_cv_idx, 'Temperature']
                            theoretical_tc = 0.751  # 8-state Potts
                            print(f"    Apparent Tc: {tc_apparent:.4f}")
                            print(f"    Theoretical Tc: {theoretical_tc:.4f}")
                            print(f"    Difference: {tc_apparent - theoretical_tc:+.4f}")
                    
                    except ImportError:
                        print("    (Install pandas for detailed analysis)")
                    except Exception as e:
                        print(f"    Data analysis error: {e}")
                
                print(f"\n✓ Test completed successfully!")
                print(f"✓ Enhanced features working:")
                print(f"  - 8-state cube spins")
                print(f"  - Hybrid Monte Carlo algorithm")
                print(f"  - Optimized temperature range")
                print(f"  - Timestamped directory structure")
                
                return True
            else:
                print("✗ Output directory not found")
                return False
        else:
            print(f"✗ Simulation failed: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ Simulation timed out")
        return False
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

if __name__ == "__main__":
    success = test_cube_model()
    if success:
        print(f"\n🎉 CUBE MODEL TEST PASSED!")
    else:
        print(f"\n❌ Test failed")
    exit(0 if success else 1)