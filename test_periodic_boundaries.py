#!/usr/bin/env python3
"""
Test and verify periodic boundary conditions for the 3-layer simulation
"""

import numpy as np
import matplotlib.pyplot as plt

def test_periodic_boundaries(nx=4, ny=4, nz=3):
    """Test the periodic boundary condition implementation"""
    
    layer_size = nx * ny
    total_sites = nx * ny * nz
    
    print(f"Testing periodic boundaries for {nx}×{ny}×{nz} lattice")
    print(f"Total sites: {total_sites}")
    print()
    
    # Test function to convert 3D coordinates to linear index
    def coord_to_index(ix, iy, iz):
        return ix + iy * nx + iz * layer_size
    
    # Test function to convert linear index to 3D coordinates
    def index_to_coord(la):
        iz = la // layer_size
        layer_pos = la % layer_size
        iy = layer_pos // nx
        ix = layer_pos % nx
        return ix, iy, iz
    
    # Test coordinate conversion
    print("Testing coordinate conversion:")
    for la in range(min(12, total_sites)):
        ix, iy, iz = index_to_coord(la)
        la_reconstructed = coord_to_index(ix, iy, iz)
        print(f"Site {la:2d}: ({ix},{iy},{iz}) -> {la_reconstructed}")
    print()
    
    # Calculate neighbors for each site
    neighbors = {}
    
    for la in range(total_sites):
        ix, iy, iz = index_to_coord(la)
        
        # In-plane neighbors (4 neighbors within the same layer)
        ix_right = (ix + 1) % nx
        ix_left = (ix - 1 + nx) % nx
        iy_up = (iy + 1) % ny
        iy_down = (iy - 1 + ny) % ny
        
        # Inter-layer neighbors (2 neighbors in adjacent layers)
        iz_above = (iz + 1) % nz
        iz_below = (iz - 1 + nz) % nz
        
        neighbors[la] = {
            'right': coord_to_index(ix_right, iy, iz),
            'left': coord_to_index(ix_left, iy, iz),
            'up': coord_to_index(ix, iy_up, iz),
            'down': coord_to_index(ix, iy_down, iz),
            'above': coord_to_index(ix, iy, iz_above),
            'below': coord_to_index(ix, iy, iz_below)
        }
    
    # Test boundary wrapping
    print("Testing boundary wrapping:")
    print("=" * 50)
    
    # Test a few edge cases
    test_sites = [
        (0, 0, 0),      # corner site
        (nx-1, 0, 0),   # right edge
        (0, ny-1, 0),   # top edge
        (nx-1, ny-1, 0), # top-right corner
        (0, 0, nz-1),   # corner on last layer
    ]
    
    for ix, iy, iz in test_sites:
        la = coord_to_index(ix, iy, iz)
        neigh = neighbors[la]
        
        print(f"Site ({ix},{iy},{iz}) [index {la}]:")
        for direction, neighbor_la in neigh.items():
            nix, niy, niz = index_to_coord(neighbor_la)
            print(f"  {direction:>5}: site {neighbor_la:2d} at ({nix},{niy},{niz})")
        print()
    
    # Verify reciprocal connections
    print("Verifying reciprocal connections:")
    print("=" * 40)
    
    errors = 0
    for la in range(total_sites):
        neigh = neighbors[la]
        
        # Check if right neighbor's left neighbor is this site
        right_neighbor = neigh['right']
        if neighbors[right_neighbor]['left'] != la:
            print(f"ERROR: Site {la} -> right -> left != {la}")
            errors += 1
        
        # Check if up neighbor's down neighbor is this site
        up_neighbor = neigh['up']
        if neighbors[up_neighbor]['down'] != la:
            print(f"ERROR: Site {la} -> up -> down != {la}")
            errors += 1
        
        # Check if above neighbor's below neighbor is this site
        above_neighbor = neigh['above']
        if neighbors[above_neighbor]['below'] != la:
            print(f"ERROR: Site {la} -> above -> below != {la}")
            errors += 1
    
    if errors == 0:
        print("✓ All reciprocal connections verified successfully!")
    else:
        print(f"✗ Found {errors} connection errors")
    
    print()
    
    # Test layer connectivity
    print("Testing inter-layer connectivity:")
    print("=" * 35)
    
    for iz in range(nz):
        above_layer = (iz + 1) % nz
        below_layer = (iz - 1 + nz) % nz
        print(f"Layer {iz}: connects to layer {above_layer} (above) and {below_layer} (below)")
    
    print()
    
    # Visualize a small section
    if nx <= 4 and ny <= 4:
        visualize_connections(nx, ny, nz, neighbors)
    
    return neighbors

def visualize_connections(nx, ny, nz, neighbors):
    """Visualize the neighbor connections for a small lattice"""
    
    print("Neighbor connections for each site:")
    print("=" * 40)
    
    layer_size = nx * ny
    
    for iz in range(nz):
        print(f"\nLayer {iz}:")
        print("  " + "  ".join([f"{i:2d}" for i in range(nx)]))
        
        for iy in range(ny):
            row = f"{iy:1d} "
            for ix in range(nx):
                la = ix + iy * nx + iz * layer_size
                row += f"{la:2d} "
            print(row)
        
        # Show some neighbor examples for this layer
        if iz == 0:  # Only show for first layer to avoid clutter
            print("\n  Sample neighbor connections:")
            for sample_site in [0, nx-1, layer_size-nx, layer_size-1]:
                if sample_site < layer_size:
                    neigh = neighbors[sample_site]
                    print(f"    Site {sample_site}: R={neigh['right']}, L={neigh['left']}, U={neigh['up']}, D={neigh['down']}, A={neigh['above']}, B={neigh['below']}")

def run_boundary_test():
    """Run comprehensive boundary condition tests"""
    
    print("Periodic Boundary Condition Test")
    print("=" * 50)
    print()
    
    # Test different lattice sizes
    test_cases = [
        (4, 4, 3),  # Standard test case
        (2, 2, 3),  # Minimal case
        (6, 4, 3),  # Rectangular
    ]
    
    for nx, ny, nz in test_cases:
        print(f"\nTesting {nx}×{ny}×{nz} lattice:")
        print("-" * 30)
        neighbors = test_periodic_boundaries(nx, ny, nz)
        print()

if __name__ == "__main__":
    run_boundary_test()