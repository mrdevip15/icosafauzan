#!/usr/bin/env python3
"""
3D Visualization of Quasi-3D Lattice with 8-State Potts Spins
Visualizes the 4-layer stacked lattice structure and spin interactions
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def get_spin_directions():
    """Define the 8 spin directions (vertices of a cube, normalized)"""
    invsqrt3 = 1.0/np.sqrt(3.0)
    
    # 8 vertices of a cube: (±1,±1,±1) normalized
    spins = np.array([
        [ invsqrt3,  invsqrt3,  invsqrt3],  # (1,1,1)
        [-invsqrt3,  invsqrt3,  invsqrt3],  # (-1,1,1)
        [-invsqrt3, -invsqrt3,  invsqrt3],  # (-1,-1,1)
        [ invsqrt3, -invsqrt3,  invsqrt3],  # (1,-1,1)
        [ invsqrt3,  invsqrt3, -invsqrt3],  # (1,1,-1)
        [-invsqrt3,  invsqrt3, -invsqrt3],  # (-1,1,-1)
        [-invsqrt3, -invsqrt3, -invsqrt3],  # (-1,-1,-1)
        [ invsqrt3, -invsqrt3, -invsqrt3],  # (1,-1,-1)
    ])
    
    return spins

def create_lattice_positions(nx, ny, nz):
    """Create 3D positions for lattice sites"""
    positions = []
    site_indices = {}
    
    layer_size = nx * ny
    
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                # Linear index
                la = ix + iy * nx + iz * layer_size
                
                # 3D position (with some spacing between layers)
                x = ix
                y = iy  
                z = iz * 2.0  # Space layers apart for better visualization
                
                positions.append([x, y, z])
                site_indices[(ix, iy, iz)] = la
    
    return np.array(positions), site_indices

def get_neighbors(nx, ny, nz):
    """Calculate neighbor connections for the 4-layer lattice"""
    layer_size = nx * ny
    neighbors = {}
    
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                la = ix + iy * nx + iz * layer_size
                
                # In-plane neighbors (4 neighbors within the same layer)
                neighbors[la] = []
                
                # Right neighbor
                ix_right = (ix + 1) % nx
                la_right = ix_right + iy * nx + iz * layer_size
                neighbors[la].append(la_right)
                
                # Left neighbor  
                ix_left = (ix - 1 + nx) % nx
                la_left = ix_left + iy * nx + iz * layer_size
                neighbors[la].append(la_left)
                
                # Up neighbor
                iy_up = (iy + 1) % ny
                la_up = ix + iy_up * nx + iz * layer_size
                neighbors[la].append(la_up)
                
                # Down neighbor
                iy_down = (iy - 1 + ny) % ny
                la_down = ix + iy_down * nx + iz * layer_size
                neighbors[la].append(la_down)
                
                # Inter-layer neighbors (2 neighbors in adjacent layers)
                # Layer above
                iz_above = (iz + 1) % nz
                la_above = ix + iy * nx + iz_above * layer_size
                neighbors[la].append(la_above)
                
                # Layer below
                iz_below = (iz - 1 + nz) % nz
                la_below = ix + iy * nx + iz_below * layer_size
                neighbors[la].append(la_below)
    
    return neighbors

def visualize_lattice_structure(nx=4, ny=4, nz=4):
    """Visualize the 3D lattice structure"""
    
    # Create lattice positions
    positions, site_indices = create_lattice_positions(nx, ny, nz)
    neighbors = get_neighbors(nx, ny, nz)
    
    # Create figure
    fig = plt.figure(figsize=(15, 12))
    
    # Plot 1: Lattice structure with connections
    ax1 = fig.add_subplot(221, projection='3d')
    
    # Plot lattice sites
    colors = plt.cm.tab10(np.linspace(0, 1, nz))
    
    for iz in range(nz):
        layer_positions = positions[iz*nx*ny:(iz+1)*nx*ny]
        ax1.scatter(layer_positions[:, 0], layer_positions[:, 1], layer_positions[:, 2], 
                   c=[colors[iz]], s=100, alpha=0.8, label=f'Layer {iz}')
    
    # Draw connections (sample some to avoid clutter)
    sample_sites = list(range(0, len(positions), 4))  # Sample every 4th site
    
    for la in sample_sites:
        pos1 = positions[la]
        for neighbor_la in neighbors[la]:
            pos2 = positions[neighbor_la]
            
            # Different line styles for in-plane vs inter-layer connections
            if abs(pos1[2] - pos2[2]) < 0.1:  # In-plane connection
                ax1.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                        'b-', alpha=0.3, linewidth=0.5)
            else:  # Inter-layer connection
                ax1.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                        'r-', alpha=0.6, linewidth=1.0)
    
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z (Layer)')
    ax1.set_title(f'Quasi-3D Lattice Structure\n{nx}×{ny}×{nz} with Periodic BC')
    ax1.legend()
    
    # Plot 2: Spin directions (8-state Potts model)
    ax2 = fig.add_subplot(222, projection='3d')
    
    spin_directions = get_spin_directions()
    spin_colors = plt.cm.Set3(np.linspace(0, 1, 8))
    
    # Plot spin directions as arrows from origin
    for i, (spin, color) in enumerate(zip(spin_directions, spin_colors)):
        ax2.quiver(0, 0, 0, spin[0], spin[1], spin[2], 
                  color=color, arrow_length_ratio=0.1, linewidth=3, 
                  label=f'Spin {i}')
    
    ax2.set_xlabel('X Component')
    ax2.set_ylabel('Y Component') 
    ax2.set_zlabel('Z Component')
    ax2.set_title('8-State Potts Spin Directions\n(Cube Vertices)')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot 3: Single layer with neighbor connections
    ax3 = fig.add_subplot(223, projection='3d')
    
    # Show just one layer in detail
    layer = 0
    layer_positions = positions[layer*nx*ny:(layer+1)*nx*ny]
    
    # Plot sites
    ax3.scatter(layer_positions[:, 0], layer_positions[:, 1], layer_positions[:, 2], 
               c='blue', s=150, alpha=0.8)
    
    # Add site labels
    for i, pos in enumerate(layer_positions):
        ax3.text(pos[0], pos[1], pos[2]+0.1, str(i), fontsize=8)
    
    # Draw all in-plane connections for this layer
    for ix in range(nx):
        for iy in range(ny):
            la = ix + iy * nx + layer * nx * ny
            pos1 = positions[la]
            
            # Right connection
            ix_right = (ix + 1) % nx
            la_right = ix_right + iy * nx + layer * nx * ny
            pos2 = positions[la_right]
            ax3.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                    'b-', linewidth=1.5)
            
            # Up connection
            iy_up = (iy + 1) % ny
            la_up = ix + iy_up * nx + layer * nx * ny
            pos2 = positions[la_up]
            ax3.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                    'g-', linewidth=1.5)
    
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    ax3.set_title(f'Single Layer Detail\nLayer {layer} with In-Plane Connections')
    
    # Plot 4: Inter-layer connections
    ax4 = fig.add_subplot(224, projection='3d')
    
    # Show a few sites and their inter-layer connections
    sample_positions = positions[::4]  # Every 4th site
    sample_indices = list(range(0, len(positions), 4))
    
    ax4.scatter(sample_positions[:, 0], sample_positions[:, 1], sample_positions[:, 2], 
               c=sample_positions[:, 2], s=100, alpha=0.8, cmap='viridis')
    
    # Draw inter-layer connections
    for la in sample_indices:
        pos1 = positions[la]
        for neighbor_la in neighbors[la]:
            pos2 = positions[neighbor_la]
            
            # Only draw inter-layer connections
            if abs(pos1[2] - pos2[2]) > 0.1:
                ax4.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], [pos1[2], pos2[2]], 
                        'r-', linewidth=2.0, alpha=0.7)
    
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z (Layer)')
    ax4.set_title('Inter-Layer Connections\n(Red lines between layers)')
    
    plt.tight_layout()
    plt.savefig('quasi_3d_lattice_visualization.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig

def visualize_spin_configuration(nx=4, ny=4, nz=4, temperature=1.0):
    """Visualize a sample spin configuration"""
    
    # Create lattice positions
    positions, _ = create_lattice_positions(nx, ny, nz)
    spin_directions = get_spin_directions()
    
    # Generate a sample spin configuration (random for demonstration)
    np.random.seed(42)  # For reproducible results
    
    # At low temperature, prefer ordered state; at high temperature, more random
    if temperature < 0.5:
        # Mostly aligned spins (ferromagnetic state)
        spin_states = np.random.choice([0, 1], size=len(positions), p=[0.8, 0.2])
    elif temperature < 1.0:
        # Mixed state
        spin_states = np.random.choice(8, size=len(positions))
    else:
        # High temperature - random
        spin_states = np.random.choice(8, size=len(positions))
    
    # Create figure
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Color map for different spin states
    colors = plt.cm.Set3(np.linspace(0, 1, 8))
    
    # Plot each site with its spin direction
    for i, (pos, spin_state) in enumerate(zip(positions, spin_states)):
        # Plot site
        ax.scatter(pos[0], pos[1], pos[2], c=[colors[spin_state]], s=100, alpha=0.8)
        
        # Plot spin direction as arrow
        spin_dir = spin_directions[spin_state] * 0.3  # Scale arrow size
        ax.quiver(pos[0], pos[1], pos[2], 
                 spin_dir[0], spin_dir[1], spin_dir[2], 
                 color=colors[spin_state], alpha=0.7, arrow_length_ratio=0.3)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z (Layer)')
    ax.set_title(f'Sample Spin Configuration\nT = {temperature:.1f}, 8-State Potts Model')
    
    # Create legend for spin states
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=colors[i], markersize=10, 
                                 label=f'Spin {i}') for i in range(8)]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f'spin_configuration_T{temperature:.1f}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig

if __name__ == "__main__":
    print("Creating 3D lattice visualization...")
    
    # Visualize lattice structure
    fig1 = visualize_lattice_structure(nx=4, ny=4, nz=4)
    
    print("Creating spin configuration visualizations...")
    
    # Visualize spin configurations at different temperatures
    fig2 = visualize_spin_configuration(nx=4, ny=4, nz=4, temperature=0.3)  # Low T
    fig3 = visualize_spin_configuration(nx=4, ny=4, nz=4, temperature=1.5)  # High T
    
    print("Visualizations saved as:")
    print("- quasi_3d_lattice_visualization.png")
    print("- spin_configuration_T0.3.png") 
    print("- spin_configuration_T1.5.png")