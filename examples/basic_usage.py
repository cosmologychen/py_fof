"""
CUBE FoF - Basic Usage Example

This example demonstrates the basic usage of the py_fof package
for finding halos in cosmological simulations.
"""

import numpy as np
from py_fof import fof


def example_single_halo():
    """
    Example 1: Single halo detection
    """
    print("="*60)
    print("Example 1: Single Halo Detection")
    print("="*60)
    
    # Generate test particles (clustered distribution)
    np.random.seed(42)
    n_particles = 1000
    box_size = 100.0
    
    # Create a clustered distribution
    xp = np.random.randn(3, n_particles).astype(np.float64) * 10 + 50
    
    print(f"Generated {n_particles} particles")
    print(f"Box size: {box_size}")
    
    # Run FOF
    in_halo = fof(
        xp,
        L_Box=512.0,
        linking_length=2.0,
        np_halo_min=10,
        verbose=1
    )
    
    n_in_halo = np.sum(in_halo)
    print(f"\nParticles in halo: {n_in_halo} / {n_particles}")
    print(f"Fraction: {n_in_halo / n_particles * 100:.1f}%")
    
    return xp, in_halo


def example_multiple_halos():
    """
    Example 2: Multiple halos detection
    """
    print("\n" + "="*60)
    print("Example 2: Multiple Halos Detection")
    print("="*60)
    
    np.random.seed(123)
    
    # Generate three separate halos
    halo1 = np.random.randn(3, 300).astype(np.float64) * 5 + np.array([[50], [50], [50]])
    halo2 = np.random.randn(3, 300).astype(np.float64) * 5 + np.array([[150], [150], [150]])
    halo3 = np.random.randn(3, 300).astype(np.float64) * 5 + np.array([[250], [250], [250]])
    
    xp_list = [halo1, halo2, halo3]
    
    print(f"Generated 3 halos with 300 particles each")
    
    # Run FOF
    in_halo = fof(
        xp_list,
        L_Box=512.0,
        linking_length=2.0,
        np_halo_min=10,
        verbose=1
    )
    
    total_particles = sum(xp.shape[1] for xp in xp_list)
    n_in_halo = np.sum(in_halo)
    
    print(f"\nTotal particles: {total_particles}")
    print(f"Particles in halos: {n_in_halo}")
    print(f"Fraction: {n_in_halo / total_particles * 100:.1f}%")
    
    return xp_list, in_halo


def example_with_binding_check():
    """
    Example 3: Halo detection with velocity-based binding check
    """
    print("\n" + "="*60)
    print("Example 3: With Velocity Binding Check")
    print("="*60)
    
    np.random.seed(456)
    
    # Generate halo particles
    n_particles = 500
    xp = np.random.randn(3, n_particles).astype(np.float64) * 8 + 100
    
    # Generate velocities (some high velocity particles)
    vp = np.random.randn(3, n_particles).astype(np.float64) * 0.5
    
    # Add some high velocity particles (will be excluded by binding check)
    n_high_v = 50
    high_v_indices = np.random.choice(n_particles, n_high_v, replace=False)
    vp[:, high_v_indices] *= 50  # High velocity
    
    print(f"Generated {n_particles} particles")
    print(f"High velocity particles: {n_high_v}")
    
    # Run FOF without binding check
    in_halo_no_check = fof(
        xp,
        L_Box=512.0,
        linking_length=2.0,
        np_halo_min=10,
        bound_check=False
    )
    
    # Run FOF with binding check
    in_halo_with_check = fof(
        xp,
        vp_list=vp,
        L_Box=512.0,
        linking_length=2.0,
        np_halo_min=10,
        bound_check=True,
        verbose=1
    )
    
    print(f"\nWithout binding check: {np.sum(in_halo_no_check)} particles")
    print(f"With binding check: {np.sum(in_halo_with_check)} particles")
    print(f"Excluded by binding check: {np.sum(in_halo_no_check) - np.sum(in_halo_with_check)}")
    
    return xp, vp, in_halo_with_check


def main():
    """
    Run all examples
    """
    print("\n" + "#"*60)
    print("# CUBE FoF - Basic Usage Examples")
    print("#"*60)
    
    # Example 1: Single halo
    xp1, in_halo1 = example_single_halo()
    
    # Example 2: Multiple halos
    xp_list, in_halo2 = example_multiple_halos()
    
    # Example 3: With binding check
    xp3, vp3, in_halo3 = example_with_binding_check()
    
    print("\n" + "#"*60)
    print("# All examples completed successfully!")
    print("#"*60 + "\n")


if __name__ == "__main__":
    main()
