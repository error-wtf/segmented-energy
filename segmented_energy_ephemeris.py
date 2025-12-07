#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmented Energy Model with Real Ephemeris Data

Combines Special Relativity (SR) and General Relativity (GR) energy contributions
using segmented spacetime approach with real astronomical data.

Theory:
    E_tot = E_rest + E_GR + E_SR
    
    where:
        E_rest = m * c²
        E_GR = Σ_{n=1..N} (γ_GR_n - 1) * m * c² * w_GR_n
        E_SR = Σ_{n=1..N} (γ_SR_n - 1) * m * c² * w_SR_n
        
        γ_SR_n = 1 / sqrt(1 - v_n²/c²)
        γ_GR_n = 1 / sqrt(1 - 2GM/(r_n*c²))

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import os
import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
import matplotlib.pyplot as plt
from typing import Tuple, Dict

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'


# =============================================================================
# Lorentz Factors
# =============================================================================

def gamma_sr(v: u.Quantity) -> np.ndarray:
    """
    Special Relativistic Lorentz factor.
    
    γ_SR = 1 / sqrt(1 - v²/c²)
    
    Parameters
    ----------
    v : astropy.Quantity
        Velocity (with units, e.g. m/s)
        
    Returns
    -------
    gamma_SR : numpy.ndarray
        Lorentz factor (dimensionless)
        
    Notes
    -----
    For orbital velocities v << c, γ_SR ≈ 1 + (1/2)v²/c²
    """
    beta_sq = (v / c)**2
    gamma = 1.0 / np.sqrt(1 - beta_sq)
    return gamma.value


def gamma_gr(M: u.Quantity, r: u.Quantity) -> np.ndarray:
    """
    General Relativistic gamma factor (Schwarzschild time dilation).
    
    γ_GR = 1 / sqrt(1 - 2GM/(rc²))
    
    Parameters
    ----------
    M : astropy.Quantity
        Central mass (with units, e.g. kg)
    r : astropy.Quantity
        Radial distance (with units, e.g. m)
        
    Returns
    -------
    gamma_GR : numpy.ndarray
        GR gamma factor (dimensionless)
        
    Notes
    -----
    This is the time dilation factor for a static observer
    in the Schwarzschild metric.
    """
    factor = (2 * G * M / (r * c**2)).decompose()
    gamma = 1.0 / np.sqrt(1 - factor)
    return gamma.value


# =============================================================================
# Segmentation
# =============================================================================

def compute_segment_arrays(M: u.Quantity, 
                            m: u.Quantity,
                            r_min: u.Quantity, 
                            r_max: u.Quantity,
                            N_SEG: int,
                            v_current: u.Quantity,
                            segmentation: str = "linear") -> Dict:
    """
    Compute arrays for segmented energy calculation.
    
    Parameters
    ----------
    M : astropy.Quantity
        Central mass
    m : astropy.Quantity
        Test mass
    r_min : astropy.Quantity
        Minimum radius
    r_max : astropy.Quantity
        Maximum radius (current position)
    N_SEG : int
        Number of segments
    v_current : astropy.Quantity
        Current orbital velocity
    segmentation : str
        "linear" or "logarithmic"
        
    Returns
    -------
    result : dict
        Dictionary containing:
        - r_n: segment radii
        - v_n: segment velocities
        - gamma_SR_n: SR Lorentz factors
        - gamma_GR_n: GR gamma factors
        - w_SR_n: SR weights
        - w_GR_n: GR weights
    """
    # Generate segment radii
    if segmentation == "linear":
        # Midpoint rule: r_n = r_min + (n - 0.5) * Δr
        n = np.arange(1, N_SEG + 1)
        delta_r = (r_max - r_min) / N_SEG
        r_n = r_min + (n - 0.5) * delta_r
    
    elif segmentation == "logarithmic":
        # Logarithmic spacing (better for large ranges)
        r_n = np.geomspace(r_min.value, r_max.value, N_SEG) * r_min.unit
    
    else:
        raise ValueError(f"Unknown segmentation: {segmentation}")
    
    # Velocity: use constant v_current for all segments (baseline)
    # Future: could implement v_n = sqrt(G*M/r_n) for Keplerian orbits
    v_n = np.full(N_SEG, v_current.value) * v_current.unit
    
    # Compute gamma factors
    gamma_SR_n = gamma_sr(v_n)
    gamma_GR_n = gamma_gr(M, r_n)
    
    # Uniform weights (can be modified later for phi-segmentation)
    w_SR_n = np.ones(N_SEG) / N_SEG
    w_GR_n = np.ones(N_SEG) / N_SEG
    
    return {
        'r_n': r_n,
        'v_n': v_n,
        'gamma_SR_n': gamma_SR_n,
        'gamma_GR_n': gamma_GR_n,
        'w_SR_n': w_SR_n,
        'w_GR_n': w_GR_n,
        'N_SEG': N_SEG,
        'segmentation': segmentation,
    }


# =============================================================================
# Energy Computation
# =============================================================================

def compute_segment_energies(m: u.Quantity,
                               gamma_SR_n: np.ndarray,
                               gamma_GR_n: np.ndarray,
                               w_SR_n: np.ndarray,
                               w_GR_n: np.ndarray) -> Tuple[u.Quantity, u.Quantity, 
                                                              u.Quantity, u.Quantity]:
    """
    Compute energy contributions from segments.
    
    E_SR_n = (γ_SR_n - 1) * m * c² * w_SR_n
    E_GR_n = (γ_GR_n - 1) * m * c² * w_GR_n
    
    Parameters
    ----------
    m : astropy.Quantity
        Test mass
    gamma_SR_n : numpy.ndarray
        SR Lorentz factors
    gamma_GR_n : numpy.ndarray
        GR gamma factors
    w_SR_n : numpy.ndarray
        SR weights
    w_GR_n : numpy.ndarray
        GR weights
        
    Returns
    -------
    E_SR_n : astropy.Quantity
        SR energy per segment
    E_GR_n : astropy.Quantity
        GR energy per segment
    E_SR : astropy.Quantity
        Total SR energy
    E_GR : astropy.Quantity
        Total GR energy
    """
    # Energy per segment
    E_SR_n = (gamma_SR_n - 1.0) * m * c**2 * w_SR_n
    E_GR_n = (gamma_GR_n - 1.0) * m * c**2 * w_GR_n
    
    # Total energies
    E_SR = np.sum(E_SR_n)
    E_GR = np.sum(E_GR_n)
    
    return E_SR_n, E_GR_n, E_SR, E_GR


def compute_total_energy(m: u.Quantity,
                          E_SR: u.Quantity,
                          E_GR: u.Quantity) -> Tuple[u.Quantity, u.Quantity]:
    """
    Compute total energy.
    
    E_tot = E_rest + E_SR + E_GR
    
    Parameters
    ----------
    m : astropy.Quantity
        Test mass
    E_SR : astropy.Quantity
        Total SR energy
    E_GR : astropy.Quantity
        Total GR energy
        
    Returns
    -------
    E_rest : astropy.Quantity
        Rest energy
    E_tot : astropy.Quantity
        Total energy
    """
    E_rest = m * c**2
    E_tot = E_rest + E_SR + E_GR
    
    return E_rest, E_tot


# =============================================================================
# Real Ephemeris Data
# =============================================================================

def get_earth_ephemeris(epoch: str = "2025-01-01T00:00:00") -> Dict:
    """
    Get Earth's position and velocity at a given epoch.
    
    Parameters
    ----------
    epoch : str
        ISO format date-time string
        
    Returns
    -------
    ephemeris : dict
        Dictionary containing:
        - time: astropy Time object
        - position: barycentric position
        - velocity: barycentric velocity
        - r: distance from barycenter (≈ Sun)
        - v: speed
    """
    # Create time object
    time = Time(epoch, scale='tdb')
    
    # Get Earth's barycentric position and velocity
    pos, vel = get_body_barycentric_posvel('earth', time)
    
    # Convert to Cartesian
    pos_cart = pos.xyz
    vel_cart = vel.xyz
    
    # Compute magnitudes
    r = np.sqrt(np.sum(pos_cart**2))
    v = np.sqrt(np.sum(vel_cart**2))
    
    return {
        'time': time,
        'position': pos_cart,
        'velocity': vel_cart,
        'r': r.to(u.m),
        'v': v.to(u.m/u.s),
    }


# =============================================================================
# Reporting
# =============================================================================

def print_summary(ephemeris: Dict, segments: Dict, 
                   E_rest: u.Quantity, E_SR: u.Quantity, 
                   E_GR: u.Quantity, E_tot: u.Quantity,
                   m: u.Quantity):
    """Print summary of computation."""
    
    print("="*80)
    print("SEGMENTED ENERGY MODEL WITH REAL EPHEMERIS DATA")
    print("="*80)
    print()
    
    print("SCENARIO:")
    print("-"*80)
    print(f"  Body: Earth")
    print(f"  Epoch: {ephemeris['time'].iso}")
    print(f"  Central Mass: {M_sun.to(u.kg):.6e}")
    print(f"  Test Mass: {m}")
    print()
    
    print("EPHEMERIS:")
    print("-"*80)
    print(f"  Distance (r): {ephemeris['r'].to(u.AU):.6f}")
    print(f"  Velocity (v): {ephemeris['v'].to(u.km/u.s):.6f}")
    print()
    
    print("SEGMENTATION:")
    print("-"*80)
    print(f"  Type: {segments['segmentation']}")
    print(f"  Number of segments (N): {segments['N_SEG']}")
    print(f"  r_min: {segments['r_n'][0].to(u.AU):.6f}")
    print(f"  r_max: {segments['r_n'][-1].to(u.AU):.6f}")
    print()
    
    print("ENERGIES:")
    print("-"*80)
    print(f"  E_rest = {E_rest.to(u.J):.12e} J")
    print(f"  E_SR   = {E_SR.to(u.J):.12e} J")
    print(f"  E_GR   = {E_GR.to(u.J):.12e} J")
    print(f"  E_tot  = {E_tot.to(u.J):.12e} J")
    print()
    print(f"  E_SR / E_rest = {(E_SR/E_rest).decompose().value:.12e}")
    print(f"  E_GR / E_rest = {(E_GR/E_rest).decompose().value:.12e}")
    print(f"  E_tot / E_rest = {(E_tot/E_rest).decompose().value:.12f}")
    print()


def print_segment_table(segments: Dict, E_SR_n: u.Quantity, E_GR_n: u.Quantity):
    """Print table of first and last segments."""
    
    print("SEGMENT TABLE (first 5 and last 5):")
    print("-"*80)
    print(f"{'n':>5} {'r_n [AU]':>12} {'v_n [km/s]':>12} "
          f"{'gamma_SR':>12} {'gamma_GR':>12} "
          f"{'E_SR_n [J]':>15} {'E_GR_n [J]':>15}")
    print("-"*80)
    
    N = segments['N_SEG']
    
    # First 5
    for i in range(min(5, N)):
        print(f"{i+1:5d} "
              f"{segments['r_n'][i].to(u.AU).value:12.6f} "
              f"{segments['v_n'][i].to(u.km/u.s).value:12.6f} "
              f"{segments['gamma_SR_n'][i]:12.9f} "
              f"{segments['gamma_GR_n'][i]:12.9f} "
              f"{E_SR_n[i].to(u.J).value:15.6e} "
              f"{E_GR_n[i].to(u.J).value:15.6e}")
    
    if N > 10:
        print("  ...")
        
        # Last 5
        for i in range(max(N-5, 5), N):
            print(f"{i+1:5d} "
                  f"{segments['r_n'][i].to(u.AU).value:12.6f} "
                  f"{segments['v_n'][i].to(u.km/u.s).value:12.6f} "
                  f"{segments['gamma_SR_n'][i]:12.9f} "
                  f"{segments['gamma_GR_n'][i]:12.9f} "
                  f"{E_SR_n[i].to(u.J).value:15.6e} "
                  f"{E_GR_n[i].to(u.J).value:15.6e}")
    
    print("="*80)


# =============================================================================
# Plotting
# =============================================================================

def plot_results(segments: Dict, E_SR_n: u.Quantity, E_GR_n: u.Quantity):
    """Create plots of segment properties."""
    
    r_AU = segments['r_n'].to(u.AU).value
    
    # Figure 1: Gamma factors
    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    ax1.plot(r_AU, segments['gamma_SR_n'], 'b-', linewidth=2, label='gamma_SR')
    ax1.set_xlabel('r [AU]')
    ax1.set_ylabel('gamma_SR')
    ax1.set_title('Special Relativistic Lorentz Factor')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.plot(r_AU, segments['gamma_GR_n'], 'r-', linewidth=2, label='gamma_GR')
    ax2.set_xlabel('r [AU]')
    ax2.set_ylabel('gamma_GR')
    ax2.set_title('General Relativistic Gamma Factor')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('gammas_vs_r.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: gammas_vs_r.png")
    
    # Figure 2: Energies per segment
    fig2, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(r_AU, E_SR_n.to(u.J).value, 'b-', linewidth=2, label='E_SR_n')
    ax.plot(r_AU, E_GR_n.to(u.J).value, 'r-', linewidth=2, label='E_GR_n')
    ax.set_xlabel('r [AU]')
    ax.set_ylabel('Energy per segment [J]')
    ax.set_title('Energy Contributions per Segment')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('energies_vs_r.png', dpi=150, bbox_inches='tight')
    print("Plot saved: energies_vs_r.png")
    
    plt.show()


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Main function to run the segmented energy computation."""
    
    # =========================================================================
    # Configuration
    # =========================================================================
    
    # Epoch for computation
    epoch = "2025-01-01T00:00:00"
    
    # Test mass
    m = 1.0 * u.kg
    
    # Central mass (Sun)
    M = M_sun
    
    # Number of segments
    N_SEG = 100
    
    # Segmentation type
    segmentation = "linear"  # or "logarithmic"
    
    # =========================================================================
    # Get real ephemeris data
    # =========================================================================
    
    print("Fetching Earth ephemeris data...")
    ephemeris = get_earth_ephemeris(epoch)
    
    r_current = ephemeris['r']
    v_current = ephemeris['v']
    
    # Define segmentation range
    r_min = 0.1 * r_current
    r_max = r_current
    
    # =========================================================================
    # Compute segments
    # =========================================================================
    
    print("Computing segment arrays...")
    segments = compute_segment_arrays(
        M=M,
        m=m,
        r_min=r_min,
        r_max=r_max,
        N_SEG=N_SEG,
        v_current=v_current,
        segmentation=segmentation
    )
    
    # =========================================================================
    # Compute energies
    # =========================================================================
    
    print("Computing energies...")
    E_SR_n, E_GR_n, E_SR, E_GR = compute_segment_energies(
        m=m,
        gamma_SR_n=segments['gamma_SR_n'],
        gamma_GR_n=segments['gamma_GR_n'],
        w_SR_n=segments['w_SR_n'],
        w_GR_n=segments['w_GR_n']
    )
    
    E_rest, E_tot = compute_total_energy(m, E_SR, E_GR)
    
    # =========================================================================
    # Output results
    # =========================================================================
    
    print("\n")
    print_summary(ephemeris, segments, E_rest, E_SR, E_GR, E_tot, m)
    print_segment_table(segments, E_SR_n, E_GR_n)
    
    # =========================================================================
    # Create plots
    # =========================================================================
    
    print("\nGenerating plots...")
    plot_results(segments, E_SR_n, E_GR_n)
    
    print("\n" + "="*80)
    print("COMPUTATION COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()
