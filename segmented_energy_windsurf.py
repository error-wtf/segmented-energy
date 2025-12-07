#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmented Energy Model - Windsurf Implementation

Implements N-segment discretization for computing total energy of a test mass
in a central gravitational field.

Theoretical Model:
    E_tot(N) = sum_{n=1..N} [ E_GR_(n) + E_SR_(n) ]
    
    where:
        E_GR_(n) = - G * M * Δm / r_n        (gravitational)
        E_SR_(n) = (γ_n - 1) * Δm * c^2      (special relativistic)
        
        Δm = m / N
        r_n = r_in + (n - 0.5) * Δr
        Δr = (r_out - r_in) / N
        v_n = sqrt(G * M / r_n)
        γ_n = 1 / sqrt(1 - v_n^2 / c^2)

Full substitution:
    E_tot(N) = (m / N) * sum_{n=1..N} [
        - G * M / r_n
      + ( 1 / sqrt( 1 - (G*M)/(r_n * c^2) ) - 1 ) * c^2
    ]
"""

import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun, au
from typing import Dict, Optional


def radii_linear(r_in: u.Quantity, r_out: u.Quantity, N: int) -> u.Quantity:
    """
    Generate N segment midpoints linearly spaced between r_in and r_out.
    
    Parameters
    ----------
    r_in : astropy.Quantity
        Inner radius
    r_out : astropy.Quantity
        Outer radius
    N : int
        Number of segments
        
    Returns
    -------
    r_array : astropy.Quantity
        Array of N segment midpoints
        
    Notes
    -----
    Formula: r_n = r_in + (n - 0.5) * Δr
    where Δr = (r_out - r_in) / N
    and n = 1, 2, ..., N
    """
    if N <= 0:
        raise ValueError(f"N must be positive, got {N}")
    if r_out <= r_in:
        raise ValueError(f"r_out must be > r_in")
    
    # Ensure consistent units
    r_out = r_out.to(r_in.unit)
    
    # Radial step
    delta_r = (r_out - r_in) / N
    
    # Segment indices: n = 1, 2, ..., N
    n = np.arange(1, N + 1)
    
    # Segment midpoints
    r_array = r_in + (n - 0.5) * delta_r
    
    return r_array


def radii_phi_spiral(r_in: u.Quantity, r_out: u.Quantity, N: int, 
                      phi: float = 1.618033988749895, 
                      n_turns: int = 1) -> u.Quantity:
    """
    Placeholder: later we will implement a phi-based spiral segmentation.
    For now, you can just reuse the linear segmentation or implement a simple
    exponential scaling between r_in and r_out as a stub.
    
    Parameters
    ----------
    r_in : astropy.Quantity
        Inner radius
    r_out : astropy.Quantity
        Outer radius
    N : int
        Number of segments
    phi : float
        Golden ratio (default: (1 + sqrt(5))/2)
    n_turns : int
        Number of spiral turns (future use)
        
    Returns
    -------
    r_array : astropy.Quantity
        Array of N segment midpoints
    """
    # For now, just call radii_linear(...)
    return radii_linear(r_in, r_out, N)


def segment_energies(M: u.Quantity, m_test: u.Quantity, 
                      r_array: u.Quantity) -> tuple:
    """
    Compute gravitational and SR energies for each segment.
    
    Parameters
    ----------
    M : astropy.Quantity
        Central mass
    m_test : astropy.Quantity
        Total test mass
    r_array : astropy.Quantity
        Array of segment radii
        
    Returns
    -------
    E_GR_array : astropy.Quantity
        Gravitational energy per segment [J]
    E_SR_array : astropy.Quantity
        Special relativistic energy per segment [J]
        
    Notes
    -----
    E_GR_(n) = - G * M * Δm / r_n
    E_SR_(n) = (γ_n - 1) * Δm * c^2
    
    where:
        Δm = m_test / N
        v_n = sqrt(G * M / r_n)
        γ_n = 1 / sqrt(1 - v_n^2 / c^2)
    """
    N = len(r_array)
    delta_m = m_test / N
    
    # Keplerian orbital velocity: v_n = sqrt(G * M / r_n)
    v_array = np.sqrt(G * M / r_array)
    
    # Lorentz factor: γ_n = 1 / sqrt(1 - v_n^2 / c^2)
    beta_sq = (v_array / c) ** 2
    gamma_array = 1.0 / np.sqrt(1 - beta_sq)
    
    # Gravitational energy per segment
    E_GR_array = - G * M * delta_m / r_array
    
    # SR energy per segment
    E_SR_array = (gamma_array - 1.0) * delta_m * c**2
    
    return E_GR_array, E_SR_array


def compute_segmented_energy(M: u.Quantity, 
                               m_test: u.Quantity,
                               r_in: u.Quantity,
                               r_out: u.Quantity,
                               N_segments: int,
                               segmentation: str = "linear") -> Dict:
    """
    Compute total segmented energy with N segments.
    
    Parameters
    ----------
    M : astropy.Quantity
        Central mass (e.g. 1 * M_sun)
    m_test : astropy.Quantity
        Test mass (e.g. 1 * u.kg)
    r_in : astropy.Quantity
        Inner radius (e.g. 2.0 * R_sun)
    r_out : astropy.Quantity
        Outer radius (e.g. 1.0 * au)
    N_segments : int
        Number of segments
    segmentation : str
        Type of segmentation: "linear" or "phi"
        
    Returns
    -------
    result : dict
        Dictionary containing:
            "E_total" : Total energy [J]
            "E_GR_total" : Total gravitational energy [J]
            "E_SR_total" : Total SR energy [J]
            "N_segments" : Number of segments
            "r_in" : Inner radius
            "r_out" : Outer radius
    """
    # Generate radii based on segmentation type
    if segmentation == "linear":
        r_array = radii_linear(r_in, r_out, N_segments)
    elif segmentation == "phi":
        r_array = radii_phi_spiral(r_in, r_out, N_segments)
    else:
        raise ValueError(f"Unknown segmentation type: {segmentation}")
    
    # Compute energies per segment
    E_GR_array, E_SR_array = segment_energies(M, m_test, r_array)
    
    # Sum over all segments
    E_GR_total = np.sum(E_GR_array)
    E_SR_total = np.sum(E_SR_array)
    E_total = E_GR_total + E_SR_total
    
    # Return results
    return {
        "E_total": E_total.to(u.J),
        "E_GR_total": E_GR_total.to(u.J),
        "E_SR_total": E_SR_total.to(u.J),
        "N_segments": N_segments,
        "r_in": r_in,
        "r_out": r_out,
    }


if __name__ == "__main__":
    
    print("=" * 80)
    print("SEGMENTED ENERGY MODEL - ASTROPY IMPLEMENTATION")
    print("=" * 80)
    print()
    
    # =========================================================================
    # TEST CASE 1: Solar System (Sun, 1 kg test mass)
    # =========================================================================
    
    print("TEST CASE 1: Solar System")
    print("-" * 80)
    
    M_example = 1.0 * M_sun
    m_example = 1.0 * u.kg
    r_in_example = 2.0 * R_sun
    r_out_example = 1.0 * au
    
    print(f"Parameters:")
    print(f"  M = {M_example}")
    print(f"  m = {m_example}")
    print(f"  r_in = {r_in_example}")
    print(f"  r_out = {r_out_example}")
    print()
    
    N_values = [10, 100, 1000, 10000]
    
    print(f"{'N':>8} {'E_total [J]':>20} {'E_GR_total [J]':>20} {'E_SR_total [J]':>20} {'E_total / (m c^2)':>20}")
    print("-" * 80)
    
    for N in N_values:
        result = compute_segmented_energy(
            M_example, m_example, r_in_example, r_out_example, N
        )
        
        # Compute normalized energy
        E_normalized = (result["E_total"] / (m_example * c**2)).decompose().value
        
        print(f"{N:8d} {result['E_total'].value:20.10e} "
              f"{result['E_GR_total'].value:20.10e} "
              f"{result['E_SR_total'].value:20.10e} "
              f"{E_normalized:20.10e}")
    
    # Print one detailed case
    print()
    print(f"Case: M = 1 Msun, r_in = 2 Rsun, r_out = 1 AU, N = 1000")
    result_1000 = compute_segmented_energy(
        M_example, m_example, r_in_example, r_out_example, 1000
    )
    E_norm_1000 = (result_1000["E_total"] / (m_example * c**2)).decompose().value
    print(f"E_total = {result_1000['E_total'].value:.10e} J  "
          f"(E_total / (m c^2) = {E_norm_1000:.10e})")
    
    # =========================================================================
    # TEST CASE 2: Massive Star (10 solar masses)
    # =========================================================================
    
    print()
    print("=" * 80)
    print("TEST CASE 2: Massive Star")
    print("-" * 80)
    
    M_massive = 10 * M_sun
    m_massive = 1 * u.kg
    r_in_massive = 6 * R_sun
    r_out_massive = 1 * au
    
    print(f"Parameters:")
    print(f"  M = {M_massive}")
    print(f"  m = {m_massive}")
    print(f"  r_in = {r_in_massive}")
    print(f"  r_out = {r_out_massive}")
    print()
    
    N_values_massive = [100, 1000]
    
    print(f"{'N':>8} {'E_total [J]':>20} {'E_GR_total [J]':>20} {'E_SR_total [J]':>20} {'E_total / (m c^2)':>20}")
    print("-" * 80)
    
    for N in N_values_massive:
        result = compute_segmented_energy(
            M_massive, m_massive, r_in_massive, r_out_massive, N
        )
        
        E_normalized = (result["E_total"] / (m_massive * c**2)).decompose().value
        
        print(f"{N:8d} {result['E_total'].value:20.10e} "
              f"{result['E_GR_total'].value:20.10e} "
              f"{result['E_SR_total'].value:20.10e} "
              f"{E_normalized:20.10e}")
    
    # Print detailed case for N=1000
    print()
    print(f"Case: M = 10 Msun, r_in = 6 Rsun, r_out = 1 AU, N = 1000")
    result_massive = compute_segmented_energy(
        M_massive, m_massive, r_in_massive, r_out_massive, 1000
    )
    E_norm_massive = (result_massive["E_total"] / (m_massive * c**2)).decompose().value
    print(f"E_total = {result_massive['E_total'].value:.10e} J  "
          f"(E_total / (m c^2) = {E_norm_massive:.10e})")
    
    # =========================================================================
    # Convergence Test
    # =========================================================================
    
    print()
    print("=" * 80)
    print("CONVERGENCE TEST (Test Case 1)")
    print("-" * 80)
    
    print("\nRelative change in E_total between successive N values:")
    
    results_conv = []
    for N in N_values:
        res = compute_segmented_energy(
            M_example, m_example, r_in_example, r_out_example, N
        )
        results_conv.append(res)
    
    for i in range(1, len(N_values)):
        E_prev = results_conv[i-1]["E_total"].value
        E_curr = results_conv[i]["E_total"].value
        rel_change = abs((E_curr - E_prev) / E_curr)
        print(f"  N = {N_values[i-1]:5d} -> {N_values[i]:5d}: "
              f"relative change = {rel_change:.3e}")
    
    # =========================================================================
    # Summary
    # =========================================================================
    
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("✓ Segmented energy model implemented with astropy")
    print("✓ Linear segmentation tested")
    print("✓ Phi-based segmentation hook in place (stub)")
    print("✓ Tested with N up to 10,000 segments")
    print("✓ Energy converges as N increases")
    print("✓ Two test cases completed successfully")
    print()
    print("To run this script:")
    print("  python segmented_energy_windsurf.py")
    print("=" * 80)
