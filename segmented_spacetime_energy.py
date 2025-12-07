#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmented Spacetime Energy Model

Physikalisches Modell:
----------------------
Berechnet die Energie eines Testteilchens in einem segmentierten Raumzeitmodell
um einen Zentralkörper. Berücksichtigt sowohl spezielle Relativität (SR)
durch Bahngeschwindigkeit als auch allgemeine Relativität (GR) durch
Schwarzschild-Zeitdilatation.

Das Modell teilt den Radialbereich [r_min, r_max] in N Segmente und berechnet
für jedes Segment die SR- und GR-Lorentzfaktoren sowie die zugehörigen Energien.

Usage:
------
    python segmented_spacetime_energy.py

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import numpy as np
import astropy.units as u
from astropy.constants import G, c, M_sun, R_sun, M_earth, R_earth
import matplotlib.pyplot as plt
from typing import Tuple


# =============================================================================
# Hilfsfunktionen
# =============================================================================

def segmented_radii(N: int, r_min: u.Quantity, r_max: u.Quantity) -> u.Quantity:
    """
    Erzeugt N gleichmäßig verteilte Radien zwischen r_min und r_max.
    
    Parameters
    ----------
    N : int
        Anzahl der Segmente
    r_min : astropy.Quantity
        Minimaler Radius (mit Längeneinheit)
    r_max : astropy.Quantity
        Maximaler Radius (mit Längeneinheit)
        
    Returns
    -------
    r_array : astropy.Quantity
        Array von N Radien [r_0, r_1, ..., r_{N-1}]
        
    Notes
    -----
    Radialsegmentierung:
        Δr = (r_max - r_min) / (N - 1)
        r_n = r_min + n * Δr   für n = 0, 1, ..., N-1
    """
    if N < 2:
        raise ValueError("N muss mindestens 2 sein")
    if r_max <= r_min:
        raise ValueError("r_max muss größer als r_min sein")
    
    # Stelle sicher, dass beide Radien die gleiche Einheit haben
    r_max = r_max.to(r_min.unit)
    
    # Erzeuge linearen Abstand
    delta_r = (r_max - r_min) / (N - 1)
    
    # Segment-Indizes: n = 0, 1, ..., N-1
    n = np.arange(N)
    
    # Radien: r_n = r_min + n * Δr
    r_array = r_min + n * delta_r
    
    return r_array


def gamma_sr(M: u.Quantity, r: u.Quantity) -> np.ndarray:
    """
    Berechnet den spezial-relativistischen Lorentzfaktor aus Bahngeschwindigkeit.
    
    Parameters
    ----------
    M : astropy.Quantity
        Masse des Zentralkörpers
    r : astropy.Quantity
        Radius oder Array von Radien
        
    Returns
    -------
    gamma_sr : numpy.ndarray
        Spezial-relativistischer Lorentzfaktor (dimensionslos)
        
    Notes
    -----
    Formeln:
        v²(r) = G * M / r                    (Kreisbahn-Geschwindigkeit)
        γ_sr = 1 / sqrt(1 - v² / c²)
             = 1 / sqrt(1 - (G*M) / (r*c²))
    """
    # Berechne v²/c² = (G*M) / (r*c²)
    beta_sq = (G * M / (r * c**2)).decompose()
    
    # γ_sr = 1 / sqrt(1 - β²)
    gamma_sr_val = 1.0 / np.sqrt(1 - beta_sq)
    
    return gamma_sr_val.value  # Dimensionslos


def gamma_gr(M: u.Quantity, r: u.Quantity) -> np.ndarray:
    """
    Berechnet den allgemein-relativistischen Lorentzfaktor (Schwarzschild).
    
    Parameters
    ----------
    M : astropy.Quantity
        Masse des Zentralkörpers
    r : astropy.Quantity
        Radius oder Array von Radien
        
    Returns
    -------
    gamma_gr : numpy.ndarray
        Allgemein-relativistischer Lorentzfaktor (dimensionslos)
        
    Notes
    -----
    Formel:
        γ_gr = 1 / sqrt(1 - 2*G*M / (c²*r))
        
    Dies ist der Zeitdilatationsfaktor für einen stationären Beobachter
    im Schwarzschild-Feld.
    """
    # Berechne 2*G*M / (c²*r)
    factor = (2 * G * M / (c**2 * r)).decompose()
    
    # γ_gr = 1 / sqrt(1 - factor)
    gamma_gr_val = 1.0 / np.sqrt(1 - factor)
    
    return gamma_gr_val.value  # Dimensionslos


def energies_segment(M: u.Quantity, m: u.Quantity, 
                      r: u.Quantity) -> Tuple[u.Quantity, u.Quantity, 
                                               u.Quantity, u.Quantity,
                                               np.ndarray, np.ndarray]:
    """
    Berechnet die Energien und Lorentzfaktoren für Segment(e).
    
    Parameters
    ----------
    M : astropy.Quantity
        Masse des Zentralkörpers
    m : astropy.Quantity
        Masse des Testteilchens
    r : astropy.Quantity
        Radius oder Array von Radien
        
    Returns
    -------
    E_rest : astropy.Quantity
        Ruheenergie m*c²
    E_sr : astropy.Quantity
        SR-Energieanteil (γ_sr - 1) * m*c²
    E_gr : astropy.Quantity
        GR-Energieanteil (γ_gr - 1) * m*c²
    E_tot : astropy.Quantity
        Gesamtenergie E_rest + E_sr + E_gr
    gamma_sr_vals : numpy.ndarray
        SR-Lorentzfaktoren (dimensionslos)
    gamma_gr_vals : numpy.ndarray
        GR-Lorentzfaktoren (dimensionslos)
        
    Notes
    -----
    Formeln:
        E_rest = m * c²
        E_sr(n) = (γ_sr(n) - 1) * m * c²
        E_gr(n) = (γ_gr(n) - 1) * m * c²
        E_tot(n) = E_rest + E_sr(n) + E_gr(n)
                 = m*c² * [γ_sr(n) + γ_gr(n) - 1]
    """
    # Ruheenergie
    E_rest = m * c**2
    
    # Berechne Lorentzfaktoren
    gamma_sr_vals = gamma_sr(M, r)
    gamma_gr_vals = gamma_gr(M, r)
    
    # SR-Energieanteil: (γ_sr - 1) * m*c²
    E_sr = (gamma_sr_vals - 1.0) * E_rest
    
    # GR-Energieanteil: (γ_gr - 1) * m*c²
    E_gr = (gamma_gr_vals - 1.0) * E_rest
    
    # Gesamtenergie
    E_tot = E_rest + E_sr + E_gr
    
    return E_rest, E_sr, E_gr, E_tot, gamma_sr_vals, gamma_gr_vals


# =============================================================================
# Hauptprogramm
# =============================================================================

if __name__ == "__main__":
    
    print("=" * 80)
    print("SEGMENTED SPACETIME ENERGY MODEL")
    print("=" * 80)
    print()
    
    # =========================================================================
    # Parameter definieren
    # =========================================================================
    
    print("PARAMETER:")
    print("-" * 80)
    
    # Zentralkörper: Sonne
    M = M_sun
    print(f"Zentralkörper-Masse M = {M}")
    
    # Schwarzschild-Radius: r_s = 2*G*M/c²
    r_s = (2 * G * M / c**2).to(u.km)
    print(f"Schwarzschild-Radius r_s = {r_s:.3f}")
    
    # Minimaler Radius: 3 Gravitationsradien (1.5 r_s)
    r_min = 3 * G * M / c**2
    r_min = r_min.to(u.km)
    print(f"Minimaler Radius r_min = {r_min:.3f} = {r_min/r_s:.2f} r_s")
    
    # Maximaler Radius: 100 * r_min
    r_max = 100 * r_min
    print(f"Maximaler Radius r_max = {r_max:.3f} = {r_max/r_s:.2f} r_s")
    
    # Testteilchen-Masse
    m = 1.0 * u.kg
    print(f"Testteilchen-Masse m = {m}")
    
    # Anzahl Segmente
    N = 50
    print(f"Anzahl Segmente N = {N}")
    print()
    
    # =========================================================================
    # Radien berechnen
    # =========================================================================
    
    print("BERECHNUNG:")
    print("-" * 80)
    
    r_array = segmented_radii(N, r_min, r_max)
    print(f"✓ {N} Radien berechnet: {r_array[0]:.3f} ... {r_array[-1]:.3f}")
    
    # =========================================================================
    # Energien und Lorentzfaktoren berechnen
    # =========================================================================
    
    E_rest, E_sr, E_gr, E_tot, gamma_sr_vals, gamma_gr_vals = energies_segment(
        M, m, r_array
    )
    
    print(f"✓ Lorentzfaktoren und Energien berechnet")
    print()
    
    # =========================================================================
    # Plausibilitätstests
    # =========================================================================
    
    print("PLAUSIBILITÄTSTESTS:")
    print("-" * 80)
    
    # Test 1: γ >= 1
    assert np.all(gamma_sr_vals >= 1.0), "gamma_sr muss >= 1 sein"
    assert np.all(gamma_gr_vals >= 1.0), "gamma_gr muss >= 1 sein"
    print("✓ Alle γ_sr >= 1 und γ_gr >= 1")
    
    # Test 2: γ -> 1 für große Radien
    gamma_sr_max = gamma_sr_vals[-1]
    gamma_gr_max = gamma_gr_vals[-1]
    print(f"✓ Bei r_max: γ_sr = {gamma_sr_max:.6f}, γ_gr = {gamma_gr_max:.6f}")
    print(f"  (beide nähern sich 1 für r → ∞)")
    
    # Test 3: Energie-Verhältnisse werden kleiner
    E_sr_normalized = (E_sr / (m * c**2)).decompose().value
    E_gr_normalized = (E_gr / (m * c**2)).decompose().value
    
    print(f"✓ Bei r_min: E_sr/(mc²) = {E_sr_normalized[0]:.6f}, "
          f"E_gr/(mc²) = {E_gr_normalized[0]:.6f}")
    print(f"✓ Bei r_max: E_sr/(mc²) = {E_sr_normalized[-1]:.6e}, "
          f"E_gr/(mc²) = {E_gr_normalized[-1]:.6e}")
    print(f"  (relative Energien nehmen mit r ab)")
    print()
    
    # =========================================================================
    # Tabellen-Ausgabe (erste und letzte 5 Segmente)
    # =========================================================================
    
    print("ERGEBNISSE (erste 5 und letzte 5 Segmente):")
    print("-" * 80)
    
    # Header
    print(f"{'n':>4} {'r_n [km]':>15} {'r_n/r_s':>10} "
          f"{'γ_sr':>10} {'γ_gr':>10} "
          f"{'E_sr/(mc²)':>12} {'E_gr/(mc²)':>12} {'E_tot/(mc²)':>12}")
    print("-" * 80)
    
    # Normalisierte Energien
    E_tot_normalized = (E_tot / (m * c**2)).decompose().value
    
    # Erste 5 Segmente
    for i in range(min(5, N)):
        r_n_km = r_array[i].to(u.km).value
        r_n_rs = (r_array[i] / r_s).decompose().value
        print(f"{i:4d} {r_n_km:15.3f} {r_n_rs:10.4f} "
              f"{gamma_sr_vals[i]:10.6f} {gamma_gr_vals[i]:10.6f} "
              f"{E_sr_normalized[i]:12.6e} {E_gr_normalized[i]:12.6e} "
              f"{E_tot_normalized[i]:12.6f}")
    
    if N > 10:
        print("...")
        
        # Letzte 5 Segmente
        for i in range(max(N-5, 5), N):
            r_n_km = r_array[i].to(u.km).value
            r_n_rs = (r_array[i] / r_s).decompose().value
            print(f"{i:4d} {r_n_km:15.3f} {r_n_rs:10.4f} "
                  f"{gamma_sr_vals[i]:10.6f} {gamma_gr_vals[i]:10.6f} "
                  f"{E_sr_normalized[i]:12.6e} {E_gr_normalized[i]:12.6e} "
                  f"{E_tot_normalized[i]:12.6f}")
    
    print()
    
    # =========================================================================
    # Plots erstellen
    # =========================================================================
    
    print("VISUALISIERUNG:")
    print("-" * 80)
    
    # Normalisierte Radien (r/r_min)
    r_normalized = (r_array / r_min).decompose().value
    
    # Figure mit 3 Subplots
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Plot 1: Lorentzfaktoren
    ax = axes[0]
    ax.plot(r_normalized, gamma_sr_vals, 'b-', label='γ_sr (SR)', linewidth=2)
    ax.plot(r_normalized, gamma_gr_vals, 'r-', label='γ_gr (GR)', linewidth=2)
    ax.set_xlabel('r / r_min')
    ax.set_ylabel('Lorentzfaktor γ')
    ax.set_title('Lorentzfaktoren vs. Radius')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # Plot 2: Normalisierte Energien (einzeln)
    ax = axes[1]
    ax.plot(r_normalized, E_sr_normalized, 'b-', 
            label='E_sr / (mc²)', linewidth=2)
    ax.plot(r_normalized, E_gr_normalized, 'r-', 
            label='E_gr / (mc²)', linewidth=2)
    ax.set_xlabel('r / r_min')
    ax.set_ylabel('E / (mc²)')
    ax.set_title('Normalisierte Energie-Anteile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Plot 3: Gesamtenergie
    ax = axes[2]
    ax.plot(r_normalized, E_tot_normalized, 'g-', 
            label='E_tot / (mc²)', linewidth=2)
    ax.axhline(y=1.0, color='k', linestyle='--', 
               linewidth=1, label='E_rest / (mc²) = 1')
    ax.set_xlabel('r / r_min')
    ax.set_ylabel('E_tot / (mc²)')
    ax.set_title('Normalisierte Gesamtenergie')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    plt.tight_layout()
    
    # Speichern
    output_file = 'segmented_spacetime_energy.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Plot gespeichert: {output_file}")
    
    plt.show()
    
    # =========================================================================
    # Zusammenfassung
    # =========================================================================
    
    print()
    print("=" * 80)
    print("ZUSAMMENFASSUNG")
    print("=" * 80)
    print(f"✓ {N} Segmente zwischen {r_min/r_s:.2f} r_s und {r_max/r_s:.2f} r_s")
    print(f"✓ SR-Effekt (γ_sr - 1) bei r_min: {gamma_sr_vals[0] - 1:.6f}")
    print(f"✓ GR-Effekt (γ_gr - 1) bei r_min: {gamma_gr_vals[0] - 1:.6f}")
    print(f"✓ Maximale Gesamtenergie: {E_tot_normalized.max():.6f} mc²")
    print(f"✓ Minimale Gesamtenergie: {E_tot_normalized.min():.6f} mc²")
    print()
    print("PHYSIKALISCHE INTERPRETATION:")
    print("  • γ_sr beschreibt relativistische Effekte durch Bahngeschwindigkeit")
    print("  • γ_gr beschreibt Zeitdilatation im Gravitationsfeld")
    print("  • Beide Effekte werden mit zunehmendem Radius schwächer")
    print("  • Bei r_min (nahe Schwarzschild-Radius) sind beide signifikant")
    print("  • Bei großen Radien dominiert die Ruheenergie mc²")
    print("=" * 80)
