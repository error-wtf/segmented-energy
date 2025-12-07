#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmented Spacetime (SSZ) Energy Model - Complete Implementation

Implementiert das vollständige SSZ-Modell mit:
- Segment Density: Xi(r) = Xi_max(1 - exp(-phi*r/r_s))
- SSZ Time Dilation: D_SSZ = 1/(1 + Xi(r))
- phi-Spiral Segmentierung
- Alle SSZ-spezifischen Observable

Basiert auf der SSZ-Theorie von Casu & Wrede.

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import os
import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun, au
from typing import Dict, Tuple
import matplotlib.pyplot as plt

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'


# =============================================================================
# SSZ Fundamental Constants
# =============================================================================

PHI = (1 + np.sqrt(5)) / 2  # Golden Ratio phi = 1.618...

# SSZ Parameter (aus 01_MATHEMATICAL_FOUNDATIONS.md)
XI_MAX_DEFAULT = 0.8        # Maximale Segment Density
ALPHA_FS = 1/137.035999      # Feinstrukturkonstante


# =============================================================================
# SSZ Core Functions
# =============================================================================

def schwarzschild_radius(M: u.Quantity) -> u.Quantity:
    """Schwarzschild-Radius: r_s = 2GM/c²"""
    return (2 * G * M / c**2).to(u.km)


def segment_density_Xi(r: u.Quantity, M: u.Quantity, 
                        Xi_max: float = XI_MAX_DEFAULT) -> np.ndarray:
    """
    SSZ Segment Density (Exponential Form).
    
    Xi(r) = Xi_max * (1 - exp(-phi * r/r_s))
    
    Parameters
    ----------
    r : astropy.Quantity
        Radius oder Array von Radien
    M : astropy.Quantity
        Zentrale Masse
    Xi_max : float
        Maximale Segment Density (default: 0.8)
        
    Returns
    -------
    Xi : numpy.ndarray
        Segment Density (dimensionslos)
        
    Notes
    -----
    - Xi(r) -> 0 fuer r -> inf (Kontinuum)
    - Xi(r) -> Xi_max fuer r -> r_s (maximal diskret)
    - Saturation verhindert Singularitaet
    """
    r_s = schwarzschild_radius(M)
    
    # Xi(r) = Xi_max * (1 - exp(-phi * r_s/r))  [KORREKTE FORM!]
    # Bei r -> inf: Xi -> 0 (Kontinuum)
    # Bei r -> r_s: Xi -> Xi_max (maximal diskret)
    rs_over_r = (r_s / r).decompose().value
    Xi = Xi_max * (1 - np.exp(-PHI * rs_over_r))
    
    return Xi


def ssz_time_dilation(r: u.Quantity, M: u.Quantity,
                       Xi_max: float = XI_MAX_DEFAULT) -> np.ndarray:
    """
    SSZ Time Dilation Factor.
    
    D_SSZ(r) = 1 / (1 + Xi(r))
    
    Parameters
    ----------
    r : astropy.Quantity
        Radius
    M : astropy.Quantity
        Zentrale Masse
    Xi_max : float
        Maximale Segment Density
        
    Returns
    -------
    D_SSZ : numpy.ndarray
        SSZ Zeitdilatationsfaktor (dimensionslos)
        
    Notes
    -----
    - D_SSZ < 1: Zeit läuft langsamer
    - D_SSZ = 1: Flacher Raum (Xi = 0)
    - Singularitaetsfrei: D_SSZ(r_s) = 1/(1 + Xi_max) ~ 0.56
    """
    Xi = segment_density_Xi(r, M, Xi_max)
    D_SSZ = 1.0 / (1 + Xi)
    
    return D_SSZ


def gr_time_dilation(r: u.Quantity, M: u.Quantity) -> np.ndarray:
    """
    GR Time Dilation (Schwarzschild).
    
    D_GR(r) = sqrt(1 - r_s/r)
    
    Parameters
    ----------
    r : astropy.Quantity
        Radius
    M : astropy.Quantity
        Zentrale Masse
        
    Returns
    -------
    D_GR : numpy.ndarray
        GR Zeitdilatationsfaktor
    """
    r_s = schwarzschild_radius(M)
    factor = (r_s / r).decompose().value
    D_GR = np.sqrt(1 - factor)
    
    return D_GR


def universal_intersection(M: u.Quantity) -> u.Quantity:
    """
    Universal Intersection Point: r* wo D_GR = D_SSZ.
    
    r* = 1.386562 * r_s (massenunabhängig!)
    
    Parameters
    ----------
    M : astropy.Quantity
        Zentrale Masse
        
    Returns
    -------
    r_star : astropy.Quantity
        Universal Intersection Radius
    """
    r_s = schwarzschild_radius(M)
    r_star = 1.386562 * r_s
    
    return r_star


# =============================================================================
# phi-Spiral Segmentierung
# =============================================================================

def radii_phi_spiral(r_in: u.Quantity, r_out: u.Quantity, N: int,
                      phi: float = PHI) -> u.Quantity:
    """
    phi-Spiral Radial Segmentierung (SSZ-optimiert).
    
    Verwendet phi-gewichtete exponentielle Verteilung fuer
    natuerliche Skalierung konsistent mit SSZ-Theorie.
    
    Parameters
    ----------
    r_in : astropy.Quantity
        Innerer Radius
    r_out : astropy.Quantity
        Äußerer Radius
    N : int
        Anzahl Segmente
    phi : float
        Golden ratio (default: phi)
        
    Returns
    -------
    r_array : astropy.Quantity
        Array von N Radien
    """
    if N <= 0:
        raise ValueError(f"N muss positiv sein, erhalten: {N}")
    if r_out <= r_in:
        raise ValueError(f"r_out muss > r_in sein")
    
    r_out = r_out.to(r_in.unit)
    
    # phi-gewichtete exponentielle Verteilung
    ratio = (r_out / r_in).decompose().value
    n_array = np.arange(0.5, N, 1.0)
    
    # Exponent moduliert mit phi
    exponents = (n_array / N) ** (1 / phi)
    
    r_array = r_in.value * (ratio ** exponents) * r_in.unit
    
    return r_array


def radii_logarithmic(r_in: u.Quantity, r_out: u.Quantity, N: int) -> u.Quantity:
    """Logarithmische Segmentierung (für Vergleich)"""
    r_out = r_out.to(r_in.unit)
    r_array = np.geomspace(r_in.value, r_out.value, N) * r_in.unit
    return r_array


# =============================================================================
# SSZ Energy Components
# =============================================================================

def compute_ssz_energy(M: u.Quantity, m: u.Quantity, r_array: u.Quantity,
                        Xi_max: float = XI_MAX_DEFAULT) -> Dict:
    """
    Berechne SSZ-Energien mit Segment Density Korrektur.
    
    SSZ modifiziert die Energie durch:
    1. Zeitdilatation: D_SSZ = 1/(1 + Xi(r))
    2. Modifizierte Lorentzfaktoren
    3. Segment-basierte Korrekturen
    
    Parameters
    ----------
    M : astropy.Quantity
        Zentrale Masse
    m : astropy.Quantity
        Testmasse
    r_array : astropy.Quantity
        Array von Radien
    Xi_max : float
        Maximale Segment Density
        
    Returns
    -------
    result : dict
        Vollständige SSZ-Ergebnisse
    """
    N = len(r_array)
    delta_m = m / N
    
    # SSZ Segment Density
    Xi = segment_density_Xi(r_array, M, Xi_max)
    
    # SSZ Time Dilation
    D_SSZ = ssz_time_dilation(r_array, M, Xi_max)
    
    # GR Time Dilation (für Vergleich)
    D_GR = gr_time_dilation(r_array, M)
    
    # Ruheenergie
    E_rest = m * c**2
    
    # Keplerian velocity
    v = np.sqrt(G * M / r_array)
    beta_sq = (v / c)**2
    
    # SSZ-modifizierte Lorentzfaktoren
    # gamma_SSZ = gamma_SR * (1 + SSZ-Korrektur)
    gamma_sr = 1.0 / np.sqrt(1 - beta_sq)
    
    # SSZ-Korrektur: Zeit laeuft langsamer -> effektiv hoeheres gamma
    gamma_ssz = gamma_sr / D_SSZ
    
    # Gravitationale Energie (unverändert in SSZ)
    E_GR_segments = -G * M * delta_m / r_array
    E_GR_total = np.sum(E_GR_segments)
    
    # SSZ-modifizierte SR-Energie
    # E_SR_SSZ = (gamma_SSZ - 1) * dm * c^2
    E_SR_SSZ_segments = (gamma_ssz - 1.0) * delta_m * c**2
    E_SR_SSZ_total = np.sum(E_SR_SSZ_segments)
    
    # Standard SR (für Vergleich)
    E_SR_std_segments = (gamma_sr - 1.0) * delta_m * c**2
    E_SR_std_total = np.sum(E_SR_std_segments)
    
    # Gesamtenergie
    E_total_SSZ = E_rest + E_GR_total + E_SR_SSZ_total
    E_total_GR = E_rest + E_GR_total + E_SR_std_total
    
    return {
        'E_total_SSZ': E_total_SSZ.to(u.J),
        'E_total_GR': E_total_GR.to(u.J),
        'E_rest': E_rest.to(u.J),
        'E_GR_total': E_GR_total.to(u.J),
        'E_SR_SSZ_total': E_SR_SSZ_total.to(u.J),
        'E_SR_std_total': E_SR_std_total.to(u.J),
        'E_normalized_SSZ': (E_total_SSZ / E_rest).decompose().value,
        'E_normalized_GR': (E_total_GR / E_rest).decompose().value,
        'Xi': Xi,
        'D_SSZ': D_SSZ,
        'D_GR': D_GR,
        'gamma_ssz': gamma_ssz,
        'gamma_sr': gamma_sr,
        'r_array': r_array,
        'N': N,
        'M': M,
        'm': m,
        'Xi_max': Xi_max,
    }


# =============================================================================
# SSZ Observable
# =============================================================================

def predict_ssz_observables(result: Dict) -> Dict:
    """
    Berechne SSZ-spezifische Observable.
    
    Parameters
    ----------
    result : dict
        Output von compute_ssz_energy
        
    Returns
    -------
    obs : dict
        SSZ Observable
    """
    # SSZ Redshift: z_SSZ = 1/D_SSZ - 1
    z_SSZ = 1.0 / result['D_SSZ'] - 1
    
    # GR Redshift (für Vergleich)
    z_GR = 1.0 / result['D_GR'] - 1
    
    # Shapiro Delay (SSZ-modifiziert)
    r_array = result['r_array']
    M = result['M']
    
    if len(r_array) > 1:
        # SSZ: Delay ist groesser wegen Xi(r)
        shapiro_base = (4 * G * M / c**3 * np.log(r_array[-1] / r_array[0])).to(u.s)
        
        # SSZ-Korrektur: Mittelwert von 1 + Xi(r)
        Xi_mean = np.mean(result['Xi'])
        shapiro_SSZ = shapiro_base * (1 + Xi_mean)
        shapiro_GR = shapiro_base
    else:
        shapiro_SSZ = None
        shapiro_GR = None
    
    return {
        'z_SSZ': z_SSZ,
        'z_GR': z_GR,
        'D_SSZ': result['D_SSZ'],
        'D_GR': result['D_GR'],
        'shapiro_SSZ': shapiro_SSZ,
        'shapiro_GR': shapiro_GR,
        'Xi': result['Xi'],
    }


# =============================================================================
# Unified SSZ Computation
# =============================================================================

def compute_ssz_unified(M: u.Quantity, m: u.Quantity,
                         r_in: u.Quantity, r_out: u.Quantity,
                         N: int = 1000,
                         segmentation: str = "phi",
                         Xi_max: float = XI_MAX_DEFAULT,
                         verbose: bool = False) -> Dict:
    """
    Unified SSZ Energy Computation.
    
    Parameters
    ----------
    M : astropy.Quantity
        Zentrale Masse
    m : astropy.Quantity
        Testmasse
    r_in : astropy.Quantity
        Innerer Radius
    r_out : astropy.Quantity
        Äußerer Radius
    N : int
        Anzahl Segmente
    segmentation : str
        "phi" (SSZ-optimiert) oder "logarithmic"
    Xi_max : float
        Maximale Segment Density
    verbose : bool
        Detaillierte Ausgabe
        
    Returns
    -------
    result : dict
        Vollständige SSZ-Ergebnisse inkl. GR-Vergleich
    """
    # Radien generieren
    if segmentation == "phi":
        r_array = radii_phi_spiral(r_in, r_out, N)
    elif segmentation == "logarithmic":
        r_array = radii_logarithmic(r_in, r_out, N)
    else:
        raise ValueError(f"Unbekannte Segmentierung: {segmentation}")
    
    # SSZ-Energie berechnen
    result = compute_ssz_energy(M, m, r_array, Xi_max)
    result['segmentation'] = segmentation
    result['r_in'] = r_in
    result['r_out'] = r_out
    
    # Schwarzschild-Radius
    r_s = schwarzschild_radius(M)
    result['r_s'] = r_s
    
    # Universal Intersection
    r_star = universal_intersection(M)
    result['r_star'] = r_star
    
    # Observable
    obs = predict_ssz_observables(result)
    result['observables'] = obs
    
    if verbose:
        print(f"\nSSZ Unified Energy Computation")
        print(f"Segmentierung: {segmentation}, N = {N}, Xi_max = {Xi_max}")
        print(f"M = {M}, r_s = {r_s:.3f}")
        print(f"r* = {r_star:.3f} (Universal Intersection)")
        print(f"r_in = {r_in} = {(r_in/r_s).decompose().value:.2f} r_s")
        print(f"r_out = {r_out} = {(r_out/r_s).decompose().value:.2f} r_s")
        print(f"\nSSZ vs GR:")
        print(f"  E_total_SSZ    = {result['E_total_SSZ'].value:.6e} J")
        print(f"  E_total_GR     = {result['E_total_GR'].value:.6e} J")
        print(f"  SSZ/GR Ratio   = {(result['E_total_SSZ']/result['E_total_GR']).decompose().value:.6f}")
        print(f"  E_SSZ/E_rest   = {result['E_normalized_SSZ']:.12f}")
        print(f"  E_GR/E_rest    = {result['E_normalized_GR']:.12f}")
        print(f"\nSegment Density:")
        print(f"  Xi(r_in)  = {result['Xi'][0]:.6f}")
        print(f"  Xi(r_out) = {result['Xi'][-1]:.6f}")
        print(f"  Xi_mean   = {np.mean(result['Xi']):.6f}")
    
    return result


# =============================================================================
# Main / Tests
# =============================================================================

if __name__ == "__main__":
    
    print("="*80)
    print("SEGMENTED SPACETIME (SSZ) ENERGY MODEL")
    print("Vollstaendige Implementation mit Xi(r) und phi-Spiral")
    print("="*80)
    print()
    
    # Test 1: Sonne (Standard-Fall)
    print("TEST 1: Sonne (SSZ vs GR)")
    print("-"*80)
    
    result = compute_ssz_unified(
        M=M_sun,
        m=1.0 * u.kg,
        r_in=10.0 * R_sun,
        r_out=1.0 * au,
        N=1000,
        segmentation="phi",
        Xi_max=0.8,
        verbose=True
    )
    
    obs = result['observables']
    
    print(f"\nObservable (bei r_in):")
    print(f"  z_SSZ = {obs['z_SSZ'][0]:.6e}")
    print(f"  z_GR  = {obs['z_GR'][0]:.6e}")
    print(f"  D_SSZ = {obs['D_SSZ'][0]:.6f}")
    print(f"  D_GR  = {obs['D_GR'][0]:.6f}")
    
    if obs['shapiro_SSZ'] is not None:
        print(f"\nShapiro Delay:")
        print(f"  SSZ: {obs['shapiro_SSZ'].to(u.us).value:.3f} us")
        print(f"  GR:  {obs['shapiro_GR'].to(u.us).value:.3f} us")
        print(f"  SSZ/GR: {(obs['shapiro_SSZ']/obs['shapiro_GR']).decompose().value:.3f}")
    
    print("\n" + "="*80)
    print("SSZ-Modell bereit für vollständige Tests!")
    print("="*80)
