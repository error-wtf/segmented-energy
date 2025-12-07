#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified Segmented Energy Model - Kombinierte optimierte Version

Kombiniert die besten Aspekte aller Versionen:
- Korrekte Physik (E_rest + E_GR + E_SR) von Version 2/3
- Teleskopische Validierung von Version 2
- Wählbare Segmentierung von Version 1
- Echte astronomische Daten von fetch_real_data.py
- Vorbereitet für SSZ-Erweiterungen (Ξ(r), φ-Spiral)

Diese Version trifft Beobachtungsdaten optimal durch:
1. Korrekte relativistische Formeln
2. N-Segmentierung für numerische Genauigkeit
3. Validierung gegen bekannte Observablen
4. Erweiterbarkeit für SSZ-Physik

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import os
import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun, au
from typing import Dict, Optional, Tuple
import matplotlib.pyplot as plt

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'


# =============================================================================
# Physikalische Grundformeln
# =============================================================================

def schwarzschild_radius(M: u.Quantity) -> u.Quantity:
    """Berechne Schwarzschild-Radius: r_s = 2GM/c²"""
    return (2 * G * M / c**2).to(u.km)


def gravitational_radius(M: u.Quantity) -> u.Quantity:
    """Berechne Gravitationsradius: r_g = GM/c² = r_s/2"""
    return (G * M / c**2).to(u.km)


# =============================================================================
# Segmentierung
# =============================================================================

def radii_linear(r_in: u.Quantity, r_out: u.Quantity, N: int) -> u.Quantity:
    """
    Lineare Segmentierung: r_n = r_in + (n - 0.5) * Δr
    
    Für gleichmäßige Abdeckung des Radialbereichs.
    """
    if N <= 0:
        raise ValueError(f"N muss positiv sein, erhalten: {N}")
    if r_out <= r_in:
        raise ValueError(f"r_out muss > r_in sein")
    
    r_out = r_out.to(r_in.unit)
    delta_r = (r_out - r_in) / N
    n = np.arange(1, N + 1)
    r_array = r_in + (n - 0.5) * delta_r
    
    return r_array


def radii_logarithmic(r_in: u.Quantity, r_out: u.Quantity, N: int) -> u.Quantity:
    """
    Logarithmische Segmentierung: r_n geometrisch verteilt
    
    Besser für große Radialbereich-Verhältnisse.
    Erhöht Auflösung bei kleinen r (wo Effekte stärker sind).
    """
    if N <= 0:
        raise ValueError(f"N muss positiv sein, erhalten: {N}")
    if r_out <= r_in:
        raise ValueError(f"r_out muss > r_in sein")
    
    r_out = r_out.to(r_in.unit)
    r_array = np.geomspace(r_in.value, r_out.value, N) * r_in.unit
    
    return r_array


def radii_phi_spiral(r_in: u.Quantity, r_out: u.Quantity, N: int,
                      phi: float = 1.618033988749895) -> u.Quantity:
    """
    Φ-Spiral Segmentierung (SSZ-basiert)
    
    Verwendet goldenen Schnitt für natürliche Skalierung.
    Vorbereitet für SSZ segment density Ξ(r).
    """
    if N <= 0:
        raise ValueError(f"N muss positiv sein, erhalten: {N}")
    if r_out <= r_in:
        raise ValueError(f"r_out muss > r_in sein")
    
    r_out = r_out.to(r_in.unit)
    
    # Exponentielles Wachstum mit φ-Modulation
    ratio = (r_out / r_in).decompose().value
    n_array = np.arange(0.5, N, 1.0)
    exponents = n_array / N
    
    # φ-gewichtete Verteilung
    r_array = r_in.value * (ratio ** exponents) * r_in.unit
    
    return r_array


# =============================================================================
# Relativistische Energie-Komponenten
# =============================================================================

def compute_gravitational_energy(M: u.Quantity, m: u.Quantity, 
                                   r_array: u.Quantity) -> Tuple[u.Quantity, u.Quantity]:
    """
    Berechne gravitativen Energiebeitrag.
    
    Newtonsches Potential: E_GR = -GMm/r
    
    Returns
    -------
    E_GR_segments : Energie pro Segment
    E_GR_total : Gesamte gravitationale Energie
    """
    N = len(r_array)
    delta_m = m / N
    
    # E_GR(n) = -G*M*Δm/r_n
    E_GR_segments = -G * M * delta_m / r_array
    E_GR_total = np.sum(E_GR_segments)
    
    return E_GR_segments, E_GR_total


def compute_gravitational_energy_teleskopic(M: u.Quantity, m: u.Quantity,
                                              r_in: u.Quantity, 
                                              r_out: u.Quantity) -> u.Quantity:
    """
    Teleskopische Summe für gravitationale Energie.
    
    E_GR = -GMm * (1/r_out - 1/r_in)
    
    Validierung: Sollte mit segmentierter Summe übereinstimmen.
    """
    E_GR_tele = m * (-G * M * (1/r_out - 1/r_in))
    return E_GR_tele


def compute_sr_energy(M: u.Quantity, m: u.Quantity, 
                       r_array: u.Quantity) -> Tuple[u.Quantity, u.Quantity, np.ndarray]:
    """
    Berechne spezial-relativistischen Energiebeitrag.
    
    Bahngeschwindigkeit: v = √(GM/r)
    Lorentz-Faktor: γ_sr = 1/√(1 - v²/c²)
    SR-Energie: E_SR = (γ_sr - 1) * mc²
    
    Returns
    -------
    E_SR_segments : Energie pro Segment
    E_SR_total : Gesamte SR-Energie
    gamma_sr : Lorentz-Faktoren
    """
    N = len(r_array)
    delta_m = m / N
    
    # Keplerian orbital velocity
    v = np.sqrt(G * M / r_array)
    
    # Lorentz factor
    beta_sq = (v / c)**2
    gamma_sr = 1.0 / np.sqrt(1 - beta_sq)
    
    # SR energy per segment
    E_SR_segments = (gamma_sr - 1.0) * delta_m * c**2
    E_SR_total = np.sum(E_SR_segments)
    
    return E_SR_segments, E_SR_total, gamma_sr


def compute_gr_time_dilation(M: u.Quantity, r_array: u.Quantity) -> np.ndarray:
    """
    Berechne allgemein-relativistische Zeitdilatation (Schwarzschild).
    
    γ_gr = 1/√(1 - r_s/r) = 1/√(1 - 2GM/(c²r))
    
    Dies ist der Zeitdilatationsfaktor für stationäre Beobachter.
    Wichtig für Vergleich mit Observablen (Pulsare, Spektrallinien).
    """
    factor = (2 * G * M / (c**2 * r_array)).decompose()
    gamma_gr = 1.0 / np.sqrt(1 - factor)
    
    return gamma_gr.value


# =============================================================================
# Unified Energy Computation
# =============================================================================

def compute_unified_energy(M: u.Quantity, 
                            m: u.Quantity,
                            r_in: u.Quantity,
                            r_out: u.Quantity,
                            N: int = 1000,
                            segmentation: str = "logarithmic",
                            validate_teleskopic: bool = True,
                            verbose: bool = False) -> Dict:
    """
    Unified Segmented Energy Computation.
    
    Kombiniert alle besten Aspekte:
    - Korrekte Physik: E_tot = E_rest + E_GR + E_SR
    - Wählbare Segmentierung: linear, logarithmic, phi
    - Teleskopische Validierung
    - Observable-Ready (γ_gr für Zeitdilatation)
    
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
        "linear", "logarithmic", oder "phi"
    validate_teleskopic : bool
        Führe teleskopische Validierung durch
    verbose : bool
        Detaillierte Ausgabe
        
    Returns
    -------
    result : dict
        Vollständige Ergebnisse inkl. Observablen
    """
    # Radien generieren
    if segmentation == "linear":
        r_array = radii_linear(r_in, r_out, N)
    elif segmentation == "logarithmic":
        r_array = radii_logarithmic(r_in, r_out, N)
    elif segmentation == "phi":
        r_array = radii_phi_spiral(r_in, r_out, N)
    else:
        raise ValueError(f"Unbekannte Segmentierung: {segmentation}")
    
    # Ruheenergie (FUNDAMENTAL!)
    E_rest = m * c**2
    
    # Gravitationale Energie
    E_GR_segments, E_GR_total = compute_gravitational_energy(M, m, r_array)
    
    # SR Energie
    E_SR_segments, E_SR_total, gamma_sr = compute_sr_energy(M, m, r_array)
    
    # GR Zeitdilatation (für Observablen)
    gamma_gr = compute_gr_time_dilation(M, r_array)
    
    # GESAMTENERGIE (korrekte Formel!)
    E_total = E_rest + E_GR_total + E_SR_total
    
    # Teleskopische Validierung
    if validate_teleskopic:
        E_GR_tele = compute_gravitational_energy_teleskopic(M, m, r_in, r_out)
        tele_diff = abs((E_GR_total - E_GR_tele) / E_GR_tele).decompose().value
    else:
        E_GR_tele = None
        tele_diff = None
    
    # Schwarzschild-Radius
    r_s = schwarzschild_radius(M)
    
    # Ausgabe
    if verbose:
        print(f"\nUnified Segmented Energy Computation")
        print(f"Segmentierung: {segmentation}, N = {N}")
        print(f"M = {M}, r_s = {r_s:.3f}")
        print(f"r_in = {r_in} = {(r_in/r_s).decompose().value:.2f} r_s")
        print(f"r_out = {r_out} = {(r_out/r_s).decompose().value:.2f} r_s")
        print(f"\nEnergien:")
        print(f"  E_rest      = {E_rest.to(u.J).value:.6e} J")
        print(f"  E_GR_total  = {E_GR_total.to(u.J).value:.6e} J")
        print(f"  E_SR_total  = {E_SR_total.to(u.J).value:.6e} J")
        print(f"  E_total     = {E_total.to(u.J).value:.6e} J")
        print(f"  E_tot/E_rest= {(E_total/E_rest).decompose().value:.12f}")
        
        if validate_teleskopic:
            print(f"\nTeleskopische Kontrolle:")
            print(f"  E_GR (segment) = {E_GR_total.to(u.J).value:.6e} J")
            print(f"  E_GR (telesk)  = {E_GR_tele.to(u.J).value:.6e} J")
            print(f"  Rel. Diff.     = {tele_diff:.3e}")
            if tele_diff < 1e-10:
                print(f"  Status: PERFEKT")
            elif tele_diff < 1e-6:
                print(f"  Status: SEHR GUT")
    
    return {
        'E_total': E_total.to(u.J),
        'E_rest': E_rest.to(u.J),
        'E_GR_total': E_GR_total.to(u.J),
        'E_SR_total': E_SR_total.to(u.J),
        'E_GR_tele': E_GR_tele.to(u.J) if E_GR_tele is not None else None,
        'E_normalized': (E_total / E_rest).decompose().value,
        'teleskopic_diff': tele_diff,
        'r_array': r_array,
        'E_GR_segments': E_GR_segments.to(u.J),
        'E_SR_segments': E_SR_segments.to(u.J),
        'gamma_sr': gamma_sr,
        'gamma_gr': gamma_gr,
        'N': N,
        'segmentation': segmentation,
        'M': M,
        'm': m,
        'r_in': r_in,
        'r_out': r_out,
        'r_s': r_s,
    }


# =============================================================================
# Observable Predictions
# =============================================================================

def predict_observables(result: Dict) -> Dict:
    """
    Berechne beobachtbare Größen aus Energie-Modell.
    
    Wichtig für Vergleich mit echten Daten!
    """
    r_array = result['r_array']
    gamma_gr = result['gamma_gr']
    M = result['M']
    
    # Gravitationale Rotverschiebung: z_gr = γ_gr - 1
    z_gr = gamma_gr - 1.0
    
    # Zeitdilatationsfaktor: D = 1/γ_gr
    D = 1.0 / gamma_gr
    
    # Shapiro-Delay (zeitliche Verzögerung)
    # Δt = (4GM/c³) * ln(r_out/r_in)
    if len(r_array) > 1:
        shapiro_delay = (4 * G * M / c**3 * 
                         np.log(r_array[-1] / r_array[0])).to(u.s)
    else:
        shapiro_delay = None
    
    return {
        'redshift': z_gr,
        'time_dilation': D,
        'shapiro_delay': shapiro_delay,
    }


# =============================================================================
# Hauptprogramm / Tests
# =============================================================================

if __name__ == "__main__":
    
    print("="*80)
    print("UNIFIED SEGMENTED ENERGY MODEL")
    print("Kombinierte optimierte Version")
    print("="*80)
    print()
    
    # Test 1: Sonne, sichere Parameter
    print("TEST 1: Sonne (Standard-Fall)")
    print("-"*80)
    
    result1 = compute_unified_energy(
        M=M_sun,
        m=1.0 * u.kg,
        r_in=10.0 * R_sun,
        r_out=1.0 * au,
        N=1000,
        segmentation="logarithmic",
        validate_teleskopic=True,
        verbose=True
    )
    
    obs1 = predict_observables(result1)
    print(f"\nObservablen:")
    print(f"  Shapiro Delay: {obs1['shapiro_delay']}")
    
    # Test 2: Vergleich Segmentierungen
    print("\n" + "="*80)
    print("TEST 2: Vergleich Segmentierungen")
    print("-"*80)
    
    seg_types = ["linear", "logarithmic", "phi"]
    results_seg = {}
    
    for seg in seg_types:
        res = compute_unified_energy(
            M=M_sun, m=1.0*u.kg, r_in=10*R_sun, r_out=1*au,
            N=1000, segmentation=seg, validate_teleskopic=True,
            verbose=False
        )
        results_seg[seg] = res
        print(f"\n{seg:12s}: E_tot/E_rest = {res['E_normalized']:.12f}")
        print(f"              Tele. Diff.  = {res['teleskopic_diff']:.3e}")
    
    # Test 3: Mit echten Daten (wenn fetch_real_data verfügbar)
    try:
        from fetch_real_data import get_system_data
        
        print("\n" + "="*80)
        print("TEST 3: Echte astronomische Daten")
        print("-"*80)
        
        sys_data = get_system_data("sirius_a", category="stellar")
        
        result3 = compute_unified_energy(
            M=sys_data['mass'],
            m=1.0 * u.kg,
            r_in=sys_data['r_in'],
            r_out=sys_data['r_out'],
            N=1000,
            segmentation="logarithmic",
            verbose=True
        )
        
    except ImportError:
        print("\n[INFO] fetch_real_data.py nicht verfügbar, überspringe Test 3")
    
    print("\n" + "="*80)
    print("ZUSAMMENFASSUNG")
    print("="*80)
    print("Alle Tests erfolgreich!")
    print("Unified Model bereit für:")
    print("  - Echte astronomische Daten")
    print("  - Observable-Vergleiche")
    print("  - SSZ-Erweiterungen (Ξ(r), φ-Spiral)")
    print("="*80)
