#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmentiertes Energie-Modell mit SR und GR

Implementiert ein vollständiges Energie-Modell mit:
- Ruheenergie E_rest = m * c²
- Spezieller Relativität (SR): E_SR = (γ_SR - 1) * m * c²
- Allgemeiner Relativität (GR): E_GR = (γ_GR - 1) * m * c²

Gesamtenergie: E_tot = E_rest + E_GR + E_SR = Σ_n [E_rest(n) + E_GR(n) + E_SR(n)]

Mathematisches Modell:
======================

Konstanten:
- c: Lichtgeschwindigkeit
- G: Gravitationskonstante

Segmentierung:
- N = Anzahl der Segmente
- Für jedes Segment n = 1..N:
  - m_n: Masse im Segment
  - v_n: lokale Geschwindigkeit (SR)
  - M_n: eingeschlossene gravitative Masse
  - r_n: Radius des Segments vom Zentrum

Lorentz-Faktoren:
- γ_SR(n) = 1 / sqrt(1 - v_n²/c²)
- r_s(n) = 2 * G * M_n / c²
- γ_GR(n) = 1 / sqrt(1 - r_s(n)/r_n) = 1 / sqrt(1 - 2*G*M_n/(r_n*c²))

Energien pro Segment:
- E_rest(n) = m_n * c²
- E_SR(n) = (γ_SR(n) - 1) * m_n * c²
- E_GR(n) = (γ_GR(n) - 1) * m_n * c²

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import G, c, M_sun, R_sun
from astropy import units as u

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'


# =============================================================================
# Lorentz-Faktoren
# =============================================================================

def gamma_sr(v):
    """
    Special Relativistic Lorentz-Faktor.
    
    γ_SR = 1 / sqrt(1 - v²/c²)
    
    Parameters
    ----------
    v : astropy.Quantity
        Geschwindigkeit (mit Einheit)
        
    Returns
    -------
    gamma_SR : numpy.ndarray
        Lorentz-Faktor (dimensionslos)
    """
    beta_sq = (v / c)**2
    gamma = 1.0 / np.sqrt(1 - beta_sq)
    return gamma.value


def gamma_gr(M, r):
    """
    General Relativistic Gamma-Faktor (Schwarzschild).
    
    r_s = 2 * G * M / c²
    γ_GR = 1 / sqrt(1 - r_s/r) = 1 / sqrt(1 - 2*G*M/(r*c²))
    
    Parameters
    ----------
    M : astropy.Quantity
        Zentrale Masse (mit Einheit)
    r : astropy.Quantity
        Radius (mit Einheit)
        
    Returns
    -------
    gamma_GR : numpy.ndarray
        GR Gamma-Faktor (dimensionslos)
    """
    r_s = 2 * G * M / c**2
    factor = (r_s / r).decompose()
    gamma = 1.0 / np.sqrt(1 - factor)
    return gamma.value


# =============================================================================
# Segment-Energien
# =============================================================================

def segment_energies(m_segments, M_segments, r_segments, v_segments):
    """
    Berechne Energien für alle Segmente.
    
    Parameters
    ----------
    m_segments : astropy.Quantity (array)
        Masse pro Segment
    M_segments : astropy.Quantity (array)
        Eingeschlossene Masse pro Segment
    r_segments : astropy.Quantity (array)
        Radius pro Segment
    v_segments : astropy.Quantity (array)
        Geschwindigkeit pro Segment
        
    Returns
    -------
    result : dict
        Dictionary mit allen Segment-Energien und Summen
    """
    N = len(m_segments)
    
    # Lorentz-Faktoren für alle Segmente
    gamma_SR_segments = gamma_sr(v_segments)
    gamma_GR_segments = gamma_gr(M_segments, r_segments)
    
    # Energien pro Segment
    E_rest_segments = m_segments * c**2
    E_SR_segments = (gamma_SR_segments - 1.0) * m_segments * c**2
    E_GR_segments = (gamma_GR_segments - 1.0) * m_segments * c**2
    
    # Summen über alle Segmente
    E_rest_total = np.sum(E_rest_segments)
    E_SR_total = np.sum(E_SR_segments)
    E_GR_total = np.sum(E_GR_segments)
    
    # Gesamtenergie
    E_tot = E_rest_total + E_GR_total + E_SR_total
    
    return {
        'E_rest_segments': E_rest_segments,
        'E_SR_segments': E_SR_segments,
        'E_GR_segments': E_GR_segments,
        'gamma_SR_segments': gamma_SR_segments,
        'gamma_GR_segments': gamma_GR_segments,
        'E_rest_total': E_rest_total,
        'E_SR_total': E_SR_total,
        'E_GR_total': E_GR_total,
        'E_tot': E_tot,
        'N': N,
    }


def total_energy_from_segments(m_segments, M_segments, r_segments, v_segments):
    """
    Wrapper für segment_energies mit kompakter Ausgabe.
    
    Returns
    -------
    result : dict
        Kompaktes Dictionary mit allen Ergebnissen
    """
    return segment_energies(m_segments, M_segments, r_segments, v_segments)


# =============================================================================
# Segmentierungs-Schema
# =============================================================================

def create_segments(M_total, R_obj, N=100, r_min_factor=3.0):
    """
    Erstelle einfaches Segmentierungs-Schema.
    
    EINFACHES MODELL: Kann später durch komplexere Schemas ersetzt werden.
    
    Annahmen:
    - Gleichmäßige Massenverteilung: m_n = M_total / N
    - Kepler-Geschwindigkeiten: v_n = sqrt(G * M_total / r_n)
    - Eingeschlossene Masse: M_n = M_total (Näherung)
    - Radien: linear zwischen r_min und r_max
    
    Parameters
    ----------
    M_total : astropy.Quantity
        Gesamtmasse des Objekts
    R_obj : astropy.Quantity
        Radius des Objekts
    N : int
        Anzahl Segmente
    r_min_factor : float
        Faktor für minimalen Radius (in Schwarzschild-Radien)
        
    Returns
    -------
    segments : dict
        Dictionary mit m_n, M_n, r_n, v_n
    """
    # Schwarzschild-Radius
    r_s = 2 * G * M_total / c**2
    
    # Radien: linear zwischen r_min und R_obj
    r_min = r_min_factor * r_s
    r_max = R_obj
    r_segments = np.linspace(r_min.value, r_max.to(r_min.unit).value, N) * r_min.unit
    
    # Massen: gleichmäßig verteilt
    m_segments = np.full(N, (M_total / N).value) * M_total.unit
    
    # Eingeschlossene Masse: vereinfacht M_total für alle
    M_segments = np.full(N, M_total.value) * M_total.unit
    
    # Kepler-Geschwindigkeiten
    v_segments = np.sqrt(G * M_total / r_segments)
    
    return {
        'm_segments': m_segments,
        'M_segments': M_segments,
        'r_segments': r_segments,
        'v_segments': v_segments,
        'r_s': r_s,
        'N': N,
    }


# =============================================================================
# Astrophysikalische Test-Objekte
# =============================================================================

def test_objects():
    """
    Definiere realistische astrophysikalische Test-Objekte.
    
    Returns
    -------
    objects : dict
        Dictionary mit Objekten und ihren Parametern
    """
    return {
        'main_sequence': {
            'name': 'Hauptreihenstern (Sonne)',
            'M': 1.0 * M_sun,
            'R': 1.0 * R_sun,
            'category': 'main_sequence',
        },
        'white_dwarf': {
            'name': 'Weißer Zwerg',
            'M': 0.6 * M_sun,
            'R': 0.012 * R_sun,
            'category': 'white_dwarf',
        },
        'neutron_star': {
            'name': 'Neutronenstern',
            'M': 1.4 * M_sun,
            'R': 12.0 * u.km,
            'category': 'neutron_star',
        },
    }


# =============================================================================
# Berechnungen und Analyse
# =============================================================================

def analyze_object(obj_data, N=100):
    """
    Analysiere ein astrophysikalisches Objekt.
    
    Parameters
    ----------
    obj_data : dict
        Objekt-Daten (M, R, name, category)
    N : int
        Anzahl Segmente
        
    Returns
    -------
    result : dict
        Vollständige Analyse-Ergebnisse
    """
    M = obj_data['M']
    R = obj_data['R']
    name = obj_data['name']
    category = obj_data['category']
    
    # Erstelle Segmente
    segments = create_segments(M, R, N=N)
    
    # Berechne Energien
    energies = total_energy_from_segments(
        segments['m_segments'],
        segments['M_segments'],
        segments['r_segments'],
        segments['v_segments']
    )
    
    # Normalisierungen
    E_GR_norm = (energies['E_GR_total'] / energies['E_rest_total']).decompose().value
    E_SR_norm = (energies['E_SR_total'] / energies['E_rest_total']).decompose().value
    E_tot_norm = (energies['E_tot'] / energies['E_rest_total']).decompose().value
    
    # r/r_s Verhältnis
    r_over_rs = (R / segments['r_s']).decompose().value
    
    return {
        'name': name,
        'category': category,
        'M': M,
        'R': R,
        'r_s': segments['r_s'],
        'r_over_rs': r_over_rs,
        'N': N,
        'E_rest_total': energies['E_rest_total'],
        'E_GR_total': energies['E_GR_total'],
        'E_SR_total': energies['E_SR_total'],
        'E_tot': energies['E_tot'],
        'E_GR_norm': E_GR_norm,
        'E_SR_norm': E_SR_norm,
        'E_tot_norm': E_tot_norm,
        'segments': segments,
        'energies': energies,
    }


# =============================================================================
# Ausgabe und Tabellen
# =============================================================================

def print_results(results):
    """Drucke Ergebnisse als Tabelle."""
    
    print("="*80)
    print("SEGMENTIERTES ENERGIE-MODELL - ERGEBNISSE")
    print("="*80)
    print()
    
    for res in results:
        print(f"\n{res['name']}")
        print("-"*80)
        print(f"  Masse:              {res['M'].to(M_sun):.4f}")
        print(f"  Radius:             {res['R'].to(u.km):.2f}")
        print(f"  Schwarzschild r_s:  {res['r_s'].to(u.km):.6f}")
        print(f"  R / r_s:            {res['r_over_rs']:.2e}")
        print(f"  Segmente:           {res['N']}")
        print()
        print(f"  E_rest_total:       {res['E_rest_total'].to(u.J):.6e}")
        print(f"  E_GR_total:         {res['E_GR_total'].to(u.J):.6e}")
        print(f"  E_SR_total:         {res['E_SR_total'].to(u.J):.6e}")
        print(f"  E_tot:              {res['E_tot'].to(u.J):.6e}")
        print()
        print(f"  E_GR / E_rest:      {res['E_GR_norm']:.6e}")
        print(f"  E_SR / E_rest:      {res['E_SR_norm']:.6e}")
        print(f"  E_tot / E_rest:     {res['E_tot_norm']:.12f}")
    
    print("\n" + "="*80)


# =============================================================================
# Visualisierungen
# =============================================================================

def plot_energy_normalization(results):
    """Plot 1: Energie-Normalisierung vs. Masse."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Daten extrahieren
    masses = [res['M'].to(M_sun).value for res in results]
    E_GR_norms = [res['E_GR_norm'] for res in results]
    E_SR_norms = [res['E_SR_norm'] for res in results]
    E_tot_norms = [res['E_tot_norm'] for res in results]
    categories = [res['category'] for res in results]
    names = [res['name'] for res in results]
    
    # Farben pro Kategorie
    colors = {
        'main_sequence': 'blue',
        'white_dwarf': 'orange',
        'neutron_star': 'red',
    }
    
    # Plot 1: GR und SR Beiträge
    for i, cat in enumerate(categories):
        ax1.scatter(masses[i], abs(E_GR_norms[i]), 
                   color=colors[cat], marker='o', s=200, 
                   label=f'{names[i]} (GR)' if i < 3 else None)
        ax1.scatter(masses[i], E_SR_norms[i], 
                   color=colors[cat], marker='s', s=200, alpha=0.6,
                   label=f'{names[i]} (SR)' if i < 3 else None)
    
    ax1.set_xlabel('Masse [M_sun]', fontsize=12)
    ax1.set_ylabel('|E_GR| / E_rest, E_SR / E_rest', fontsize=12)
    ax1.set_title('Relativistische Energiebeiträge', fontsize=14)
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    
    # Plot 2: Totale Energie
    for i, cat in enumerate(categories):
        ax2.scatter(masses[i], E_tot_norms[i], 
                   color=colors[cat], s=200, label=names[i])
    
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='E_rest')
    ax2.set_xlabel('Masse [M_sun]', fontsize=12)
    ax2.set_ylabel('E_tot / E_rest', fontsize=12)
    ax2.set_title('Totale Energie (normalisiert)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    
    plt.tight_layout()
    plt.savefig('energy_normalization.png', dpi=150, bbox_inches='tight')
    print("Plot gespeichert: energy_normalization.png")
    plt.show()


def plot_gamma_factors(result):
    """Plot 2: Gamma-Faktoren für ein Objekt."""
    
    segments = result['segments']
    energies = result['energies']
    
    # r/r_s für alle Segmente
    r_over_rs = (segments['r_segments'] / segments['r_s']).decompose().value
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Gamma-Faktoren
    ax1.plot(r_over_rs, energies['gamma_GR_segments'], 
            'r-', linewidth=2, label='γ_GR')
    ax1.plot(r_over_rs, energies['gamma_SR_segments'], 
            'b-', linewidth=2, label='γ_SR')
    ax1.set_xlabel('r / r_s', fontsize=12)
    ax1.set_ylabel('Gamma-Faktor', fontsize=12)
    ax1.set_title(f'Lorentz-Faktoren: {result["name"]}', fontsize=14)
    ax1.set_xscale('log')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Energien pro Segment
    ax2.plot(r_over_rs, energies['E_GR_segments'].to(u.J).value, 
            'r-', linewidth=2, label='E_GR(n)')
    ax2.plot(r_over_rs, energies['E_SR_segments'].to(u.J).value, 
            'b-', linewidth=2, label='E_SR(n)')
    ax2.set_xlabel('r / r_s', fontsize=12)
    ax2.set_ylabel('Energie pro Segment [J]', fontsize=12)
    ax2.set_title('Energie-Beiträge pro Segment', fontsize=14)
    ax2.set_xscale('log')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'gamma_factors_{result["category"]}.png', dpi=150, bbox_inches='tight')
    print(f"Plot gespeichert: gamma_factors_{result['category']}.png")
    plt.show()


def plot_segment_energies(result):
    """Plot 3: Balkendiagramm der Segment-Energien."""
    
    energies = result['energies']
    N = result['N']
    
    # Nur erste 20 Segmente für Übersichtlichkeit
    n_show = min(20, N)
    n_indices = np.arange(n_show)
    
    E_rest = energies['E_rest_segments'][:n_show].to(u.J).value
    E_GR = energies['E_GR_segments'][:n_show].to(u.J).value
    E_SR = energies['E_SR_segments'][:n_show].to(u.J).value
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    width = 0.25
    ax.bar(n_indices - width, E_rest, width, label='E_rest(n)', alpha=0.8)
    ax.bar(n_indices, E_GR, width, label='E_GR(n)', alpha=0.8)
    ax.bar(n_indices + width, E_SR, width, label='E_SR(n)', alpha=0.8)
    
    ax.set_xlabel('Segment n', fontsize=12)
    ax.set_ylabel('Energie [J]', fontsize=12)
    ax.set_title(f'Energien der ersten {n_show} Segmente: {result["name"]}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f'segment_energies_{result["category"]}.png', dpi=150, bbox_inches='tight')
    print(f"Plot gespeichert: segment_energies_{result['category']}.png")
    plt.show()


# =============================================================================
# Erklärung und Interpretation
# =============================================================================

def print_interpretation():
    """Drucke Interpretation der Ergebnisse."""
    
    print("\n" + "="*80)
    print("INTERPRETATION DER ERGEBNISSE")
    print("="*80)
    print("""
TYPISCHE GROESSENORDNUNGEN:

1. HAUPTREIHENSTERN (Sonne):
   - E_GR/E_rest ~ 10^-6 (sehr klein)
   - E_SR/E_rest ~ 10^-6 (sehr klein)
   -> Schwaches Gravitationsfeld, niedrige Geschwindigkeiten
   -> Relativistische Effekte vernachlaessigbar

2. WEISSER ZWERG:
   - E_GR/E_rest ~ 10^-3 bis 10^-2 (moderat)
   - E_SR/E_rest ~ 10^-4 (klein)
   -> Kompaktes Objekt: GR-Effekte werden wichtig
   -> R/r_s ~ 1000-10000 (noch weit vom Schwarzschild-Radius)

3. NEUTRONENSTERN:
   - E_GR/E_rest ~ 10^-2 bis 10^-1 (dominant!)
   - E_SR/E_rest ~ 10^-3 (signifikant)
   -> Extrem kompakt: R/r_s ~ 2-5
   -> GR-Effekte dominieren, ~10% der Ruheenergie!
   -> Nahe am Schwarzschild-Radius: gamma_GR >> 1

EINFLUSS DER SEGMENTIERUNG:

Die Anzahl der Segmente N beeinflusst die numerische Genauigkeit:
- Kleine N (z.B. 10): Grobe Approximation, schnell
- Mittlere N (z.B. 100): Gute Balance zwischen Genauigkeit und Geschwindigkeit
- Grosse N (z.B. 1000): Hohe Genauigkeit, langsamer

Fuer kompakte Objekte (Neutronensterne) benoetigt man mehr Segmente, da sich
gamma_GR und gamma_SR nahe am Schwarzschild-Radius stark aendern.

Die Segmentierung zeigt, wie sich relativistische Effekte radial verteilen:
- Innere Segmente (kleine r): Hohe gamma_GR, grosse GR-Beitraege
- Aeussere Segmente (grosse r): gamma_GR ~ 1, kleine GR-Beitraege
""")
    print("="*80)


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Hauptfunktion: Führe alle Berechnungen und Visualisierungen durch."""
    
    print("="*80)
    print("SEGMENTIERTES ENERGIE-MODELL MIT SR UND GR")
    print("="*80)
    print("\nStarte Berechnungen...\n")
    
    # Definiere Test-Objekte
    objects = test_objects()
    
    # Analysiere alle Objekte
    N = 100  # Anzahl Segmente
    results = []
    
    for key, obj_data in objects.items():
        print(f"Analysiere: {obj_data['name']}...")
        result = analyze_object(obj_data, N=N)
        results.append(result)
    
    # Drucke Ergebnisse
    print_results(results)
    
    # Visualisierungen
    print("\nErstelle Visualisierungen...")
    
    # Plot 1: Energie-Normalisierung für alle Objekte
    plot_energy_normalization(results)
    
    # Plot 2: Gamma-Faktoren für Neutronenstern (interessantester Fall)
    ns_result = [r for r in results if r['category'] == 'neutron_star'][0]
    plot_gamma_factors(ns_result)
    
    # Plot 3: Segment-Energien für Neutronenstern
    plot_segment_energies(ns_result)
    
    # Interpretation
    print_interpretation()
    
    print("\n" + "="*80)
    print("ALLE BERECHNUNGEN ABGESCHLOSSEN!")
    print("="*80)
    print("\nDateien erstellt:")
    print("  - energy_normalization.png")
    print("  - gamma_factors_neutron_star.png")
    print("  - segment_energies_neutron_star.png")
    print("="*80)


if __name__ == "__main__":
    main()
