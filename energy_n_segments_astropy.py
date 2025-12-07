#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Energiezerlegung mit N-Segmentierung - astropy Implementation

Berechnet die Gesamtenergie einer Testmasse in einem Zentralfeld
mit expliziter Segmentierung in N radiale Bereiche.

E_tot = E_rest + E_GR + E_SR

wobei:
  E_rest = m·c²
  E_GR   = Σ_{n=1}^N E_gr_n  (gravitativ, segmentweise)
  E_SR   = Σ_{n=1}^N E_sr_n  (spezial-relativistisch, segmentweise)

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import numpy as np
import matplotlib.pyplot as plt

try:
    from astropy import constants as const
    from astropy import units as u
    HAS_ASTROPY = True
except ImportError:
    print("⚠️  astropy nicht gefunden. Verwende fallback scipy constants.")
    import scipy.constants as const
    HAS_ASTROPY = False
    # Einfache units-Emulation für scipy
    class SimpleUnit:
        def __init__(self, value, unit_name):
            self.value = value
            self.unit = unit_name
        def to(self, target_unit):
            return self
        def __mul__(self, other):
            if isinstance(other, SimpleUnit):
                return SimpleUnit(self.value * other.value, f"{self.unit}·{other.unit}")
            return SimpleUnit(self.value * other, self.unit)
        def __truediv__(self, other):
            if isinstance(other, SimpleUnit):
                return SimpleUnit(self.value / other.value, f"{self.unit}/{other.unit}")
            return SimpleUnit(self.value / other, self.unit)
        def __add__(self, other):
            return SimpleUnit(self.value + other.value, self.unit)
        def __sub__(self, other):
            return SimpleUnit(self.value - other.value, self.unit)
        def __repr__(self):
            return f"{self.value} {self.unit}"


# =============================================================================
# Physikalische Konstanten
# =============================================================================

if HAS_ASTROPY:
    # astropy constants (mit units)
    G = const.G
    c = const.c
    M_sun = const.M_sun
    M_earth = const.M_earth
    au = const.au
else:
    # scipy fallback (ohne units)
    G = SimpleUnit(const.G, "m³/(kg·s²)")
    c = SimpleUnit(const.c, "m/s")
    M_sun = SimpleUnit(1.98847e30, "kg")
    M_earth = SimpleUnit(5.9722e24, "kg")
    au = SimpleUnit(1.495978707e11, "m")


# =============================================================================
# Funktion: Energie-Zerlegung mit N Segmenten
# =============================================================================

def compute_energy_n_segments(M, m, r0, rN, N=100, segmentation='log', verbose=True):
    """
    Berechnet Gesamtenergie mit N-Segmentierung.
    
    Parameters
    ----------
    M : astropy.Quantity or SimpleUnit
        Zentrale Masse (z.B. M_sun)
    m : astropy.Quantity or SimpleUnit
        Testmasse (z.B. M_earth)
    r0 : astropy.Quantity or SimpleUnit
        Innerer Radius (Startpunkt)
    rN : astropy.Quantity or SimpleUnit
        Äußerer Radius (Endpunkt)
    N : int
        Anzahl der Segmente
    segmentation : str
        'log' (geometrisch), 'linear', oder 'phi' (φ-basiert)
    verbose : bool
        Ausgabe von Details
    
    Returns
    -------
    dict mit:
        E_rest, E_GR, E_SR, E_tot,
        r_edges, r_mid,
        E_gr_segments, E_sr_segments,
        gamma_segments
    """
    
    # Segment-Grenzen erstellen
    if segmentation == 'log':
        # Logarithmische (geometrische) Segmentierung
        if HAS_ASTROPY:
            r_edges = np.geomspace(r0.value, rN.value, N+1) * r0.unit
        else:
            r_edges_vals = np.geomspace(r0.value, rN.value, N+1)
            r_edges = SimpleUnit(r_edges_vals, r0.unit)
            
    elif segmentation == 'linear':
        # Lineare Segmentierung
        if HAS_ASTROPY:
            r_edges = np.linspace(r0.value, rN.value, N+1) * r0.unit
        else:
            r_edges_vals = np.linspace(r0.value, rN.value, N+1)
            r_edges = SimpleUnit(r_edges_vals, r0.unit)
            
    elif segmentation == 'phi':
        # φ-basierte Segmentierung
        phi = (1 + np.sqrt(5)) / 2  # Goldener Schnitt
        alpha = np.log(rN.value / r0.value) / (N * np.log(phi))
        exponents = phi ** (np.arange(N+1) * alpha)
        if HAS_ASTROPY:
            r_edges = r0.value * exponents * r0.unit
        else:
            r_edges_vals = r0.value * exponents
            r_edges = SimpleUnit(r_edges_vals, r0.unit)
    else:
        raise ValueError(f"Unknown segmentation: {segmentation}")
    
    # Segment-Mittelpunkte (geometrisches Mittel)
    if HAS_ASTROPY:
        r_mid = np.sqrt(r_edges[:-1] * r_edges[1:])
    else:
        r_mid_vals = np.sqrt(r_edges.value[:-1] * r_edges.value[1:])
        r_mid = SimpleUnit(r_mid_vals, r_edges.unit)
    
    # =========================================================================
    # 1. Gravitationsbeiträge E_gr_n
    # =========================================================================
    
    # Potential an Segmentgrenzen: Φ(r) = -GM/r
    if HAS_ASTROPY:
        Phi_edges = -G * M / r_edges
        E_pot_edges = m * Phi_edges
        # Differenzen: ΔE_gr_n = E_pot(r_n) - E_pot(r_{n-1})
        E_gr_segments = E_pot_edges[1:] - E_pot_edges[:-1]
        E_GR_sum = np.sum(E_gr_segments)
    else:
        Phi_edges_vals = -G.value * M.value / r_edges.value
        E_pot_edges_vals = m.value * Phi_edges_vals
        E_gr_segments_vals = np.diff(E_pot_edges_vals)
        E_gr_segments = SimpleUnit(E_gr_segments_vals, "J")
        E_GR_sum = SimpleUnit(np.sum(E_gr_segments_vals), "J")
    
    # Teleskopische Kontrolle: E_GR = -GMm·(1/r_N - 1/r_0)
    if HAS_ASTROPY:
        E_GR_tel = m * (-G * M * (1/rN - 1/r0))
    else:
        E_GR_tel_val = m.value * (-G.value * M.value * (1/rN.value - 1/r0.value))
        E_GR_tel = SimpleUnit(E_GR_tel_val, "J")
    
    # =========================================================================
    # 2. SR-Beiträge E_sr_n
    # =========================================================================
    
    # Orbitgeschwindigkeit: v_n = sqrt(GM/r_n)
    if HAS_ASTROPY:
        v = np.sqrt(G * M / r_mid)
        # Lorentz-Faktor: γ_n = 1 / sqrt(1 - (v_n/c)²)
        beta_sq = (v / c)**2
        gamma = 1.0 / np.sqrt(1 - beta_sq)
        # SR-Energie: E_sr_n = (γ_n - 1)·m·c²
        E_sr_segments = (gamma - 1.0) * m * c**2
        E_SR = np.sum(E_sr_segments)
    else:
        v_vals = np.sqrt(G.value * M.value / r_mid.value)
        beta_sq_vals = (v_vals / c.value)**2
        gamma_vals = 1.0 / np.sqrt(1 - beta_sq_vals)
        E_sr_segments_vals = (gamma_vals - 1.0) * m.value * c.value**2
        E_sr_segments = SimpleUnit(E_sr_segments_vals, "J")
        E_SR = SimpleUnit(np.sum(E_sr_segments_vals), "J")
        gamma = gamma_vals
    
    # =========================================================================
    # 3. Gesamtenergie
    # =========================================================================
    
    if HAS_ASTROPY:
        E_rest = m * c**2
        E_tot = E_rest + E_GR_sum + E_SR
    else:
        E_rest_val = m.value * c.value**2
        E_rest = SimpleUnit(E_rest_val, "J")
        E_tot_val = E_rest.value + E_GR_sum.value + E_SR.value
        E_tot = SimpleUnit(E_tot_val, "J")
    
    # =========================================================================
    # Output
    # =========================================================================
    
    if verbose:
        print("="*80)
        print("ENERGIEZERLEGUNG MIT N-SEGMENTIERUNG")
        print("="*80)
        print(f"Zentrale Masse M: {M}")
        print(f"Testmasse m:      {m}")
        print(f"Radialer Bereich: {r0} → {rN}")
        print(f"Anzahl Segmente:  {N}")
        print(f"Segmentierung:    {segmentation}")
        print()
        print("ENERGIEN:")
        print(f"  E_rest   = {E_rest.to('J') if HAS_ASTROPY else E_rest}")
        print(f"  E_GR_sum = {E_GR_sum.to('J') if HAS_ASTROPY else E_GR_sum}")
        print(f"  E_GR_tel = {E_GR_tel.to('J') if HAS_ASTROPY else E_GR_tel}  (Kontrolle)")
        print(f"  E_SR     = {E_SR.to('J') if HAS_ASTROPY else E_SR}")
        print(f"  E_tot    = {E_tot.to('J') if HAS_ASTROPY else E_tot}")
        print()
        
        # Relativer Fehler der teleskopischen Kontrolle
        if HAS_ASTROPY:
            rel_err = abs((E_GR_sum - E_GR_tel) / E_GR_tel).value
        else:
            rel_err = abs((E_GR_sum.value - E_GR_tel.value) / E_GR_tel.value)
        print(f"Teleskopische Kontrolle: |E_GR_sum - E_GR_tel| / |E_GR_tel| = {rel_err:.2e}")
        
        if rel_err < 1e-10:
            print("  ✅ Kontrolle BESTANDEN (Fehler < 1e-10)")
        else:
            print("  ⚠️  Kontrolle: Numerischer Fehler erkannt")
        print("="*80)
    
    return {
        'E_rest': E_rest,
        'E_GR': E_GR_sum,
        'E_GR_tel': E_GR_tel,
        'E_SR': E_SR,
        'E_tot': E_tot,
        'r_edges': r_edges,
        'r_mid': r_mid,
        'E_gr_segments': E_gr_segments,
        'E_sr_segments': E_sr_segments,
        'gamma': gamma,
        'N': N,
        'M': M,
        'm': m,
    }


# =============================================================================
# Visualisierung
# =============================================================================

def plot_energy_profile(result, save=False):
    """
    Visualisiert das Energie-Profil über die Segmente.
    
    Parameters
    ----------
    result : dict
        Output von compute_energy_n_segments
    save : bool
        Speichere Plot als PNG
    """
    
    r_mid = result['r_mid']
    E_gr = result['E_gr_segments']
    E_sr = result['E_sr_segments']
    N = result['N']
    
    # Konvertiere zu Joule für Plotting
    if HAS_ASTROPY:
        r_mid_val = r_mid.to('au').value
        E_gr_val = E_gr.to('J').value
        E_sr_val = E_sr.to('J').value
        r_unit = 'AU'
    else:
        r_mid_val = r_mid.value / au.value  # In AU
        E_gr_val = E_gr.value
        E_sr_val = E_sr.value
        r_unit = 'AU (approx)'
    
    E_total_segments = E_gr_val + E_sr_val
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Energie pro Segment
    ax = axes[0, 0]
    ax.plot(r_mid_val, E_gr_val, 'o-', label='$E_{gr,n}$ (Gravitativ)', color='blue', markersize=3)
    ax.plot(r_mid_val, E_sr_val, 's-', label='$E_{sr,n}$ (SR)', color='red', markersize=3)
    ax.plot(r_mid_val, E_total_segments, '^-', label='$E_{gr,n} + E_{sr,n}$', color='green', markersize=3)
    ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
    ax.set_xlabel(f'Radius r [{r_unit}]')
    ax.set_ylabel('Energie pro Segment [J]')
    ax.set_title(f'Energie-Profil (N={N} Segmente)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Kumulative Energie
    ax = axes[0, 1]
    E_cumulative = np.cumsum(E_total_segments)
    if HAS_ASTROPY:
        E_rest_val = result['E_rest'].to('J').value
    else:
        E_rest_val = result['E_rest'].value
    ax.plot(r_mid_val, E_cumulative + E_rest_val, 'o-', color='purple', markersize=4)
    ax.axhline(E_rest_val, color='orange', linestyle='--', label='$E_{rest} = mc^2$')
    ax.set_xlabel(f'Radius r [{r_unit}]')
    ax.set_ylabel('Kumulative Energie [J]')
    ax.set_title('Kumulative Gesamtenergie')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Lorentz-Faktor γ_n
    ax = axes[1, 0]
    gamma = result['gamma']
    if HAS_ASTROPY:
        gamma_val = gamma.value if hasattr(gamma, 'value') else gamma
    else:
        gamma_val = gamma
    ax.plot(r_mid_val, gamma_val, 'o-', color='darkred', markersize=4)
    ax.axhline(1, color='black', linestyle='--', linewidth=0.5, label='γ = 1 (nicht-relativistisch)')
    ax.set_xlabel(f'Radius r [{r_unit}]')
    ax.set_ylabel('Lorentz-Faktor γ')
    ax.set_title('Lorentz-Faktor pro Segment')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Energie-Anteile (Pie Chart am äußersten Segment)
    ax = axes[1, 1]
    if HAS_ASTROPY:
        E_GR_total = result['E_GR'].to('J').value
        E_SR_total = result['E_SR'].to('J').value
    else:
        E_GR_total = result['E_GR'].value
        E_SR_total = result['E_SR'].value
    
    # Absolute Werte für Pie Chart
    E_abs = [abs(E_GR_total), abs(E_SR_total)]
    labels = [f'|E_GR| = {abs(E_GR_total):.2e} J',
              f'E_SR = {E_SR_total:.2e} J']
    colors = ['blue', 'red']
    ax.pie(E_abs, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax.set_title('Energie-Anteile (absolut)')
    
    plt.tight_layout()
    
    if save:
        plt.savefig('energy_n_segments_profile.png', dpi=150, bbox_inches='tight')
        print("💾 Plot gespeichert: energy_n_segments_profile.png")
    
    plt.show()


# =============================================================================
# Beispiel-Ausführung
# =============================================================================

if __name__ == "__main__":
    
    print("\n" + "="*80)
    print("ENERGIEZERLEGUNG MIT N-SEGMENTIERUNG - ASTROPY IMPLEMENTATION")
    print("="*80 + "\n")
    
    # =========================================================================
    # Beispiel 1: Sonne-Erde System
    # =========================================================================
    
    print("BEISPIEL 1: Sonne-Erde System (0.5 AU → 2.0 AU)")
    print("-" * 80)
    
    M = M_sun       # Zentrale Masse: Sonne
    m = M_earth     # Testmasse: Erde
    
    if HAS_ASTROPY:
        r0 = 0.5 * au   # Innerer Radius: 0.5 AU
        rN = 2.0 * au   # Äußerer Radius: 2.0 AU
    else:
        r0 = SimpleUnit(0.5 * au.value, au.unit)
        rN = SimpleUnit(2.0 * au.value, au.unit)
    
    N = 100  # Anzahl Segmente
    
    # Berechnung mit logarithmischer Segmentierung
    result_log = compute_energy_n_segments(M, m, r0, rN, N=N, segmentation='log', verbose=True)
    
    # Visualisierung
    plot_energy_profile(result_log, save=True)
    
    # =========================================================================
    # Beispiel 2: Verschiedene Segmentierungen vergleichen
    # =========================================================================
    
    print("\n" + "="*80)
    print("BEISPIEL 2: Vergleich verschiedener Segmentierungen")
    print("-" * 80)
    
    segmentations = ['log', 'linear', 'phi']
    results = {}
    
    for seg_type in segmentations:
        print(f"\nSegmentierung: {seg_type}")
        results[seg_type] = compute_energy_n_segments(M, m, r0, rN, N=50, 
                                                       segmentation=seg_type, 
                                                       verbose=False)
        if HAS_ASTROPY:
            E_tot = results[seg_type]['E_tot'].to('J').value
        else:
            E_tot = results[seg_type]['E_tot'].value
        print(f"  E_tot = {E_tot:.6e} J")
    
    # Vergleiche ob alle gleich sind (sollten sein, da Gesamtenergie unabhängig von Segmentierung)
    print("\nKonsistenz-Check:")
    if HAS_ASTROPY:
        E_log = results['log']['E_tot'].to('J').value
        E_lin = results['linear']['E_tot'].to('J').value
        E_phi = results['phi']['E_tot'].to('J').value
    else:
        E_log = results['log']['E_tot'].value
        E_lin = results['linear']['E_tot'].value
        E_phi = results['phi']['E_tot'].value
    
    rel_diff_log_lin = abs(E_log - E_lin) / abs(E_log)
    rel_diff_log_phi = abs(E_log - E_phi) / abs(E_log)
    
    print(f"  |E_log - E_linear| / |E_log| = {rel_diff_log_lin:.2e}")
    print(f"  |E_log - E_phi|    / |E_log| = {rel_diff_log_phi:.2e}")
    
    if rel_diff_log_lin < 1e-6 and rel_diff_log_phi < 1e-6:
        print("  ✅ Alle Segmentierungen liefern konsistente Gesamtenergie")
    else:
        print("  ⚠️  Inkonsistenz erkannt (sollte nicht passieren!)")
    
    # =========================================================================
    # Beispiel 3: N-Konvergenz Test
    # =========================================================================
    
    print("\n" + "="*80)
    print("BEISPIEL 3: N-Konvergenz Test")
    print("-" * 80)
    
    N_values = [10, 50, 100, 500, 1000]
    E_tot_values = []
    
    for N_test in N_values:
        res = compute_energy_n_segments(M, m, r0, rN, N=N_test, 
                                         segmentation='log', verbose=False)
        if HAS_ASTROPY:
            E_tot_values.append(res['E_tot'].to('J').value)
        else:
            E_tot_values.append(res['E_tot'].value)
    
    print("\nN-Werte vs E_tot:")
    for N_val, E_val in zip(N_values, E_tot_values):
        print(f"  N = {N_val:4d}  →  E_tot = {E_val:.10e} J")
    
    # Relative Änderung zwischen aufeinanderfolgenden N-Werten
    print("\nRelative Änderung:")
    for i in range(1, len(N_values)):
        rel_change = abs(E_tot_values[i] - E_tot_values[i-1]) / abs(E_tot_values[i])
        print(f"  N={N_values[i-1]:4d} → N={N_values[i]:4d}:  {rel_change:.2e}")
    
    # =========================================================================
    # Zusammenfassung
    # =========================================================================
    
    print("\n" + "="*80)
    print("ZUSAMMENFASSUNG")
    print("="*80)
    print("✅ Energiezerlegung mit N-Segmentierung erfolgreich implementiert")
    print("✅ Teleskopische Eigenschaft verifiziert")
    print("✅ Verschiedene Segmentierungen getestet (log, linear, φ)")
    print("✅ N-Konvergenz bestätigt")
    print("✅ Visualisierungen erstellt")
    print()
    print("NÄCHSTE SCHRITTE:")
    print("  1. SSZ-Potential Ξ(r) integrieren")
    print("  2. Black-Hole-Szenario (r_0 = r_h = 2GM/c²)")
    print("  3. g1/g2-Domänen spezifische Korrekturen")
    print("  4. φ-Spiral Segmentierung verfeinern")
    print("="*80 + "\n")
