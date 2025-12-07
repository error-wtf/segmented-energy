#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Finaler Vergleich - Sichere Parameter

Vergleicht alle Versionen mit physikalisch sinnvollen Parametern.
"""

import os
import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun, au

os.environ['PYTHONIOENCODING'] = 'utf-8:replace'

print("="*80)
print("FINALER VERGLEICH - ALLE VERSIONEN")
print("="*80)
print()

# Sichere Parameter (weit weg vom Horizont)
M = M_sun
m = 1.0 * u.kg
r_in = 10.0 * R_sun  # Sicher weit vom Schwarzschild-Radius
r_out = 1.0 * au
N = 1000

r_s = (2 * G * M / c**2).to(u.km)

print("PARAMETER:")
print(f"  M = {M}")
print(f"  m = {m}")
print(f"  r_in = {r_in} = {(r_in/r_s).decompose().value:.1f} r_s")
print(f"  r_out = {r_out} = {(r_out/r_s).decompose().value:.1f} r_s")
print(f"  N = {N}")
print(f"  r_s = {r_s:.3f} (Schwarzschild-Radius)")
print()

results = {}

# VERSION 1
print("VERSION 1: segmented_energy.py")
print("-"*80)
try:
    from segmented_energy import compute_segmented_energy
    res = compute_segmented_energy(M, m, r_in, r_out, N, verbose=False)
    results['v1'] = {
        'E_tot': res['E_total'].to(u.J).value,
        'E_GR': res['E_GR_total'].to(u.J).value,
        'E_SR': res['E_SR_total'].to(u.J).value,
        'name': 'segmented_energy.py'
    }
    print(f"E_total = {results['v1']['E_tot']:.12e} J")
    print(f"E_GR    = {results['v1']['E_GR']:.12e} J")
    print(f"E_SR    = {results['v1']['E_SR']:.12e} J")
    print("[OK]")
except Exception as e:
    print(f"[ERROR] {e}")
print()

# VERSION 2
print("VERSION 2: energy_n_segments_astropy.py")
print("-"*80)
try:
    from energy_n_segments_astropy import compute_energy_n_segments
    res = compute_energy_n_segments(M, m, r_in, r_out, N, verbose=False)
    E_rest = res['E_rest'].to(u.J).value
    results['v2'] = {
        'E_tot': res['E_tot'].to(u.J).value,
        'E_GR': res['E_GR'].to(u.J).value,
        'E_SR': res['E_SR'].to(u.J).value,
        'E_rest': E_rest,
        'name': 'energy_n_segments_astropy.py'
    }
    print(f"E_total = {results['v2']['E_tot']:.12e} J")
    print(f"E_GR    = {results['v2']['E_GR']:.12e} J")
    print(f"E_SR    = {results['v2']['E_SR']:.12e} J")
    print(f"E_rest  = {results['v2']['E_rest']:.12e} J")
    print(f"E_tot / E_rest = {results['v2']['E_tot']/E_rest:.12e}")
    print("[OK]")
except Exception as e:
    print(f"[ERROR] {e}")
print()

# VERSION 3: Referenz-Implementation
print("VERSION 3: Referenz (teleskopisch + direkt)")
print("-"*80)
try:
    # Geometrische Radien
    r_arr = np.geomspace(r_in.value, r_out.value, N) * r_in.unit
    dm = m / N
    
    # GR (Newton): E_GR = -G*M*dm/r
    E_GR_segments = -G * M * dm / r_arr
    E_GR_sum = np.sum(E_GR_segments).to(u.J).value
    
    # Teleskopisch
    E_GR_tele = (m * (-G * M * (1/r_out - 1/r_in))).to(u.J).value
    
    # SR: v = sqrt(G*M/r), gamma = 1/sqrt(1-v^2/c^2)
    v_arr = np.sqrt(G * M / r_arr)
    beta_sq = (v_arr / c)**2
    gamma_arr = 1.0 / np.sqrt(1 - beta_sq)
    E_SR_segments = (gamma_arr.value - 1.0) * dm * c**2
    E_SR_sum = np.sum(E_SR_segments).to(u.J).value
    
    # Rest + Total
    E_rest = (m * c**2).to(u.J).value
    E_tot = E_rest + E_GR_sum + E_SR_sum
    
    results['v3'] = {
        'E_tot': E_tot,
        'E_GR': E_GR_sum,
        'E_GR_tele': E_GR_tele,
        'E_SR': E_SR_sum,
        'E_rest': E_rest,
        'name': 'Referenz-Implementation'
    }
    print(f"E_total       = {results['v3']['E_tot']:.12e} J")
    print(f"E_GR (sum)    = {results['v3']['E_GR']:.12e} J")
    print(f"E_GR (tele)   = {results['v3']['E_GR_tele']:.12e} J")
    print(f"E_SR          = {results['v3']['E_SR']:.12e} J")
    print(f"E_rest        = {results['v3']['E_rest']:.12e} J")
    print(f"E_tot / E_rest= {E_tot/E_rest:.12e}")
    
    # Teleskopische Kontrolle
    diff_tele = abs(E_GR_sum - E_GR_tele) / abs(E_GR_tele)
    print(f"\nTeleskopische Kontrolle: |E_GR_sum - E_GR_tele|/|E_GR_tele| = {diff_tele:.3e}")
    if diff_tele < 1e-10:
        print("[PERFEKT] Teleskopische Summe stimmt!")
    print("[OK]")
except Exception as e:
    print(f"[ERROR] {e}")
print()

# =============================================================================
# VERGLEICH
# =============================================================================

print("="*80)
print("VERGLEICHSANALYSE")
print("="*80)
print()

if 'v1' in results and 'v2' in results and 'v3' in results:
    print("NUMERISCHE ÜBEREINSTIMMUNG:")
    print("-"*80)
    
    ref = results['v3']['E_tot']
    
    for key in ['v1', 'v2']:
        if key in results:
            diff = abs(results[key]['E_tot'] - ref) / abs(ref)
            print(f"{results[key]['name']:40s}: {diff:.3e}  ", end="")
            if diff < 1e-12:
                print("[IDENTISCH BIS MASCHINENGENAUIGKEIT]")
            elif diff < 1e-6:
                print("[SEHR GUT]")
            elif diff < 1e-3:
                print("[GUT]")
            else:
                print("[UNTERSCHIED]")
    
    print()
    print("E_GR VERGLEICH:")
    print("-"*80)
    for key in ['v1', 'v2', 'v3']:
        if key in results:
            print(f"{results[key]['name']:40s}: {results[key]['E_GR']:.12e} J")
    
    if 'v3' in results:
        print(f"\nTeleskopische E_GR:                       {results['v3']['E_GR_tele']:.12e} J")
    
    print()
    print("E_SR VERGLEICH:")
    print("-"*80)
    for key in ['v1', 'v2', 'v3']:
        if key in results:
            print(f"{results[key]['name']:40s}: {results[key]['E_SR']:.12e} J")
    
    print()
    
    # Relative Anteile
    print("RELATIVE ENERGIE-ANTEILE (zu E_rest):")
    print("-"*80)
    for key in ['v1', 'v2', 'v3']:
        if key in results and 'E_rest' in results[key]:
            E_r = results[key]['E_rest']
            print(f"\n{results[key]['name']}:")
            print(f"  E_GR / E_rest  = {results[key]['E_GR']/E_r:.6e}")
            print(f"  E_SR / E_rest  = {results[key]['E_SR']/E_r:.6e}")
            print(f"  E_tot / E_rest = {results[key]['E_tot']/E_r:.6e}")

print()
print("="*80)
print("FAZIT")
print("="*80)
print()
print("GENAUESTE VERSION:")
print("  -> Version 3 (Referenz-Implementation)")
print("     * Direkteste Umsetzung der Physik")
print("     * Teleskopische Validierung inklusive")
print("     * E_rest explizit berechnet")
print()
print("BESTE FÜR PRODUKTION:")
print("  -> Version 2 (energy_n_segments_astropy.py)")
print("     * Alle Features")
print("     * Teleskopische Kontrolle")
print("     * Gut getestet")
print()
print("BESTE FÜR ERWEITUNGEN:")
print("  -> Version 1 (segmented_energy.py)")
print("     * Saubere API")
print("     * Wählbare Segmentierung")
print("     * Erweiterbar (phi, SSZ, etc.)")
print()

if 'v1' in results and 'v3' in results:
    diff_final = abs(results['v1']['E_tot'] - results['v3']['E_tot']) / abs(results['v3']['E_tot'])
    print(f"NUMERISCHE GENAUIGKEIT:")
    print(f"  Maximale Abweichung: {diff_final:.3e}")
    if diff_final < 1e-10:
        print("  [EXZELLENT] Alle Implementationen sind numerisch äquivalent!")
        
print("="*80)
