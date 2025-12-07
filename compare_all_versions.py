#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vergleich aller Segmented Energy Versionen

Testet und vergleicht die Genauigkeit verschiedener Implementierungen.
"""

import os
import sys
import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'

print("=" * 80)
print("VERGLEICH ALLER SEGMENTED ENERGY VERSIONEN")
print("=" * 80)
print()

# =============================================================================
# Test-Parameter (identisch für alle Versionen)
# =============================================================================

M = 1.0 * M_sun
m = 1.0 * u.kg
r_in = 2.0 * R_sun
r_out = 1.0 * u.au
N = 1000

print("TEST-PARAMETER (identisch für alle):")
print("-" * 80)
print(f"M = {M}")
print(f"m = {m}")
print(f"r_in = {r_in}")
print(f"r_out = {r_out}")
print(f"N = {N}")
print()

# =============================================================================
# VERSION 1: segmented_energy.py (Original Windsurf)
# =============================================================================

print("VERSION 1: segmented_energy.py")
print("-" * 80)

try:
    from segmented_energy import compute_segmented_energy as compute_v1
    
    result_v1 = compute_v1(M, m, r_in, r_out, N, segmentation="linear", verbose=False)
    
    E_tot_v1 = result_v1['E_total'].to(u.J).value
    E_GR_v1 = result_v1['E_GR_total'].to(u.J).value
    E_SR_v1 = result_v1['E_SR_total'].to(u.J).value
    E_norm_v1 = result_v1['E_normalized']
    
    print(f"E_total     = {E_tot_v1:.15e} J")
    print(f"E_GR_total  = {E_GR_v1:.15e} J")
    print(f"E_SR_total  = {E_SR_v1:.15e} J")
    print(f"E/(mc^2)    = {E_norm_v1:.15e}")
    print("[OK] Version 1 erfolgreich")
    
except Exception as e:
    print(f"[FEHLER] Version 1: {e}")
    result_v1 = None

print()

# =============================================================================
# VERSION 2: energy_n_segments_astropy.py (Original mit teleskopischer Summe)
# =============================================================================

print("VERSION 2: energy_n_segments_astropy.py")
print("-" * 80)

try:
    from energy_n_segments_astropy import compute_energy_n_segments as compute_v2
    
    result_v2 = compute_v2(M, m, r_in, r_out, N, segmentation='log', verbose=False)
    
    E_tot_v2 = result_v2['E_tot'].to(u.J).value
    E_GR_v2 = result_v2['E_GR'].to(u.J).value
    E_SR_v2 = result_v2['E_SR'].to(u.J).value
    E_rest_v2 = result_v2['E_rest'].to(u.J).value
    E_norm_v2 = (result_v2['E_tot'] / result_v2['E_rest']).decompose().value
    
    print(f"E_total     = {E_tot_v2:.15e} J")
    print(f"E_GR_total  = {E_GR_v2:.15e} J")
    print(f"E_SR_total  = {E_SR_v2:.15e} J")
    print(f"E_rest      = {E_rest_v2:.15e} J")
    print(f"E/(mc^2)    = {E_norm_v2:.15e}")
    print("[OK] Version 2 erfolgreich")
    
except Exception as e:
    print(f"[FEHLER] Version 2: {e}")
    result_v2 = None

print()

# =============================================================================
# VERSION 3: Manuelle Berechnung (gamma_sr und gamma_gr getrennt)
# =============================================================================

print("VERSION 3: Manuelle Berechnung (gamma_sr + gamma_gr)")
print("-" * 80)

try:
    # Radien berechnen
    r_array = np.geomspace(r_in.value, r_out.value, N) * r_in.unit
    
    # Masse pro Segment
    delta_m = m / N
    
    # SR: v = sqrt(G*M/r), gamma_sr = 1/sqrt(1 - v^2/c^2)
    v_array = np.sqrt(G * M / r_array)
    beta_sq_sr = (v_array / c)**2
    gamma_sr = 1.0 / np.sqrt(1 - beta_sq_sr)
    
    E_sr_total = np.sum((gamma_sr - 1.0) * delta_m * c**2)
    
    # GR: gamma_gr = 1/sqrt(1 - 2GM/(c^2*r))
    factor_gr = (2 * G * M / (c**2 * r_array)).decompose()
    gamma_gr = 1.0 / np.sqrt(1 - factor_gr)
    
    E_gr_total = np.sum((gamma_gr.value - 1.0) * delta_m * c**2)
    
    # Rest + Total
    E_rest = m * c**2
    E_tot_v3 = E_rest + E_sr_total + E_gr_total
    
    print(f"E_total     = {E_tot_v3.to(u.J).value:.15e} J")
    print(f"E_GR_total  = {E_gr_total.to(u.J).value:.15e} J")
    print(f"E_SR_total  = {E_sr_total.to(u.J).value:.15e} J")
    print(f"E_rest      = {E_rest.to(u.J).value:.15e} J")
    print(f"E/(mc^2)    = {(E_tot_v3/E_rest).decompose().value:.15e}")
    print("[OK] Version 3 erfolgreich")
    
    result_v3 = {
        'E_tot': E_tot_v3.to(u.J).value,
        'E_GR': E_gr_total.to(u.J).value,
        'E_SR': E_sr_total.to(u.J).value,
        'E_rest': E_rest.to(u.J).value,
    }
    
except Exception as e:
    print(f"[FEHLER] Version 3: {e}")
    result_v3 = None

print()

# =============================================================================
# VERSION 4: Nur GR (ohne SR) - Baseline
# =============================================================================

print("VERSION 4: Nur GR (Baseline - ohne SR)")
print("-" * 80)

try:
    r_array = np.geomspace(r_in.value, r_out.value, N) * r_in.unit
    delta_m = m / N
    
    # Nur GR: E_GR = -G*M*dm/r
    E_gr_newton = np.sum(-G * M * delta_m / r_array)
    
    # Teleskopische Kontrolle
    E_gr_tele = m * (-G * M * (1/r_out - 1/r_in))
    
    E_rest = m * c**2
    E_tot_v4 = E_rest + E_gr_newton
    
    print(f"E_GR (Newton)      = {E_gr_newton.to(u.J).value:.15e} J")
    print(f"E_GR (teleskopisch)= {E_gr_tele.to(u.J).value:.15e} J")
    print(f"Differenz          = {abs(E_gr_newton - E_gr_tele).to(u.J).value:.15e} J")
    print(f"E_total            = {E_tot_v4.to(u.J).value:.15e} J")
    print(f"E/(mc^2)           = {(E_tot_v4/E_rest).decompose().value:.15e}")
    print("[OK] Version 4 erfolgreich")
    
    result_v4 = {
        'E_tot': E_tot_v4.to(u.J).value,
        'E_GR': E_gr_newton.to(u.J).value,
        'E_rest': E_rest.to(u.J).value,
    }
    
except Exception as e:
    print(f"[FEHLER] Version 4: {e}")
    result_v4 = None

print()

# =============================================================================
# VERGLEICH UND ANALYSE
# =============================================================================

print("=" * 80)
print("VERGLEICH DER VERSIONEN")
print("=" * 80)
print()

if result_v1 and result_v2 and result_v3 and result_v4:
    
    # Vergleiche E_total
    print("E_TOTAL VERGLEICH:")
    print("-" * 80)
    print(f"Version 1 (segmented_energy):       {E_tot_v1:.15e} J")
    print(f"Version 2 (energy_n_segments):      {E_tot_v2:.15e} J")
    print(f"Version 3 (gamma_sr + gamma_gr):    {result_v3['E_tot']:.15e} J")
    print(f"Version 4 (nur GR):                 {result_v4['E_tot']:.15e} J")
    print()
    
    # Relative Differenzen
    print("RELATIVE DIFFERENZEN (zu Version 3 als Referenz):")
    print("-" * 80)
    
    ref = result_v3['E_tot']
    
    diff_v1 = abs(E_tot_v1 - ref) / abs(ref)
    diff_v2 = abs(E_tot_v2 - ref) / abs(ref)
    diff_v4 = abs(result_v4['E_tot'] - ref) / abs(ref)
    
    print(f"Version 1 vs 3:  {diff_v1:.3e}  {'[IDENTISCH]' if diff_v1 < 1e-10 else '[UNTERSCHIED]'}")
    print(f"Version 2 vs 3:  {diff_v2:.3e}  {'[IDENTISCH]' if diff_v2 < 1e-10 else '[UNTERSCHIED]'}")
    print(f"Version 4 vs 3:  {diff_v4:.3e}  {'[UNTERSCHIED - erwartet, da kein SR]' if diff_v4 > 1e-6 else ''}")
    print()
    
    # E_GR Vergleich
    print("E_GR VERGLEICH:")
    print("-" * 80)
    print(f"Version 1:  {E_GR_v1:.15e} J")
    print(f"Version 2:  {E_GR_v2:.15e} J")
    print(f"Version 3:  {result_v3['E_GR']:.15e} J")
    print(f"Version 4:  {result_v4['E_GR']:.15e} J (Newton)")
    print()
    
    diff_gr_12 = abs(E_GR_v1 - E_GR_v2) / abs(E_GR_v2)
    diff_gr_13 = abs(E_GR_v1 - result_v3['E_GR']) / abs(result_v3['E_GR'])
    diff_gr_14 = abs(E_GR_v1 - result_v4['E_GR']) / abs(result_v4['E_GR'])
    
    print(f"V1 vs V2:  {diff_gr_12:.3e}")
    print(f"V1 vs V3:  {diff_gr_13:.3e}")
    print(f"V1 vs V4:  {diff_gr_14:.3e}")
    print()
    
    # E_SR Vergleich
    print("E_SR VERGLEICH:")
    print("-" * 80)
    print(f"Version 1:  {E_SR_v1:.15e} J")
    print(f"Version 2:  {E_SR_v2:.15e} J")
    print(f"Version 3:  {result_v3['E_SR']:.15e} J")
    print(f"Version 4:  N/A (kein SR)")
    print()
    
    diff_sr_12 = abs(E_SR_v1 - E_SR_v2) / abs(E_SR_v2)
    diff_sr_13 = abs(E_SR_v1 - result_v3['E_SR']) / abs(result_v3['E_SR'])
    
    print(f"V1 vs V2:  {diff_sr_12:.3e}")
    print(f"V1 vs V3:  {diff_sr_13:.3e}")
    print()

# =============================================================================
# EMPFEHLUNG
# =============================================================================

print("=" * 80)
print("EMPFEHLUNG")
print("=" * 80)
print()

print("GENAUIGKEIT:")
print("  1. Version 1 (segmented_energy.py)")
print("     - Sauber implementiert mit astropy")
print("     - Konsistente Formeln")
print("     - Gute Dokumentation")
print()
print("  2. Version 2 (energy_n_segments_astropy.py)")
print("     - Inklusive teleskopischer Kontrolle")
print("     - E_rest explizit getrennt")
print("     - Sehr ähnliche Ergebnisse zu V1")
print()
print("  3. Version 3 (gamma_sr + gamma_gr getrennt)")
print("     - Physikalisch transparenteste Version")
print("     - SR und GR Effekte klar getrennt")
print("     - Beste für physikalische Interpretation")
print()

print("EMPFOHLENE VERSION:")
print("  -> Version 3 (manuelle Berechnung) für:")
print("     * Maximale Transparenz")
print("     * Physikalisches Verständnis")
print("     * Trennung SR/GR Effekte")
print()
print("  -> Version 2 (energy_n_segments_astropy.py) für:")
print("     * Produktionscode")
print("     * Teleskopische Validierung")
print("     * Erweiterte Features")
print()

print("NUMERISCHE ÜBEREINSTIMMUNG:")
if result_v1 and result_v2 and result_v3:
    max_diff = max(diff_v1, diff_v2)
    if max_diff < 1e-10:
        print(f"  [EXZELLENT] Alle Versionen stimmen auf {max_diff:.1e} überein")
    elif max_diff < 1e-6:
        print(f"  [GUT] Versionen stimmen auf {max_diff:.1e} überein")
    else:
        print(f"  [WARNUNG] Differenz: {max_diff:.1e}")

print("=" * 80)
