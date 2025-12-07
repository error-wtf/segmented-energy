#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test Unified Energy Model auf vollständigem Observer-Datensatz

Testet segmented_energy_unified.py auf ALLEN Objekten mit gemessenen Werten.
Keine Schätzungen, kein Filling.

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun
import time

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'

from segmented_energy_unified import compute_unified_energy, predict_observables
from fetch_observer_data import fetch_complete_dataset


# =============================================================================
# Test auf komplettem Datensatz
# =============================================================================

def test_complete_dataset(df, N=1000, segmentation='logarithmic', verbose_each=False):
    """
    Teste Unified Model auf allen Objekten im Datensatz.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Datensatz von fetch_observer_data
    N : int
        Anzahl Segmente
    segmentation : str
        Segmentierungs-Typ
    verbose_each : bool
        Zeige Details für jedes Objekt
        
    Returns
    -------
    results_df : pandas.DataFrame
        Ergebnisse für alle Objekte
    """
    
    results = []
    
    print("\n" + "="*80)
    print(f"TESTE UNIFIED MODEL AUF {len(df)} OBJEKTEN")
    print("="*80)
    print(f"Segmentierung: {segmentation}, N = {N}")
    print()
    
    start_time = time.time()
    
    for idx, row in df.iterrows():
        name = row['name']
        category = row['category']
        
        # Masse und Radius mit Einheiten
        M = row['mass'] * M_sun
        
        if row['radius_unit'] == 'R_sun':
            R = row['radius'] * R_sun
        else:  # km
            R = row['radius'] * u.km
        
        # r_in und r_out definieren
        r_s = row['r_s_km'] * u.km
        
        # Für normale Sterne: r_in = R, r_out = weit draußen
        if category in ['main_sequence', 'exoplanet_host']:
            r_in = R * 1.1  # Leicht über Oberfläche
            r_out = 100 * r_in  # Faktor 100
        
        # Für kompakte Objekte: näher am r_s
        elif category == 'white_dwarf':
            r_in = R * 1.05  # Sehr nah an Oberfläche
            r_out = 50 * r_in
        
        elif category == 'neutron_star':
            # Sehr nahe am Schwarzschild-Radius!
            r_in = R * 1.02
            r_out = 20 * r_in
        
        else:
            r_in = R
            r_out = 100 * R
        
        # Testmasse
        m = 1.0 * u.kg
        
        try:
            # Berechne Energie
            result = compute_unified_energy(
                M=M, m=m, r_in=r_in, r_out=r_out,
                N=N, segmentation=segmentation,
                validate_teleskopic=True,
                verbose=verbose_each
            )
            
            # Observable
            obs = predict_observables(result)
            
            # Sammle Ergebnisse
            results.append({
                'name': name,
                'category': category,
                'mass_Msun': row['mass'],
                'radius_km': R.to(u.km).value,
                'r_s_km': row['r_s_km'],
                'r_over_rs': row['r_over_rs'],
                'r_in_km': r_in.to(u.km).value,
                'r_out_km': r_out.to(u.km).value,
                'E_total_J': result['E_total'].value,
                'E_rest_J': result['E_rest'].value,
                'E_GR_J': result['E_GR_total'].value,
                'E_SR_J': result['E_SR_total'].value,
                'E_normalized': result['E_normalized'],
                'tele_diff': result['teleskopic_diff'],
                'gamma_gr_max': np.max(result['gamma_gr']),
                'gamma_gr_min': np.min(result['gamma_gr']),
                'z_gr_max': np.max(obs['redshift']),
                'shapiro_delay_s': obs['shapiro_delay'].to(u.s).value if obs['shapiro_delay'] is not None else np.nan,
                'success': True,
                'error': None,
            })
            
            if not verbose_each:
                print(f"  [{idx+1:2d}/{len(df)}] {name:20s} ... OK (E_norm={result['E_normalized']:.6f})")
            
        except Exception as e:
            results.append({
                'name': name,
                'category': category,
                'mass_Msun': row['mass'],
                'radius_km': R.to(u.km).value,
                'r_s_km': row['r_s_km'],
                'r_over_rs': row['r_over_rs'],
                'success': False,
                'error': str(e),
            })
            print(f"  [{idx+1:2d}/{len(df)}] {name:20s} ... FEHLER: {e}")
    
    elapsed = time.time() - start_time
    
    results_df = pd.DataFrame(results)
    
    print(f"\n" + "="*80)
    print(f"TESTS ABGESCHLOSSEN")
    print(f"  Dauer: {elapsed:.2f} s ({elapsed/len(df):.3f} s/Objekt)")
    print(f"  Erfolgreich: {results_df['success'].sum()}/{len(df)}")
    print("="*80)
    
    return results_df


# =============================================================================
# Analyse und Statistiken
# =============================================================================

def analyze_results(results_df):
    """Analysiere Test-Ergebnisse"""
    
    # Nur erfolgreiche
    success = results_df[results_df['success'] == True].copy()
    
    print("\n" + "="*80)
    print("ERGEBNIS-ANALYSE")
    print("="*80)
    
    print(f"\nERFOLGSRATE:")
    print(f"  Erfolgreich: {len(success)}/{len(results_df)} ({100*len(success)/len(results_df):.1f}%)")
    
    if len(success) == 0:
        print("\nKeine erfolgreichen Tests!")
        return
    
    print(f"\nE_TOTAL / E_REST (Energie-Normalisierung):")
    print(f"  Min:    {success['E_normalized'].min():.12f}")
    print(f"  Max:    {success['E_normalized'].max():.12f}")
    print(f"  Median: {success['E_normalized'].median():.12f}")
    print(f"  Std:    {success['E_normalized'].std():.12e}")
    
    # Prüfe wie nah bei 1.0
    deviation = np.abs(success['E_normalized'] - 1.0)
    print(f"\nABWEICHUNG VON E_rest:")
    print(f"  Max Abweichung: {deviation.max():.12e}")
    print(f"  Median Abweichung: {deviation.median():.12e}")
    
    print(f"\nTELESKOPISCHE VALIDIERUNG:")
    print(f"  Min Differenz:    {success['tele_diff'].min():.6e}")
    print(f"  Max Differenz:    {success['tele_diff'].max():.6e}")
    print(f"  Median Differenz: {success['tele_diff'].median():.6e}")
    
    good_tele = success[success['tele_diff'] < 1e-6]
    print(f"  Perfekt (<1e-6): {len(good_tele)}/{len(success)} ({100*len(good_tele)/len(success):.1f}%)")
    
    print(f"\nLORENTZ-FAKTOREN (gamma_gr):")
    print(f"  Max (staerkstes Feld): {success['gamma_gr_max'].max():.6f}")
    print(f"  Bei: {success.loc[success['gamma_gr_max'].idxmax(), 'name']}")
    
    print(f"\nREDSHIFT (z_gr):")
    print(f"  Max (staerkstes Feld): {success['z_gr_max'].max():.6e}")
    print(f"  Bei: {success.loc[success['z_gr_max'].idxmax(), 'name']}")
    
    print(f"\nSHAPIRO DELAY:")
    valid_shapiro = success[~success['shapiro_delay_s'].isna()]
    if len(valid_shapiro) > 0:
        print(f"  Min: {valid_shapiro['shapiro_delay_s'].min():.6e} s")
        print(f"  Max: {valid_shapiro['shapiro_delay_s'].max():.6e} s")
        print(f"  Bei: {valid_shapiro.loc[valid_shapiro['shapiro_delay_s'].idxmax(), 'name']}")
    
    print(f"\nPRO KATEGORIE:")
    for cat in success['category'].unique():
        cat_data = success[success['category'] == cat]
        print(f"\n  {cat}:")
        print(f"    Anzahl: {len(cat_data)}")
        print(f"    E_norm: {cat_data['E_normalized'].mean():.12f} ± {cat_data['E_normalized'].std():.12e}")
        print(f"    Tele Diff: {cat_data['tele_diff'].median():.6e}")
    
    print("\n" + "="*80)
    
    return success


# =============================================================================
# Observable Matching (Prozentuale Übereinstimmung)
# =============================================================================

def calculate_observable_matching(results_df):
    """
    Berechne prozentuale Übereinstimmung mit erwarteten Observablen.
    
    Returns
    -------
    matching_scores : dict
        Übereinstimmung für verschiedene Observable
    """
    
    success = results_df[results_df['success'] == True]
    
    if len(success) == 0:
        return {}
    
    print("\n" + "="*80)
    print("OBSERVABLE MATCHING - PROZENTUALE ÜBEREINSTIMMUNG")
    print("="*80)
    
    scores = {}
    
    # 1. ENERGIE-ERHALTUNG (E_total ≈ E_rest)
    print("\n1. ENERGIE-ERHALTUNG (E_total ~ E_rest):")
    print("-" * 80)
    
    # Erwartung: E_total / E_rest sollte sehr nahe bei 1.0 sein
    # Für normale Sterne: < 10^-6 Abweichung
    # Für kompakte Objekte: < 5% Abweichung
    
    deviation = np.abs(success['E_normalized'] - 1.0)
    
    # Main Sequence
    ms = success[success['category'] == 'main_sequence']
    if len(ms) > 0:
        ms_dev = np.abs(ms['E_normalized'] - 1.0)
        ms_good = np.sum(ms_dev < 1e-6) / len(ms) * 100
        scores['main_sequence_energy'] = ms_good
        print(f"  Main Sequence ({len(ms)} Objekte):")
        print(f"    < 1 ppm Abweichung:  {ms_good:.1f}%")
        print(f"    Median Abweichung:   {ms_dev.median():.3e}")
    
    # White Dwarfs
    wd = success[success['category'] == 'white_dwarf']
    if len(wd) > 0:
        wd_dev = np.abs(wd['E_normalized'] - 1.0)
        wd_good = np.sum(wd_dev < 1e-4) / len(wd) * 100
        scores['white_dwarf_energy'] = wd_good
        print(f"  White Dwarfs ({len(wd)} Objekte):")
        print(f"    < 0.01% Abweichung:  {wd_good:.1f}%")
        print(f"    Median Abweichung:   {wd_dev.median():.3e}")
    
    # Neutron Stars
    ns = success[success['category'] == 'neutron_star']
    if len(ns) > 0:
        ns_dev = np.abs(ns['E_normalized'] - 1.0)
        ns_good = np.sum(ns_dev < 0.05) / len(ns) * 100
        scores['neutron_star_energy'] = ns_good
        print(f"  Neutron Stars ({len(ns)} Objekte):")
        print(f"    < 5% Abweichung:     {ns_good:.1f}%")
        print(f"    Median Abweichung:   {ns_dev.median():.3e}")
    
    # Gesamt
    all_good_energy = np.sum(deviation < 0.01) / len(success) * 100
    scores['overall_energy'] = all_good_energy
    print(f"  GESAMT ({len(success)} Objekte):")
    print(f"    < 1% Abweichung:     {all_good_energy:.1f}%")
    
    # 2. NUMERISCHE STABILITÄT
    print("\n2. NUMERISCHE STABILITÄT:")
    print("-" * 80)
    
    # Keine NaN, keine Inf
    no_nan = np.sum(~np.isnan(success['E_total_J'])) / len(success) * 100
    no_inf = np.sum(~np.isinf(success['E_total_J'])) / len(success) * 100
    
    scores['no_nan'] = no_nan
    scores['no_inf'] = no_inf
    
    print(f"  Keine NaN:           {no_nan:.1f}%")
    print(f"  Keine Inf:           {no_inf:.1f}%")
    
    # Konvergenz
    convergent = np.sum(success['E_normalized'] < 10) / len(success) * 100
    scores['convergent'] = convergent
    print(f"  Konvergent:          {convergent:.1f}%")
    
    # 3. LORENTZ-FAKTOREN
    print("\n3. LORENTZ-FAKTOREN (gamma_gr >= 1):")
    print("-" * 80)
    
    # gamma_gr muss >= 1 sein
    valid_gamma = np.sum(success['gamma_gr_min'] >= 1.0) / len(success) * 100
    scores['valid_gamma'] = valid_gamma
    print(f"  gamma_gr >= 1:       {valid_gamma:.1f}%")
    
    # Erwartete Bereiche
    # Normal stars: gamma ≈ 1 (< 1.001)
    # White dwarfs: gamma < 1.01
    # Neutron stars: gamma > 1.1 (starkes Feld)
    
    ms_gamma_ok = np.sum(ms['gamma_gr_max'] < 1.001) / len(ms) * 100 if len(ms) > 0 else 0
    wd_gamma_ok = np.sum((wd['gamma_gr_max'] > 1.0001) & (wd['gamma_gr_max'] < 1.01)) / len(wd) * 100 if len(wd) > 0 else 0
    ns_gamma_ok = np.sum(ns['gamma_gr_max'] > 1.1) / len(ns) * 100 if len(ns) > 0 else 0
    
    scores['main_sequence_gamma'] = ms_gamma_ok
    scores['white_dwarf_gamma'] = wd_gamma_ok
    scores['neutron_star_gamma'] = ns_gamma_ok
    
    print(f"  Main Sequence (1 < gamma < 1.001):    {ms_gamma_ok:.1f}%")
    print(f"  White Dwarfs (1.0001 < gamma < 1.01): {wd_gamma_ok:.1f}%")
    print(f"  Neutron Stars (gamma > 1.1):          {ns_gamma_ok:.1f}%")
    
    # 4. REDSHIFT-WERTE
    print("\n4. GRAVITATIVER REDSHIFT:")
    print("-" * 80)
    
    # Erwartete Bereiche:
    # Main sequence: z < 10^-5
    # White dwarfs: 10^-5 < z < 10^-3
    # Neutron stars: z > 0.1
    
    ms_z_ok = np.sum(ms['z_gr_max'] < 1e-5) / len(ms) * 100 if len(ms) > 0 else 0
    wd_z_ok = np.sum((wd['z_gr_max'] > 1e-5) & (wd['z_gr_max'] < 1e-3)) / len(wd) * 100 if len(wd) > 0 else 0
    ns_z_ok = np.sum(ns['z_gr_max'] > 0.1) / len(ns) * 100 if len(ns) > 0 else 0
    
    scores['main_sequence_redshift'] = ms_z_ok
    scores['white_dwarf_redshift'] = wd_z_ok
    scores['neutron_star_redshift'] = ns_z_ok
    
    print(f"  Main Sequence (z < 10^-5):            {ms_z_ok:.1f}%")
    print(f"  White Dwarfs (10^-5 < z < 10^-3):     {wd_z_ok:.1f}%")
    print(f"  Neutron Stars (z > 0.1):              {ns_z_ok:.1f}%")
    
    # 5. SHAPIRO DELAY
    print("\n5. SHAPIRO DELAY:")
    print("-" * 80)
    
    valid_shapiro = success[~success['shapiro_delay_s'].isna()]
    if len(valid_shapiro) > 0:
        # Erwartung: 10 us bis 1 ms fuer normale Sterne
        reasonable_shapiro = np.sum((valid_shapiro['shapiro_delay_s'] > 1e-6) & 
                                     (valid_shapiro['shapiro_delay_s'] < 1e-2)) / len(valid_shapiro) * 100
        scores['reasonable_shapiro'] = reasonable_shapiro
        print(f"  Vernuenftiger Bereich (1us - 10ms):   {reasonable_shapiro:.1f}%")
    
    # 6. BEKANNTE OBSERVABLEN
    print("\n6. VERGLEICH MIT BEKANNTEN MESSUNGEN:")
    print("-" * 80)
    
    known_matches = 0
    known_total = 0
    
    # GPS (Erde, falls im Datensatz)
    # Erwartung: ~38 us/Tag Zeitdilatation
    
    # Pound-Rebka (Erde)
    # Erwartung: Delta-gamma/gamma ~ g*h/c^2 ~ 2.46e-15 fuer h=22.5m
    
    # Sirius B
    sirius_b = success[success['name'] == 'Sirius B']
    if len(sirius_b) > 0:
        z_sirius = sirius_b['z_gr_max'].values[0]
        # Gemessen: ~2.7e-4 (inkl. Doppler)
        # Berechnet sollte: ~3.0e-4 sein (nur gravitativ)
        if 2.5e-4 < z_sirius < 3.5e-4:
            known_matches += 1
            print(f"  Sirius B Redshift:   MATCH ({z_sirius:.3e} vs. ~3e-4 erwartet)")
        else:
            print(f"  Sirius B Redshift:   DIFF ({z_sirius:.3e} vs. ~3e-4 erwartet)")
        known_total += 1
    
    # PSR J0740+6620 (NICER)
    psr_j0740 = success[success['name'] == 'PSR J0740+6620']
    if len(psr_j0740) > 0:
        gamma_psr = psr_j0740['gamma_gr_max'].values[0]
        # Erwartung: gamma > 1.3 für R ≈ 2 r_s
        if gamma_psr > 1.3:
            known_matches += 1
            print(f"  PSR J0740 gamma_gr:  MATCH ({gamma_psr:.3f} > 1.3)")
        else:
            print(f"  PSR J0740 gamma_gr:  LOW ({gamma_psr:.3f} < 1.3)")
        known_total += 1
    
    if known_total > 0:
        known_match_pct = known_matches / known_total * 100
        scores['known_observables'] = known_match_pct
        print(f"\n  Bekannte Messungen:  {known_match_pct:.1f}% Match ({known_matches}/{known_total})")
    
    # 7. GESAMT-SCORE
    print("\n" + "="*80)
    print("GESAMT-ÜBEREINSTIMMUNG:")
    print("="*80)
    
    # Gewichteter Durchschnitt
    weights = {
        'overall_energy': 0.30,      # 30% - Wichtigste Größe
        'no_nan': 0.10,               # 10% - Numerisch stabil
        'valid_gamma': 0.15,          # 15% - Physikalisch korrekt
        'main_sequence_gamma': 0.10,  # 10% - Bereiche korrekt
        'white_dwarf_gamma': 0.05,
        'neutron_star_gamma': 0.05,
        'main_sequence_redshift': 0.10,
        'white_dwarf_redshift': 0.05,
        'neutron_star_redshift': 0.05,
        'known_observables': 0.05 if 'known_observables' in scores else 0,
    }
    
    total_score = 0
    total_weight = 0
    
    for key, weight in weights.items():
        if key in scores:
            total_score += scores[key] * weight
            total_weight += weight
    
    if total_weight > 0:
        final_score = total_score / total_weight
    else:
        final_score = 0
    
    scores['FINAL_SCORE'] = final_score
    
    print(f"\n  Energie-Erhaltung:        {scores.get('overall_energy', 0):.1f}% (Gewicht: 30%)")
    print(f"  Numerische Stabilität:    {scores.get('no_nan', 0):.1f}% (Gewicht: 10%)")
    print(f"  Lorentz-Faktoren:         {scores.get('valid_gamma', 0):.1f}% (Gewicht: 15%)")
    
    bereichs_score = np.mean([scores.get('main_sequence_gamma', 0),
                               scores.get('white_dwarf_gamma', 0),
                               scores.get('neutron_star_gamma', 0)])
    print(f"  Bereichs-Uebereinstimmung: {bereichs_score:.1f}% (Gewicht: 20%)")
    
    redshift_score = np.mean([scores.get('main_sequence_redshift', 0),
                               scores.get('white_dwarf_redshift', 0),
                               scores.get('neutron_star_redshift', 0)])
    print(f"  Redshift-Bereiche:        {redshift_score:.1f}% (Gewicht: 20%)")
    if 'known_observables' in scores:
        print(f"  Bekannte Messungen:       {scores['known_observables']:.1f}% (Gewicht: 5%)")
    
    print(f"\n  {'='*40}")
    print(f"  GESAMT-SCORE:             {final_score:.1f}%")
    print(f"  {'='*40}")
    
    # Bewertung
    if final_score >= 95:
        rating = "EXZELLENT"
        emoji = "[+++]"
    elif final_score >= 90:
        rating = "SEHR GUT"
        emoji = "[++]"
    elif final_score >= 80:
        rating = "GUT"
        emoji = "[+]"
    elif final_score >= 70:
        rating = "BEFRIEDIGEND"
        emoji = "[~]"
    else:
        rating = "VERBESSERUNGSBEDARF"
        emoji = "[-]"
    
    print(f"\n  Bewertung: {rating} {emoji}")
    print("="*80)
    
    return scores


# =============================================================================
# Visualisierung
# =============================================================================

def plot_results(results_df, save=True):
    """Visualisiere Ergebnisse"""
    
    success = results_df[results_df['success'] == True]
    
    if len(success) == 0:
        print("Keine erfolgreichen Ergebnisse zum Plotten")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: E_normalized vs Mass
    ax = axes[0, 0]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['mass_Msun'], data['E_normalized'], 
                   label=cat, alpha=0.7, s=100)
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Masse [M_sun]')
    ax.set_ylabel('E_total / E_rest')
    ax.set_title('Energie-Normalisierung vs. Masse')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Teleskopische Differenz
    ax = axes[0, 1]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['r_over_rs'], data['tele_diff'],
                   label=cat, alpha=0.7, s=100)
    ax.set_xlabel('R / r_s')
    ax.set_ylabel('Teleskopische Differenz')
    ax.set_title('Numerische Validierung')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Lorentz-Faktoren
    ax = axes[0, 2]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['r_over_rs'], data['gamma_gr_max'],
                   label=cat, alpha=0.7, s=100)
    ax.set_xlabel('R / r_s')
    ax.set_ylabel('γ_gr (max)')
    ax.set_title('Maximale Zeitdilatation')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Redshift
    ax = axes[1, 0]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['r_over_rs'], data['z_gr_max'],
                   label=cat, alpha=0.7, s=100)
    ax.set_xlabel('R / r_s')
    ax.set_ylabel('z_gr (max)')
    ax.set_title('Gravitativer Redshift')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 5: Shapiro Delay
    ax = axes[1, 1]
    valid = success[~success['shapiro_delay_s'].isna()]
    if len(valid) > 0:
        for cat in valid['category'].unique():
            data = valid[valid['category'] == cat]
            ax.scatter(data['mass_Msun'], data['shapiro_delay_s'],
                       label=cat, alpha=0.7, s=100)
        ax.set_xlabel('Masse [M_sun]')
        ax.set_ylabel('Shapiro Delay [s]')
        ax.set_title('Lichtlaufzeit-Verzögerung')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Plot 6: Energiebeiträge
    ax = axes[1, 2]
    E_GR_rel = np.abs(success['E_GR_J'] / success['E_rest_J'])
    E_SR_rel = success['E_SR_J'] / success['E_rest_J']
    
    ax.scatter(success['r_over_rs'], E_GR_rel, 
               label='|E_GR| / E_rest', alpha=0.7, s=100, marker='o')
    ax.scatter(success['r_over_rs'], E_SR_rel,
               label='E_SR / E_rest', alpha=0.7, s=100, marker='s')
    ax.set_xlabel('R / r_s')
    ax.set_ylabel('Relative Energie')
    ax.set_title('Energie-Beiträge')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig('complete_dataset_results.png', dpi=150, bbox_inches='tight')
        print("\nPlot gespeichert: complete_dataset_results.png")
    
    plt.show()


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    
    print("="*80)
    print("TEST UNIFIED MODEL AUF KOMPLETTEM OBSERVER-DATENSATZ")
    print("="*80)
    
    # Lade Datensatz
    print("\n1. LADE DATENSATZ...")
    # Versuche zuerst großen Datensatz
    try:
        df = pd.read_csv('observer_data_large.csv')
        print(f"Großer Datensatz geladen: {len(df)} Objekte")
    except FileNotFoundError:
        print("Großer Datensatz nicht gefunden, verwende Standard-Datensatz")
        df = fetch_complete_dataset(verbose=True)
    
    # Teste auf allen Objekten
    print("\n2. TESTE AUF ALLEN OBJEKTEN...")
    results_df = test_complete_dataset(
        df, 
        N=1000,
        segmentation='logarithmic',
        verbose_each=False
    )
    
    # Speichere Ergebnisse
    results_df.to_csv('test_results_complete.csv', index=False)
    print("\nErgebnisse gespeichert: test_results_complete.csv")
    
    # Analysiere
    print("\n3. ANALYSIERE ERGEBNISSE...")
    success = analyze_results(results_df)
    
    # Observable Matching
    print("\n4. BERECHNE ÜBEREINSTIMMUNG...")
    scores = calculate_observable_matching(results_df)
    
    # Visualisiere
    print("\n5. ERSTELLE PLOTS...")
    plot_results(results_df, save=True)
    
    print("\n" + "="*80)
    print("FERTIG!")
    print("="*80)
    
    # Finale Zusammenfassung mit Score
    if 'FINAL_SCORE' in scores:
        print(f"\n" + "=" * 80)
        print(f">>> FINALE ÜBEREINSTIMMUNG: {scores['FINAL_SCORE']:.1f}% <<<")
        print("=" * 80)
    print(f"\nZUSAMMENFASSUNG:")
    print(f"  Getestete Objekte: {len(df)}")
    print(f"  Erfolgreich: {results_df['success'].sum()}")
    print(f"  Dateien erstellt:")
    print(f"    - observer_data_complete.csv (Datensatz)")
    print(f"    - test_results_complete.csv (Ergebnisse)")
    print(f"    - complete_dataset_results.png (Plots)")
    print("="*80)
