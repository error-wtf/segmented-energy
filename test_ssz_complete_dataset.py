#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test SSZ Model auf vollständigem Observer-Datensatz

Testet Segmented Spacetime (SSZ) mit Xi(r) auf ALLEN gemessenen Objekten.
Vergleicht SSZ vs. GR direkt.

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

from segmented_energy_ssz import compute_ssz_unified, PHI, XI_MAX_DEFAULT
from segmented_energy_unified import compute_unified_energy
from fetch_observer_data import fetch_complete_dataset


# =============================================================================
# Test SSZ auf komplettem Datensatz
# =============================================================================

def test_ssz_dataset(df, N=1000, segmentation='phi', Xi_max=XI_MAX_DEFAULT, 
                      verbose_each=False):
    """
    Teste SSZ Model auf allen Objekten und vergleiche mit GR.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Datensatz von fetch_observer_data
    N : int
        Anzahl Segmente
    segmentation : str
        "phi" oder "logarithmic"
    Xi_max : float
        Maximale Segment Density
    verbose_each : bool
        Zeige Details für jedes Objekt
        
    Returns
    -------
    results_df : pandas.DataFrame
        Ergebnisse für SSZ und GR
    """
    
    results = []
    
    print("\n" + "="*80)
    print(f"TESTE SSZ MODEL AUF {len(df)} OBJEKTEN")
    print("="*80)
    print(f"Segmentierung: {segmentation}, N = {N}, Xi_max = {Xi_max}")
    print()
    
    start_time = time.time()
    
    for idx, row in df.iterrows():
        name = row['name']
        category = row['category']
        
        # Masse und Radius
        M = row['mass'] * M_sun
        
        if row['radius_unit'] == 'R_sun':
            R = row['radius'] * R_sun
        else:
            R = row['radius'] * u.km
        
        # r_in und r_out
        r_s = row['r_s_km'] * u.km
        
        if category in ['main_sequence', 'exoplanet_host']:
            r_in = R * 1.1
            r_out = 100 * r_in
        elif category == 'white_dwarf':
            r_in = R * 1.05
            r_out = 50 * r_in
        elif category == 'neutron_star':
            r_in = R * 1.02
            r_out = 20 * r_in
        else:
            r_in = R
            r_out = 100 * R
        
        m = 1.0 * u.kg
        
        try:
            # SSZ Berechnung
            result_ssz = compute_ssz_unified(
                M=M, m=m, r_in=r_in, r_out=r_out,
                N=N, segmentation=segmentation,
                Xi_max=Xi_max,
                verbose=verbose_each
            )
            
            # GR Berechnung (zum Vergleich)
            result_gr = compute_unified_energy(
                M=M, m=m, r_in=r_in, r_out=r_out,
                N=N, segmentation="logarithmic",  # GR nutzt log
                validate_teleskopic=False,
                verbose=False
            )
            
            obs_ssz = result_ssz['observables']
            
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
                
                # SSZ
                'E_total_SSZ_J': result_ssz['E_total_SSZ'].value,
                'E_SR_SSZ_J': result_ssz['E_SR_SSZ_total'].value,
                'E_normalized_SSZ': result_ssz['E_normalized_SSZ'],
                'Xi_mean': np.mean(result_ssz['Xi']),
                'Xi_max_val': np.max(result_ssz['Xi']),
                'D_SSZ_min': np.min(result_ssz['D_SSZ']),
                'z_SSZ_max': np.max(obs_ssz['z_SSZ']),
                'gamma_ssz_max': np.max(result_ssz['gamma_ssz']),
                
                # GR (zum Vergleich)
                'E_total_GR_J': result_gr['E_total'].value,
                'E_SR_GR_J': result_gr['E_SR_total'].value,
                'E_normalized_GR': result_gr['E_normalized'],
                'gamma_gr_max': np.max(result_gr['gamma_gr']),
                'z_GR_max': np.max(obs_ssz['z_GR']),
                'D_GR_min': np.min(result_gr['gamma_gr'])**(-1),
                
                # Vergleich
                'SSZ_vs_GR_ratio': (result_ssz['E_total_SSZ']/result_gr['E_total']).decompose().value,
                'Delta_E_percent': 100 * abs((result_ssz['E_total_SSZ'] - result_gr['E_total'])/result_gr['E_total']).decompose().value,
                
                'success': True,
                'error': None,
            })
            
            if not verbose_each:
                delta_e = 100 * abs((result_ssz['E_total_SSZ'] - result_gr['E_total'])/result_gr['E_total']).decompose().value
                print(f"  [{idx+1:2d}/{len(df)}] {name:20s} ... OK (SSZ/GR = {delta_e:+.2f}%)")
            
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
# Analyse SSZ vs GR
# =============================================================================

def analyze_ssz_vs_gr(results_df):
    """Analysiere SSZ vs GR Unterschiede"""
    
    success = results_df[results_df['success'] == True].copy()
    
    print("\n" + "="*80)
    print("SSZ vs GR ANALYSE")
    print("="*80)
    
    print(f"\nERFOLGSRATE:")
    print(f"  Erfolgreich: {len(success)}/{len(results_df)} ({100*len(success)/len(results_df):.1f}%)")
    
    if len(success) == 0:
        print("\nKeine erfolgreichen Tests!")
        return
    
    print(f"\nENERGIE-DIFFERENZ (SSZ vs GR):")
    print(f"  Min Delta:    {success['Delta_E_percent'].min():.4f}%")
    print(f"  Max Delta:    {success['Delta_E_percent'].max():.4f}%")
    print(f"  Median Delta: {success['Delta_E_percent'].median():.4f}%")
    
    print(f"\nSSZ SEGMENT DENSITY:")
    print(f"  Xi_mean (min):  {success['Xi_mean'].min():.6f}")
    print(f"  Xi_mean (max):  {success['Xi_mean'].max():.6f}")
    print(f"  Xi_max:         {success['Xi_max_val'].max():.6f}")
    
    print(f"\nZEITDILATATION (D_min = langsamste Zeit):")
    print(f"  SSZ D_min: {success['D_SSZ_min'].min():.6f} ({success.loc[success['D_SSZ_min'].idxmin(), 'name']})")
    print(f"  GR  D_min: {success['D_GR_min'].min():.6f} ({success.loc[success['D_GR_min'].idxmin(), 'name']})")
    
    print(f"\nREDSHIFT (z_max):")
    print(f"  SSZ z_max: {success['z_SSZ_max'].max():.6e} ({success.loc[success['z_SSZ_max'].idxmax(), 'name']})")
    print(f"  GR  z_max: {success['z_GR_max'].max():.6e} ({success.loc[success['z_GR_max'].idxmax(), 'name']})")
    
    print(f"\nPRO KATEGORIE:")
    for cat in success['category'].unique():
        cat_data = success[success['category'] == cat]
        print(f"\n  {cat}:")
        print(f"    Anzahl:           {len(cat_data)}")
        print(f"    Delta E (Median): {cat_data['Delta_E_percent'].median():.4f}%")
        print(f"    Xi_mean:          {cat_data['Xi_mean'].mean():.6f}")
        print(f"    SSZ/GR Ratio:     {cat_data['SSZ_vs_GR_ratio'].mean():.6f}")
    
    print("\n" + "="*80)
    
    return success


# =============================================================================
# SSZ Observable Matching
# =============================================================================

def calculate_ssz_matching(results_df):
    """Berechne prozentuale Übereinstimmung für SSZ"""
    
    success = results_df[results_df['success'] == True]
    
    if len(success) == 0:
        return {}
    
    print("\n" + "="*80)
    print("SSZ OBSERVABLE MATCHING - PROZENTUALE UEBEREINSTIMMUNG")
    print("="*80)
    
    scores = {}
    
    # 1. ENERGIE-ERHALTUNG (SSZ)
    print("\n1. ENERGIE-ERHALTUNG (E_total ~ E_rest):")
    print("-" * 80)
    
    deviation_ssz = np.abs(success['E_normalized_SSZ'] - 1.0)
    deviation_gr = np.abs(success['E_normalized_GR'] - 1.0)
    
    ms = success[success['category'] == 'main_sequence']
    wd = success[success['category'] == 'white_dwarf']
    ns = success[success['category'] == 'neutron_star']
    
    if len(ms) > 0:
        ms_dev_ssz = np.abs(ms['E_normalized_SSZ'] - 1.0)
        ms_good_ssz = np.sum(ms_dev_ssz < 1e-6) / len(ms) * 100
        scores['main_sequence_energy_ssz'] = ms_good_ssz
        print(f"  Main Sequence (SSZ, < 1 ppm):  {ms_good_ssz:.1f}%")
    
    if len(wd) > 0:
        wd_dev_ssz = np.abs(wd['E_normalized_SSZ'] - 1.0)
        wd_good_ssz = np.sum(wd_dev_ssz < 1e-4) / len(wd) * 100
        scores['white_dwarf_energy_ssz'] = wd_good_ssz
        print(f"  White Dwarfs (SSZ, < 0.01%):   {wd_good_ssz:.1f}%")
    
    if len(ns) > 0:
        ns_dev_ssz = np.abs(ns['E_normalized_SSZ'] - 1.0)
        ns_good_ssz = np.sum(ns_dev_ssz < 0.05) / len(ns) * 100
        scores['neutron_star_energy_ssz'] = ns_good_ssz
        print(f"  Neutron Stars (SSZ, < 5%):     {ns_good_ssz:.1f}%")
    
    all_good_ssz = np.sum(deviation_ssz < 0.01) / len(success) * 100
    scores['overall_energy_ssz'] = all_good_ssz
    print(f"  GESAMT (SSZ, < 1%):             {all_good_ssz:.1f}%")
    
    # 2. SSZ vs GR KONSISTENZ
    print("\n2. SSZ vs GR KONSISTENZ:")
    print("-" * 80)
    
    # Für normale Sterne: SSZ sollte nahe GR sein (schwaches Feld)
    # Für kompakte Objekte: SSZ darf abweichen (starkes Feld)
    
    ms_consistent = np.sum(ms['Delta_E_percent'] < 0.1) / len(ms) * 100 if len(ms) > 0 else 0
    wd_allows_diff = np.sum(wd['Delta_E_percent'] < 1.0) / len(wd) * 100 if len(wd) > 0 else 0
    ns_allows_diff = np.sum(ns['Delta_E_percent'] < 5.0) / len(ns) * 100 if len(ns) > 0 else 0
    
    scores['ms_ssz_gr_consistent'] = ms_consistent
    scores['wd_ssz_gr_acceptable'] = wd_allows_diff
    scores['ns_ssz_gr_acceptable'] = ns_allows_diff
    
    print(f"  Main Sequence (Delta < 0.1%):   {ms_consistent:.1f}%")
    print(f"  White Dwarfs (Delta < 1%):      {wd_allows_diff:.1f}%")
    print(f"  Neutron Stars (Delta < 5%):     {ns_allows_diff:.1f}%")
    
    # 3. SEGMENT DENSITY
    print("\n3. SEGMENT DENSITY (Xi):")
    print("-" * 80)
    
    # Xi sollte < Xi_max sein
    valid_xi = np.sum(success['Xi_max_val'] <= XI_MAX_DEFAULT * 1.01) / len(success) * 100
    scores['valid_xi'] = valid_xi
    print(f"  Xi <= Xi_max:                   {valid_xi:.1f}%")
    
    # Xi sollte > 0 sein (Diskretheit vorhanden)
    nonzero_xi = np.sum(success['Xi_mean'] > 0.001) / len(success) * 100
    scores['nonzero_xi'] = nonzero_xi
    print(f"  Xi > 0.001 (diskret):           {nonzero_xi:.1f}%")
    
    # 4. SSZ ZEITDILATATION
    print("\n4. SSZ ZEITDILATATION:")
    print("-" * 80)
    
    # D_SSZ sollte < 1 sein (Zeit läuft langsamer)
    valid_d_ssz = np.sum(success['D_SSZ_min'] < 1.0) / len(success) * 100
    scores['valid_d_ssz'] = valid_d_ssz
    print(f"  D_SSZ < 1 (langsamer):          {valid_d_ssz:.1f}%")
    
    # SSZ sollte singularitätsfrei sein
    finite_d_ssz = np.sum(success['D_SSZ_min'] > 0.1) / len(success) * 100
    scores['finite_d_ssz'] = finite_d_ssz
    print(f"  D_SSZ > 0.1 (singularitätsfrei):{finite_d_ssz:.1f}%")
    
    # 5. GESAMT-SCORE
    print("\n" + "="*80)
    print("GESAMT-UEBEREINSTIMMUNG (SSZ):")
    print("="*80)
    
    weights = {
        'overall_energy_ssz': 0.30,
        'ms_ssz_gr_consistent': 0.20,
        'wd_ssz_gr_acceptable': 0.10,
        'ns_ssz_gr_acceptable': 0.10,
        'valid_xi': 0.10,
        'nonzero_xi': 0.05,
        'valid_d_ssz': 0.10,
        'finite_d_ssz': 0.05,
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
    
    scores['FINAL_SCORE_SSZ'] = final_score
    
    ssz_gr_cons = np.mean([scores.get('ms_ssz_gr_consistent', 0),
                            scores.get('wd_ssz_gr_acceptable', 0),
                            scores.get('ns_ssz_gr_acceptable', 0)])
    seg_dens = np.mean([scores.get('valid_xi', 0),
                         scores.get('nonzero_xi', 0)])
    ssz_time = np.mean([scores.get('valid_d_ssz', 0),
                         scores.get('finite_d_ssz', 0)])
    
    print(f"\n  Energie-Erhaltung (SSZ):  {scores.get('overall_energy_ssz', 0):.1f}% (Gewicht: 30%)")
    print(f"  SSZ/GR Konsistenz:        {ssz_gr_cons:.1f}% (Gewicht: 40%)")
    print(f"  Segment Density:          {seg_dens:.1f}% (Gewicht: 15%)")
    print(f"  SSZ Zeitdilatation:       {ssz_time:.1f}% (Gewicht: 15%)")
    
    print(f"\n  {'='*40}")
    print(f"  GESAMT-SCORE (SSZ):       {final_score:.1f}%")
    print(f"  {'='*40}")
    
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
    
    print(f"\n  Bewertung (SSZ): {rating} {emoji}")
    print("="*80)
    
    return scores


# =============================================================================
# Visualisierung
# =============================================================================

def plot_ssz_vs_gr(results_df, save=True):
    """Visualisiere SSZ vs GR Vergleich"""
    
    success = results_df[results_df['success'] == True]
    
    if len(success) == 0:
        print("Keine Ergebnisse zum Plotten")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: E_normalized SSZ vs GR
    ax = axes[0, 0]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['E_normalized_GR'], data['E_normalized_SSZ'],
                   label=cat, alpha=0.7, s=100)
    ax.plot([0.96, 1.001], [0.96, 1.001], 'k--', alpha=0.3, label='1:1')
    ax.set_xlabel('E_GR / E_rest')
    ax.set_ylabel('E_SSZ / E_rest')
    ax.set_title('Energie-Normalisierung: SSZ vs GR')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Delta E vs R/r_s
    ax = axes[0, 1]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['r_over_rs'], data['Delta_E_percent'],
                   label=cat, alpha=0.7, s=100)
    ax.set_xlabel('R / r_s')
    ax.set_ylabel('|E_SSZ - E_GR| / E_GR [%]')
    ax.set_title('SSZ vs GR Energie-Differenz')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Segment Density
    ax = axes[0, 2]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['r_over_rs'], data['Xi_mean'],
                   label=cat, alpha=0.7, s=100)
    ax.axhline(y=XI_MAX_DEFAULT, color='r', linestyle='--', 
               label=f'Xi_max = {XI_MAX_DEFAULT}')
    ax.set_xlabel('R / r_s')
    ax.set_ylabel('Xi (mean)')
    ax.set_title('Segment Density vs. Kompaktheit')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Zeitdilatation SSZ vs GR
    ax = axes[1, 0]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['D_GR_min'], data['D_SSZ_min'],
                   label=cat, alpha=0.7, s=100)
    ax.plot([0.5, 1], [0.5, 1], 'k--', alpha=0.3, label='1:1')
    ax.set_xlabel('D_GR (min)')
    ax.set_ylabel('D_SSZ (min)')
    ax.set_title('Zeitdilatation: SSZ vs GR')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 5: Redshift SSZ vs GR
    ax = axes[1, 1]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['z_GR_max'], data['z_SSZ_max'],
                   label=cat, alpha=0.7, s=100)
    ax.plot([1e-8, 1], [1e-8, 1], 'k--', alpha=0.3, label='1:1')
    ax.set_xlabel('z_GR (max)')
    ax.set_ylabel('z_SSZ (max)')
    ax.set_title('Redshift: SSZ vs GR')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 6: Gamma SSZ vs GR
    ax = axes[1, 2]
    for cat in success['category'].unique():
        data = success[success['category'] == cat]
        ax.scatter(data['gamma_gr_max'], data['gamma_ssz_max'],
                   label=cat, alpha=0.7, s=100)
    ax.plot([1, 2], [1, 2], 'k--', alpha=0.3, label='1:1')
    ax.set_xlabel('gamma_GR (max)')
    ax.set_ylabel('gamma_SSZ (max)')
    ax.set_title('Lorentz-Faktoren: SSZ vs GR')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig('ssz_vs_gr_complete_dataset.png', dpi=150, bbox_inches='tight')
        print("\nPlot gespeichert: ssz_vs_gr_complete_dataset.png")
    
    plt.show()


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    
    print("="*80)
    print("TEST SSZ MODEL AUF KOMPLETTEM OBSERVER-DATENSATZ")
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
    
    # Teste SSZ
    print("\n2. TESTE SSZ AUF ALLEN OBJEKTEN...")
    results_df = test_ssz_dataset(
        df,
        N=1000,
        segmentation='phi',
        Xi_max=0.8,
        verbose_each=False
    )
    
    # Speichere
    results_df.to_csv('test_results_ssz_complete.csv', index=False)
    print("\nErgebnisse gespeichert: test_results_ssz_complete.csv")
    
    # Analysiere
    print("\n3. ANALYSIERE SSZ vs GR...")
    success = analyze_ssz_vs_gr(results_df)
    
    # Matching
    print("\n4. BERECHNE SSZ-UEBEREINSTIMMUNG...")
    scores = calculate_ssz_matching(results_df)
    
    # Visualisiere
    print("\n5. ERSTELLE PLOTS...")
    plot_ssz_vs_gr(results_df, save=True)
    
    print("\n" + "="*80)
    print("FERTIG!")
    print("="*80)
    
    # Finale Zusammenfassung
    if 'FINAL_SCORE_SSZ' in scores:
        print(f"\n" + "=" * 80)
        print(f">>> FINALE SSZ-UEBEREINSTIMMUNG: {scores['FINAL_SCORE_SSZ']:.1f}% <<<")
        print("=" * 80)
    
    print(f"\nZUSAMMENFASSUNG:")
    print(f"  Getestete Objekte: {len(df)}")
    print(f"  Erfolgreich: {results_df['success'].sum()}")
    print(f"  Dateien erstellt:")
    print(f"    - test_results_ssz_complete.csv (SSZ Ergebnisse)")
    print(f"    - ssz_vs_gr_complete_dataset.png (Vergleichs-Plots)")
    print("="*80)
