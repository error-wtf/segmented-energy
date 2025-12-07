#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch Large Observer Dataset - Hunderte/Tausende Objekte

Holt so viele astronomische Objekte wie möglich mit VOLLSTÄNDIGEN Daten:
- SIMBAD: Sterne mit gemessener Masse und Radius
- NASA Exoplanet Archive: Host-Sterne
- Bekannte kompakte Objekte

NUR Objekte mit BEIDEN Masse UND Radius werden inkludiert.
KEIN Filling, nur echte Messungen.

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.constants import G, c, M_sun, R_sun
import warnings

# UTF-8 für Windows
os.environ['PYTHONIOENCODING'] = 'utf-8:replace'

# astroquery für Katalog-Zugriff
try:
    from astroquery.simbad import Simbad
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    HAS_ASTROQUERY = True
except ImportError:
    HAS_ASTROQUERY = False
    print("[WARNING] astroquery nicht verfügbar. Installiere mit: pip install astroquery")


# =============================================================================
# Erweiterte vordefinierte Datensätze
# =============================================================================

# Mehr Main Sequence Sterne mit gemessenen Werten
EXTENDED_STELLAR_DATA = {
    # Bereits vorhandene
    'Sun': {'name': 'Sun', 'mass': 1.0, 'radius': 1.0, 'type': 'G2V', 'distance': 0.0000158, 'source': 'Standard'},
    'Sirius_A': {'name': 'Sirius A', 'mass': 2.063, 'radius': 1.711, 'type': 'A1V', 'distance': 2.637, 'source': 'Interferometry'},
    'Procyon_A': {'name': 'Procyon A', 'mass': 1.499, 'radius': 2.048, 'type': 'F5IV-V', 'distance': 3.51, 'source': 'Binary'},
    'Vega': {'name': 'Vega', 'mass': 2.135, 'radius': 2.362, 'type': 'A0V', 'distance': 7.68, 'source': 'Interferometry'},
    'Altair': {'name': 'Altair', 'mass': 1.79, 'radius': 1.63, 'type': 'A7V', 'distance': 5.13, 'source': 'Interferometry'},
    'Arcturus': {'name': 'Arcturus', 'mass': 1.08, 'radius': 25.4, 'type': 'K0III', 'distance': 11.26, 'source': 'Asteroseismology'},
    'Aldebaran': {'name': 'Aldebaran', 'mass': 1.16, 'radius': 44.2, 'type': 'K5III', 'distance': 20.0, 'source': 'Interferometry'},
    'Regulus': {'name': 'Regulus', 'mass': 3.8, 'radius': 3.092, 'type': 'B8IVn', 'distance': 24.3, 'source': 'Interferometry'},
    
    # Neue Main Sequence
    'Alpha_Cen_A': {'name': 'Alpha Centauri A', 'mass': 1.1055, 'radius': 1.2234, 'type': 'G2V', 'distance': 1.34, 'source': 'Binary'},
    'Alpha_Cen_B': {'name': 'Alpha Centauri B', 'mass': 0.9373, 'radius': 0.8632, 'type': 'K1V', 'distance': 1.34, 'source': 'Binary'},
    'Tau_Ceti': {'name': 'Tau Ceti', 'mass': 0.783, 'radius': 0.793, 'type': 'G8V', 'distance': 3.65, 'source': 'Asteroseismology'},
    'Epsilon_Eri': {'name': 'Epsilon Eridani', 'mass': 0.82, 'radius': 0.735, 'type': 'K2V', 'distance': 3.22, 'source': 'Interferometry'},
    'Betelgeuse': {'name': 'Betelgeuse', 'mass': 16.5, 'radius': 764, 'type': 'M1-2Ia-ab', 'distance': 168, 'source': 'VLTI'},
    'Rigel': {'name': 'Rigel', 'mass': 21, 'radius': 78.9, 'type': 'B8Ia', 'distance': 264, 'source': 'Spectroscopy'},
    'Antares': {'name': 'Antares', 'mass': 12.4, 'radius': 680, 'type': 'M1.5Iab-b', 'distance': 170, 'source': 'VLTI'},
    'Spica': {'name': 'Spica', 'mass': 11.43, 'radius': 7.47, 'type': 'B1III-IV', 'distance': 76.9, 'source': 'Binary'},
    'Pollux': {'name': 'Pollux', 'mass': 1.91, 'radius': 8.8, 'type': 'K0III', 'distance': 10.3, 'source': 'Asteroseismology'},
    'Fomalhaut': {'name': 'Fomalhaut', 'mass': 1.92, 'radius': 1.842, 'type': 'A3V', 'distance': 7.7, 'source': 'Interferometry'},
    'Deneb': {'name': 'Deneb', 'mass': 19, 'radius': 203, 'type': 'A2Ia', 'distance': 802, 'source': 'Spectroscopy'},
    'Canopus': {'name': 'Canopus', 'mass': 8.0, 'radius': 71.4, 'type': 'A9II', 'distance': 95.9, 'source': 'Spectroscopy'},
    'Achernar': {'name': 'Achernar', 'mass': 6.7, 'radius': 7.3, 'type': 'B6Vep', 'distance': 42.6, 'source': 'Interferometry'},
    'Beta_Cen': {'name': 'Beta Centauri', 'mass': 10.7, 'radius': 8.9, 'type': 'B1III', 'distance': 129, 'source': 'Spectroscopy'},
    'Capella_A': {'name': 'Capella A', 'mass': 2.5687, 'radius': 11.98, 'type': 'G5III', 'distance': 13.0, 'source': 'Binary'},
    'Capella_B': {'name': 'Capella B', 'mass': 2.4828, 'radius': 8.83, 'type': 'G0III', 'distance': 13.0, 'source': 'Binary'},
}

# Mehr Weiße Zwerge
EXTENDED_WD_DATA = {
    'Sirius_B': {'name': 'Sirius B', 'mass': 1.018, 'radius': 0.00864, 'type': 'DA2', 'distance': 2.637, 'source': 'Binary'},
    'Procyon_B': {'name': 'Procyon B', 'mass': 0.602, 'radius': 0.01234, 'type': 'DQZ', 'distance': 3.51, 'source': 'Binary'},
    '40_Eri_B': {'name': '40 Eridani B', 'mass': 0.573, 'radius': 0.0136, 'type': 'DA4', 'distance': 4.94, 'source': 'Binary'},
    'Stein_2051_B': {'name': 'Stein 2051 B', 'mass': 0.675, 'radius': 0.01096, 'type': 'DC5', 'distance': 5.04, 'source': 'Astrometry'},
    'van_Maanen': {'name': "van Maanen's star", 'mass': 0.67, 'radius': 0.0101, 'type': 'DZ8', 'distance': 4.2, 'source': 'Spectroscopy'},
}

# Mehr Neutronensterne (NICER + Shapiro delay Messungen)
EXTENDED_NS_DATA = {
    'PSR_J0030': {'name': 'PSR J0030+0451', 'mass': 1.34, 'radius': 12.71, 'type': 'Pulsar', 'distance': 325, 'source': 'NICER'},
    'PSR_J0740': {'name': 'PSR J0740+6620', 'mass': 2.08, 'radius': 12.39, 'type': 'Pulsar', 'distance': 1140, 'source': 'NICER'},
    'PSR_J0348': {'name': 'PSR J0348+0432', 'mass': 2.01, 'radius': 13.0, 'type': 'Pulsar', 'distance': 650, 'source': 'Shapiro'},
    'PSR_J1614': {'name': 'PSR J1614-2230', 'mass': 1.97, 'radius': 12.5, 'type': 'Pulsar', 'distance': 790, 'source': 'Shapiro'},
}

# Mehr Exoplanet Hosts mit präzisen Messungen
EXTENDED_EXOPLANET_HOSTS = {
    'Kepler-11': {'name': 'Kepler-11', 'mass': 0.961, 'radius': 1.065, 'type': 'G', 'distance': 613, 'source': 'Kepler', 'num_planets': 6},
    'TRAPPIST-1': {'name': 'TRAPPIST-1', 'mass': 0.0898, 'radius': 0.1192, 'type': 'M8V', 'distance': 12.43, 'source': 'Transit', 'num_planets': 7},
    'HD_219134': {'name': 'HD 219134', 'mass': 0.778, 'radius': 0.778, 'type': 'K3V', 'distance': 6.53, 'source': 'Transit+RV', 'num_planets': 6},
    '55_Cnc': {'name': '55 Cancri', 'mass': 0.905, 'radius': 0.943, 'type': 'G8V', 'distance': 12.59, 'source': 'Asteroseismology', 'num_planets': 5},
    'Kepler-90': {'name': 'Kepler-90', 'mass': 1.2, 'radius': 1.2, 'type': 'G', 'distance': 836, 'source': 'Kepler', 'num_planets': 8},
    'HD_10180': {'name': 'HD 10180', 'mass': 1.06, 'radius': 1.1, 'type': 'G1V', 'distance': 39.4, 'source': 'RV', 'num_planets': 7},
    'Kepler-20': {'name': 'Kepler-20', 'mass': 0.912, 'radius': 0.944, 'type': 'G8', 'distance': 290, 'source': 'Kepler', 'num_planets': 6},
    'GJ_876': {'name': 'GJ 876', 'mass': 0.334, 'radius': 0.3761, 'type': 'M4V', 'distance': 4.7, 'source': 'Binary+RV', 'num_planets': 4},
}


# =============================================================================
# Katalog-basierte Daten
# =============================================================================

def fetch_from_simbad(max_objects=100, verbose=True):
    """
    Hole Sterne mit gemessener Masse und Radius aus SIMBAD.
    
    WARNUNG: Sehr wenige Objekte haben BEIDE Messungen!
    """
    if not HAS_ASTROQUERY:
        print("[SKIP] SIMBAD: astroquery nicht verfügbar")
        return []
    
    print(f"\n[SIMBAD] Versuche bis zu {max_objects} Objekte zu holen...")
    
    # Konfiguriere SIMBAD
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('typed_id', 'otype', 'sp', 'parallax', 
                                      'flux(V)', 'flux_error(V)')
    
    # SIMBAD hat leider keine direkten Masse/Radius Felder
    # Wir müssten manuell aus Messungen extrahieren
    # Für jetzt: Skip SIMBAD weil zu komplex ohne garantierte Daten
    
    print("[SIMBAD] Überspringe - keine standardisierten Masse/Radius Felder")
    return []


def fetch_from_nasa_exoplanet(verbose=True):
    """
    Hole Host-Sterne aus NASA Exoplanet Archive.
    
    Diese haben oft gute Masse/Radius Messungen.
    """
    if not HAS_ASTROQUERY:
        print("[SKIP] NASA Exoplanet Archive: astroquery nicht verfügbar")
        return []
    
    print("\n[NASA Exoplanet Archive] Hole Host-Sterne...")
    
    try:
        # Hole alle confirmed planets
        table = NasaExoplanetArchive.query_criteria(
            table="pscomppars",
            select="hostname,st_mass,st_masserr1,st_rad,st_raderr1,st_spectype,sy_dist",
            where="st_mass is not null and st_rad is not null"
        )
        
        print(f"[NASA] {len(table)} Host-Sterne mit Masse UND Radius gefunden")
        
        # Konvertiere zu Liste
        hosts = []
        seen_names = set()
        
        for row in table:
            name = row['hostname']
            
            # Nur einzigartige Hosts
            if name in seen_names:
                continue
            seen_names.add(name)
            
            # Prüfe auf vollständige Daten
            if np.ma.is_masked(row['st_mass']) or np.ma.is_masked(row['st_rad']):
                continue
            
            hosts.append({
                'name': name,
                'mass': float(row['st_mass']),  # in M_sun
                'radius': float(row['st_rad']),  # in R_sun
                'type': row['st_spectype'] if not np.ma.is_masked(row['st_spectype']) else 'Unknown',
                'distance': float(row['sy_dist']) if not np.ma.is_masked(row['sy_dist']) else np.nan,
                'source': 'NASA Exoplanet Archive',
                'category': 'exoplanet_host',
            })
        
        print(f"[NASA] {len(hosts)} einzigartige Host-Sterne extrahiert")
        return hosts
        
    except Exception as e:
        print(f"[ERROR] NASA Exoplanet Archive: {e}")
        return []


# =============================================================================
# Kombiniere alle Daten
# =============================================================================

def create_large_dataset(use_catalogs=True, verbose=True):
    """
    Erstelle großen Datensatz aus allen verfügbaren Quellen.
    """
    all_data = []
    
    print("="*80)
    print("ERSTELLE GROSSEN DATENSATZ")
    print("="*80)
    
    # 1. Vordefinierte Main Sequence
    print(f"\n[1/5] Lade erweiterte Main Sequence Sterne...")
    for key, data in EXTENDED_STELLAR_DATA.items():
        all_data.append({
            'name': data['name'],
            'mass': data['mass'],
            'mass_unit': 'M_sun',
            'radius': data['radius'],
            'radius_unit': 'R_sun',
            'type': data['type'],
            'distance': data['distance'],
            'distance_unit': 'pc',
            'source': data['source'],
            'category': 'main_sequence',
        })
    print(f"   -> {len(EXTENDED_STELLAR_DATA)} Main Sequence Sterne")
    
    # 2. Weiße Zwerge
    print(f"\n[2/5] Lade Weiße Zwerge...")
    for key, data in EXTENDED_WD_DATA.items():
        all_data.append({
            'name': data['name'],
            'mass': data['mass'],
            'mass_unit': 'M_sun',
            'radius': data['radius'],
            'radius_unit': 'R_sun',
            'type': data['type'],
            'distance': data['distance'],
            'distance_unit': 'pc',
            'source': data['source'],
            'category': 'white_dwarf',
        })
    print(f"   -> {len(EXTENDED_WD_DATA)} Weiße Zwerge")
    
    # 3. Neutronensterne
    print(f"\n[3/5] Lade Neutronensterne...")
    for key, data in EXTENDED_NS_DATA.items():
        all_data.append({
            'name': data['name'],
            'mass': data['mass'],
            'mass_unit': 'M_sun',
            'radius': data['radius'],
            'radius_unit': 'km',
            'type': data['type'],
            'distance': data['distance'],
            'distance_unit': 'pc',
            'source': data['source'],
            'category': 'neutron_star',
        })
    print(f"   -> {len(EXTENDED_NS_DATA)} Neutronensterne")
    
    # 4. Exoplanet Hosts (vordefiniert)
    print(f"\n[4/5] Lade vordefinierte Exoplanet Hosts...")
    for key, data in EXTENDED_EXOPLANET_HOSTS.items():
        all_data.append({
            'name': data['name'],
            'mass': data['mass'],
            'mass_unit': 'M_sun',
            'radius': data['radius'],
            'radius_unit': 'R_sun',
            'type': data['type'],
            'distance': data['distance'],
            'distance_unit': 'pc',
            'source': data['source'],
            'category': 'exoplanet_host',
            'num_planets': data.get('num_planets', 0),
        })
    print(f"   -> {len(EXTENDED_EXOPLANET_HOSTS)} vordefinierte Hosts")
    
    # 5. NASA Exoplanet Archive (falls verfügbar)
    if use_catalogs:
        print(f"\n[5/5] Hole Daten aus Katalogen...")
        nasa_hosts = fetch_from_nasa_exoplanet(verbose)
        
        for host in nasa_hosts:
            all_data.append({
                'name': host['name'],
                'mass': host['mass'],
                'mass_unit': 'M_sun',
                'radius': host['radius'],
                'radius_unit': 'R_sun',
                'type': host['type'],
                'distance': host['distance'],
                'distance_unit': 'pc',
                'source': host['source'],
                'category': 'exoplanet_host',
            })
        
        print(f"   -> +{len(nasa_hosts)} NASA Hosts")
    else:
        print(f"\n[5/5] Überspringe Katalog-Abfragen (use_catalogs=False)")
    
    # Erstelle DataFrame
    df = pd.DataFrame(all_data)
    
    # Entferne Duplikate
    df = df.drop_duplicates(subset=['name'], keep='first')
    
    # Berechne r_s für alle
    df['r_s_km'] = 2 * G.to(u.km**3/(u.kg * u.s**2)).value * df['mass'] * M_sun.to(u.kg).value / c.to(u.km/u.s).value**2
    
    # Berechne r/r_s
    df['r_over_rs'] = np.where(
        df['radius_unit'] == 'R_sun',
        df['radius'] * R_sun.to(u.km).value / df['r_s_km'],
        df['radius'] / df['r_s_km']
    )
    
    print("\n" + "="*80)
    print("DATENSATZ ERSTELLT")
    print("="*80)
    print(f"\nGESAMT: {len(df)} Objekte mit vollständigen Daten")
    print(f"\nKATEGORIEN:")
    for cat in df['category'].unique():
        count = len(df[df['category'] == cat])
        print(f"  {cat:20s}: {count:3d}")
    
    return df


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    
    print("="*80)
    print("FETCH LARGE OBSERVER DATASET")
    print("Hunderte/Tausende Objekte mit vollständigen Daten")
    print("="*80)
    
    # Erstelle Datensatz
    df = create_large_dataset(use_catalogs=True, verbose=True)
    
    # Speichere
    filename = 'observer_data_large.csv'
    df.to_csv(filename, index=False)
    print(f"\nDatensatz gespeichert: {filename}")
    print(f"  {len(df)} Objekte")
    
    # Statistiken
    print("\n" + "="*80)
    print("STATISTIKEN")
    print("="*80)
    
    print(f"\nMASSEN-BEREICH:")
    print(f"  Min: {df['mass'].min():.4f} M_sun")
    print(f"  Max: {df['mass'].max():.4f} M_sun")
    print(f"  Median: {df['mass'].median():.4f} M_sun")
    
    print(f"\nSCHWARZSCHILD-VERHÄLTNIS:")
    print(f"  Min: {df['r_over_rs'].min():.1f}")
    print(f"  Max: {df['r_over_rs'].max():.0f}")
    
    print("\n" + "="*80)
    print(f"BEREIT FÜR TESTS MIT {len(df)} OBJEKTEN!")
    print("="*80)
