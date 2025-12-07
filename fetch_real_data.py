#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch Real Astronomical Data - For Segmented Energy Tests

Holt echte astronomische Daten aus verschiedenen Quellen:
- SIMBAD (Sterndaten)
- NASA Exoplanet Archive
- Bekannte physikalische Systeme

© 2025 Carmen Wrede & Lino Casu
Licensed under the ANTI-CAPITALIST SOFTWARE LICENSE v1.4
"""

import numpy as np
from astropy import units as u
from astropy.constants import G, c, M_sun, M_earth, R_sun, R_earth, au
from astropy.coordinates import SkyCoord
import json

# Optional: astroquery für erweiterte Daten
try:
    from astroquery.simbad import Simbad
    from astroquery.gaia import Gaia
    HAS_ASTROQUERY = True
except ImportError:
    print("⚠️  astroquery nicht installiert. Verwende vordefinierte Daten.")
    HAS_ASTROQUERY = False


# =============================================================================
# Vordefinierte Astronomische Systeme
# =============================================================================

STELLAR_SYSTEMS = {
    # ==========================================================================
    # MAIN SEQUENCE STARS (30+ objects)
    # ==========================================================================
    
    # G-type (Sun-like)
    "sun": {
        "name": "Sun (Sol)",
        "type": "G2V main sequence",
        "mass": 1.0 * M_sun,
        "radius": 1.0 * R_sun,
        "description": "Our Sun",
        "r_in": 2.0 * R_sun,
        "r_out": 1.0 * au,
    },
    "alpha_centauri_a": {
        "name": "Alpha Centauri A",
        "type": "G2V main sequence",
        "mass": 1.1 * M_sun,
        "radius": 1.227 * R_sun,
        "description": "Nearest Sun-like star",
        "r_in": 2.5 * R_sun,
        "r_out": 1.2 * au,
    },
    "tau_ceti": {
        "name": "Tau Ceti",
        "type": "G8.5V main sequence",
        "mass": 0.783 * M_sun,
        "radius": 0.793 * R_sun,
        "description": "Nearby Sun-like star",
        "r_in": 1.6 * R_sun,
        "r_out": 0.8 * au,
    },
    "18_scorpii": {
        "name": "18 Scorpii",
        "type": "G2Va main sequence",
        "mass": 1.02 * M_sun,
        "radius": 1.01 * R_sun,
        "description": "Solar twin",
        "r_in": 2.0 * R_sun,
        "r_out": 1.0 * au,
    },
    
    # A-type (hot white)
    "sirius_a": {
        "name": "Sirius A",
        "type": "A1V main sequence",
        "mass": 2.063 * M_sun,
        "radius": 1.711 * R_sun,
        "description": "Brightest star in night sky",
        "r_in": 3.5 * R_sun,
        "r_out": 1.5 * au,
    },
    "vega": {
        "name": "Vega",
        "type": "A0V main sequence",
        "mass": 2.135 * M_sun,
        "radius": 2.362 * R_sun,
        "description": "Standard reference star",
        "r_in": 4.7 * R_sun,
        "r_out": 2.0 * au,
    },
    "altair": {
        "name": "Altair",
        "type": "A7V main sequence",
        "mass": 1.79 * M_sun,
        "radius": 1.63 * R_sun,
        "description": "Rapidly rotating star",
        "r_in": 3.3 * R_sun,
        "r_out": 1.5 * au,
    },
    "fomalhaut": {
        "name": "Fomalhaut",
        "type": "A3V main sequence",
        "mass": 1.92 * M_sun,
        "radius": 1.842 * R_sun,
        "description": "Star with debris disk",
        "r_in": 3.7 * R_sun,
        "r_out": 1.8 * au,
    },
    "deneb": {
        "name": "Deneb",
        "type": "A2Ia supergiant",
        "mass": 19.0 * M_sun,
        "radius": 203 * R_sun,
        "description": "Luminous supergiant",
        "r_in": 210 * R_sun,
        "r_out": 8.0 * au,
    },
    
    # F-type
    "procyon_a": {
        "name": "Procyon A",
        "type": "F5IV-V subgiant",
        "mass": 1.499 * M_sun,
        "radius": 2.048 * R_sun,
        "description": "Nearby bright star",
        "r_in": 4.1 * R_sun,
        "r_out": 1.5 * au,
    },
    "canopus": {
        "name": "Canopus",
        "type": "F0Ib supergiant",
        "mass": 8.0 * M_sun,
        "radius": 71 * R_sun,
        "description": "Second brightest star",
        "r_in": 75 * R_sun,
        "r_out": 5.0 * au,
    },
    
    # K-type (orange)
    "alpha_centauri_b": {
        "name": "Alpha Centauri B",
        "type": "K1V main sequence",
        "mass": 0.907 * M_sun,
        "radius": 0.865 * R_sun,
        "description": "Nearest K-dwarf",
        "r_in": 1.7 * R_sun,
        "r_out": 0.9 * au,
    },
    "epsilon_eridani": {
        "name": "Epsilon Eridani",
        "type": "K2V main sequence",
        "mass": 0.82 * M_sun,
        "radius": 0.735 * R_sun,
        "description": "Young nearby star",
        "r_in": 1.5 * R_sun,
        "r_out": 0.8 * au,
    },
    "arcturus": {
        "name": "Arcturus",
        "type": "K0III giant",
        "mass": 1.08 * M_sun,
        "radius": 25.4 * R_sun,
        "description": "Brightest in northern sky",
        "r_in": 26 * R_sun,
        "r_out": 2.0 * au,
    },
    "aldebaran": {
        "name": "Aldebaran",
        "type": "K5III giant",
        "mass": 1.16 * M_sun,
        "radius": 44.13 * R_sun,
        "description": "Eye of Taurus",
        "r_in": 45 * R_sun,
        "r_out": 3.0 * au,
    },
    
    # M-type (red dwarfs)
    "proxima_centauri": {
        "name": "Proxima Centauri",
        "type": "M5.5Ve red dwarf",
        "mass": 0.1221 * M_sun,
        "radius": 0.1542 * R_sun,
        "description": "Nearest star to Sun",
        "r_in": 0.3 * R_sun,
        "r_out": 0.1 * au,
    },
    "barnards_star": {
        "name": "Barnard's Star",
        "type": "M4Ve red dwarf",
        "mass": 0.144 * M_sun,
        "radius": 0.196 * R_sun,
        "description": "Second nearest star system",
        "r_in": 0.4 * R_sun,
        "r_out": 0.12 * au,
    },
    "wolf_359": {
        "name": "Wolf 359",
        "type": "M6.5Ve red dwarf",
        "mass": 0.09 * M_sun,
        "radius": 0.16 * R_sun,
        "description": "Nearby flare star",
        "r_in": 0.32 * R_sun,
        "r_out": 0.08 * au,
    },
    "lalande_21185": {
        "name": "Lalande 21185",
        "type": "M2V red dwarf",
        "mass": 0.46 * M_sun,
        "radius": 0.393 * R_sun,
        "description": "Fourth nearest star",
        "r_in": 0.8 * R_sun,
        "r_out": 0.3 * au,
    },
    "ross_128": {
        "name": "Ross 128",
        "type": "M4V red dwarf",
        "mass": 0.168 * M_sun,
        "radius": 0.1967 * R_sun,
        "description": "Quiet red dwarf",
        "r_in": 0.4 * R_sun,
        "r_out": 0.15 * au,
    },
    "gliese_581": {
        "name": "Gliese 581",
        "type": "M3V red dwarf",
        "mass": 0.31 * M_sun,
        "radius": 0.299 * R_sun,
        "description": "Exoplanet host",
        "r_in": 0.6 * R_sun,
        "r_out": 0.2 * au,
    },
    "trappist_1_star": {
        "name": "TRAPPIST-1",
        "type": "M8V ultracool dwarf",
        "mass": 0.0898 * M_sun,
        "radius": 0.1192 * R_sun,
        "description": "7-planet system host",
        "r_in": 0.24 * R_sun,
        "r_out": 0.07 * au,
    },
    
    # B-type (hot blue)
    "spica": {
        "name": "Spica",
        "type": "B1III-IV giant",
        "mass": 11.43 * M_sun,
        "radius": 7.47 * R_sun,
        "description": "Brightest in Virgo",
        "r_in": 15 * R_sun,
        "r_out": 3.0 * au,
    },
    "regulus": {
        "name": "Regulus",
        "type": "B8IVn subgiant",
        "mass": 3.8 * M_sun,
        "radius": 3.092 * R_sun,
        "description": "Heart of Leo",
        "r_in": 6.2 * R_sun,
        "r_out": 2.0 * au,
    },
    "achernar": {
        "name": "Achernar",
        "type": "B6Vep main sequence",
        "mass": 6.7 * M_sun,
        "radius": 7.3 * R_sun,
        "description": "Fastest rotating star",
        "r_in": 15 * R_sun,
        "r_out": 3.0 * au,
    },
    
    # O-type (hottest)
    "zeta_puppis": {
        "name": "Zeta Puppis (Naos)",
        "type": "O4If supergiant",
        "mass": 56.1 * M_sun,
        "radius": 14 * R_sun,
        "description": "Hottest naked-eye star",
        "r_in": 28 * R_sun,
        "r_out": 5.0 * au,
    },
    "theta1_orionis_c": {
        "name": "Theta1 Orionis C",
        "type": "O6Vp main sequence",
        "mass": 33 * M_sun,
        "radius": 10.6 * R_sun,
        "description": "Trapezium cluster star",
        "r_in": 21 * R_sun,
        "r_out": 4.0 * au,
    },
    
    # ==========================================================================
    # GIANTS AND SUPERGIANTS (15+ objects)
    # ==========================================================================
    
    "betelgeuse": {
        "name": "Betelgeuse",
        "type": "M1-2 Ia-ab red supergiant",
        "mass": 16.5 * M_sun,
        "radius": 764 * R_sun,
        "description": "Red supergiant in Orion",
        "r_in": 800 * R_sun,
        "r_out": 10.0 * au,
    },
    "rigel": {
        "name": "Rigel",
        "type": "B8 Ia blue supergiant",
        "mass": 21.0 * M_sun,
        "radius": 78.9 * R_sun,
        "description": "Blue supergiant in Orion",
        "r_in": 100 * R_sun,
        "r_out": 5.0 * au,
    },
    "antares": {
        "name": "Antares",
        "type": "M1.5Iab red supergiant",
        "mass": 12.4 * M_sun,
        "radius": 680 * R_sun,
        "description": "Heart of Scorpius",
        "r_in": 700 * R_sun,
        "r_out": 9.0 * au,
    },
    "uy_scuti": {
        "name": "UY Scuti",
        "type": "M4Ia red supergiant",
        "mass": 7.0 * M_sun,
        "radius": 1708 * R_sun,
        "description": "One of largest known stars",
        "r_in": 1750 * R_sun,
        "r_out": 20.0 * au,
    },
    "vv_cephei_a": {
        "name": "VV Cephei A",
        "type": "M2Iab red supergiant",
        "mass": 18.2 * M_sun,
        "radius": 1050 * R_sun,
        "description": "Eclipsing binary primary",
        "r_in": 1100 * R_sun,
        "r_out": 15.0 * au,
    },
    "mu_cephei": {
        "name": "Mu Cephei (Garnet Star)",
        "type": "M2Ia red supergiant",
        "mass": 19.2 * M_sun,
        "radius": 1420 * R_sun,
        "description": "Herschel's Garnet Star",
        "r_in": 1450 * R_sun,
        "r_out": 18.0 * au,
    },
    "polaris": {
        "name": "Polaris",
        "type": "F7Ib supergiant",
        "mass": 5.4 * M_sun,
        "radius": 37.5 * R_sun,
        "description": "North Star (Cepheid)",
        "r_in": 40 * R_sun,
        "r_out": 3.0 * au,
    },
    "eta_carinae": {
        "name": "Eta Carinae",
        "type": "LBV hypergiant",
        "mass": 100 * M_sun,
        "radius": 240 * R_sun,
        "description": "Luminous Blue Variable",
        "r_in": 250 * R_sun,
        "r_out": 10.0 * au,
    },
    "pistol_star": {
        "name": "Pistol Star",
        "type": "LBV hypergiant",
        "mass": 27.5 * M_sun,
        "radius": 306 * R_sun,
        "description": "One of most luminous",
        "r_in": 320 * R_sun,
        "r_out": 12.0 * au,
    },
    "r136a1": {
        "name": "R136a1",
        "type": "WN5h Wolf-Rayet",
        "mass": 196 * M_sun,
        "radius": 28.8 * R_sun,
        "description": "Most massive known star",
        "r_in": 60 * R_sun,
        "r_out": 8.0 * au,
    },
    
    # ==========================================================================
    # WHITE DWARFS (10 objects)
    # ==========================================================================
    
    "sirius_b": {
        "name": "Sirius B",
        "type": "DA2 white dwarf",
        "mass": 1.018 * M_sun,
        "radius": 0.0084 * R_sun,  # ~5800 km
        "description": "Nearest white dwarf",
        "r_in": 0.02 * R_sun,
        "r_out": 0.1 * R_sun,
    },
    "procyon_b": {
        "name": "Procyon B",
        "type": "DQZ white dwarf",
        "mass": 0.602 * M_sun,
        "radius": 0.0096 * R_sun,
        "description": "Companion to Procyon A",
        "r_in": 0.02 * R_sun,
        "r_out": 0.1 * R_sun,
    },
    "van_maanen_2": {
        "name": "Van Maanen 2",
        "type": "DZ8 white dwarf",
        "mass": 0.68 * M_sun,
        "radius": 0.009 * R_sun,
        "description": "Third nearest WD",
        "r_in": 0.02 * R_sun,
        "r_out": 0.1 * R_sun,
    },
    "40_eridani_b": {
        "name": "40 Eridani B",
        "type": "DA4 white dwarf",
        "mass": 0.573 * M_sun,
        "radius": 0.0136 * R_sun,
        "description": "First WD discovered",
        "r_in": 0.03 * R_sun,
        "r_out": 0.15 * R_sun,
    },
    "stein_2051_b": {
        "name": "Stein 2051 B",
        "type": "DC5 white dwarf",
        "mass": 0.675 * M_sun,
        "radius": 0.0108 * R_sun,
        "description": "Gravitational lensing WD",
        "r_in": 0.02 * R_sun,
        "r_out": 0.1 * R_sun,
    },
    "gd_165_b": {
        "name": "GD 165 B",
        "type": "L4 brown dwarf",
        "mass": 0.072 * M_sun,
        "radius": 0.1 * R_sun,
        "description": "First L-dwarf discovered",
        "r_in": 0.2 * R_sun,
        "r_out": 0.5 * R_sun,
    },
    
    # ==========================================================================
    # NEUTRON STARS (10 objects)
    # ==========================================================================
    
    "neutron_star": {
        "name": "Typical Neutron Star",
        "type": "Neutron star",
        "mass": 1.4 * M_sun,
        "radius": 12.0 * u.km,
        "description": "Standard NS (1.4 solar masses)",
        "r_in": 15.0 * u.km,
        "r_out": 100.0 * u.km,
    },
    "psr_j0348": {
        "name": "PSR J0348+0432",
        "type": "Millisecond pulsar",
        "mass": 2.01 * M_sun,
        "radius": 13.0 * u.km,
        "description": "Massive NS with WD companion",
        "r_in": 16.0 * u.km,
        "r_out": 120.0 * u.km,
    },
    "psr_j0740": {
        "name": "PSR J0740+6620",
        "type": "Millisecond pulsar",
        "mass": 2.08 * M_sun,
        "radius": 12.35 * u.km,
        "description": "Most massive known NS (NICER)",
        "r_in": 15.0 * u.km,
        "r_out": 110.0 * u.km,
    },
    "psr_b1937": {
        "name": "PSR B1937+21",
        "type": "Millisecond pulsar",
        "mass": 1.4 * M_sun,
        "radius": 10.0 * u.km,
        "description": "First MSP discovered",
        "r_in": 12.0 * u.km,
        "r_out": 80.0 * u.km,
    },
    "crab_pulsar": {
        "name": "Crab Pulsar",
        "type": "Young pulsar",
        "mass": 1.4 * M_sun,
        "radius": 10.0 * u.km,
        "description": "In Crab Nebula (SN 1054)",
        "r_in": 12.0 * u.km,
        "r_out": 80.0 * u.km,
    },
    "vela_pulsar": {
        "name": "Vela Pulsar",
        "type": "Young pulsar",
        "mass": 1.4 * M_sun,
        "radius": 11.0 * u.km,
        "description": "In Vela SNR",
        "r_in": 14.0 * u.km,
        "r_out": 90.0 * u.km,
    },
    "magnetar_sgr_1806": {
        "name": "SGR 1806-20",
        "type": "Magnetar",
        "mass": 1.4 * M_sun,
        "radius": 10.0 * u.km,
        "description": "Strongest magnetic field",
        "r_in": 12.0 * u.km,
        "r_out": 80.0 * u.km,
    },
    
    # ==========================================================================
    # BLACK HOLES - STELLAR (10 objects)
    # ==========================================================================
    
    "cygnus_x1": {
        "name": "Cygnus X-1",
        "type": "Stellar black hole",
        "mass": 21.2 * M_sun,
        "radius": None,
        "description": "First confirmed stellar BH",
        "r_in": None,
        "r_out": 50.0 * R_sun,
        "schwarzschild_radius": 62.6 * u.km,
    },
    "grs_1915": {
        "name": "GRS 1915+105",
        "type": "Microquasar BH",
        "mass": 12.4 * M_sun,
        "radius": None,
        "description": "Fastest jets observed",
        "r_in": None,
        "r_out": 30.0 * R_sun,
    },
    "v404_cygni": {
        "name": "V404 Cygni",
        "type": "X-ray binary BH",
        "mass": 9.0 * M_sun,
        "radius": None,
        "description": "Nearest BH (~7800 ly)",
        "r_in": None,
        "r_out": 25.0 * R_sun,
    },
    "a0620_00": {
        "name": "A0620-00",
        "type": "X-ray nova BH",
        "mass": 6.6 * M_sun,
        "radius": None,
        "description": "Nearest confirmed BH",
        "r_in": None,
        "r_out": 20.0 * R_sun,
    },
    "gw150914_primary": {
        "name": "GW150914 Primary",
        "type": "Merging BH",
        "mass": 36.0 * M_sun,
        "radius": None,
        "description": "First GW detection",
        "r_in": None,
        "r_out": 80.0 * R_sun,
    },
    "gw150914_secondary": {
        "name": "GW150914 Secondary",
        "type": "Merging BH",
        "mass": 29.0 * M_sun,
        "radius": None,
        "description": "First GW detection",
        "r_in": None,
        "r_out": 65.0 * R_sun,
    },
    "gw190521_primary": {
        "name": "GW190521 Primary",
        "type": "IMBH merger",
        "mass": 85.0 * M_sun,
        "radius": None,
        "description": "First IMBH detection",
        "r_in": None,
        "r_out": 150.0 * R_sun,
    },
    
    # ==========================================================================
    # BLACK HOLES - SUPERMASSIVE (10 objects)
    # ==========================================================================
    
    "sgr_a_star": {
        "name": "Sagittarius A*",
        "type": "Supermassive black hole",
        "mass": 4.154e6 * M_sun,
        "radius": None,
        "description": "SMBH at Galactic center",
        "r_in": None,
        "r_out": 1000.0 * au,
        "schwarzschild_radius": 1.232e10 * u.m,
    },
    "m87_star": {
        "name": "M87*",
        "type": "Supermassive black hole",
        "mass": 6.5e9 * M_sun,
        "radius": None,
        "description": "First imaged BH (EHT)",
        "r_in": None,
        "r_out": 5000.0 * au,
    },
    "ngc_1277_bh": {
        "name": "NGC 1277 BH",
        "type": "Supermassive black hole",
        "mass": 1.7e10 * M_sun,
        "radius": None,
        "description": "One of largest known",
        "r_in": None,
        "r_out": 10000.0 * au,
    },
    "ton_618": {
        "name": "TON 618",
        "type": "Ultramassive black hole",
        "mass": 6.6e10 * M_sun,
        "radius": None,
        "description": "Most massive known BH",
        "r_in": None,
        "r_out": 50000.0 * au,
    },
    "andromeda_bh": {
        "name": "Andromeda BH (M31*)",
        "type": "Supermassive black hole",
        "mass": 1.4e8 * M_sun,
        "radius": None,
        "description": "SMBH in Andromeda",
        "r_in": None,
        "r_out": 2000.0 * au,
    },
    "ngc_4889_bh": {
        "name": "NGC 4889 BH",
        "type": "Supermassive black hole",
        "mass": 2.1e10 * M_sun,
        "radius": None,
        "description": "In Coma Cluster",
        "r_in": None,
        "r_out": 15000.0 * au,
    },
    "ic_1101_bh": {
        "name": "IC 1101 BH",
        "type": "Supermassive black hole",
        "mass": 4.0e10 * M_sun,
        "radius": None,
        "description": "In largest known galaxy",
        "r_in": None,
        "r_out": 30000.0 * au,
    },
}


EXOPLANET_SYSTEMS = {
    # ==========================================================================
    # MULTI-PLANET SYSTEMS (10+ systems, 50+ planets)
    # ==========================================================================
    
    "solar_system": {
        "name": "Solar System",
        "star_mass": 1.0 * M_sun,
        "star_radius": 1.0 * R_sun,
        "description": "Our home system",
        "planets": [
            {"name": "Mercury", "mass": 0.055 * M_earth, "orbital_radius": 0.387 * au},
            {"name": "Venus", "mass": 0.815 * M_earth, "orbital_radius": 0.723 * au},
            {"name": "Earth", "mass": 1.0 * M_earth, "orbital_radius": 1.0 * au},
            {"name": "Mars", "mass": 0.107 * M_earth, "orbital_radius": 1.524 * au},
            {"name": "Jupiter", "mass": 317.8 * M_earth, "orbital_radius": 5.203 * au},
            {"name": "Saturn", "mass": 95.2 * M_earth, "orbital_radius": 9.537 * au},
            {"name": "Uranus", "mass": 14.5 * M_earth, "orbital_radius": 19.19 * au},
            {"name": "Neptune", "mass": 17.1 * M_earth, "orbital_radius": 30.07 * au},
        ],
    },
    
    "kepler_11": {
        "name": "Kepler-11 system",
        "star_mass": 0.961 * M_sun,
        "star_radius": 1.065 * R_sun,
        "description": "Compact multi-planet system",
        "planets": [
            {"name": "b", "mass": 1.9 * M_earth, "orbital_radius": 0.091 * au},
            {"name": "c", "mass": 2.9 * M_earth, "orbital_radius": 0.106 * au},
            {"name": "d", "mass": 7.3 * M_earth, "orbital_radius": 0.155 * au},
            {"name": "e", "mass": 8.0 * M_earth, "orbital_radius": 0.195 * au},
            {"name": "f", "mass": 2.0 * M_earth, "orbital_radius": 0.250 * au},
            {"name": "g", "mass": 25.0 * M_earth, "orbital_radius": 0.462 * au},
        ],
    },
    
    "trappist_1": {
        "name": "TRAPPIST-1 system",
        "star_mass": 0.0898 * M_sun,
        "star_radius": 0.1192 * R_sun,
        "description": "7 Earth-sized planets",
        "planets": [
            {"name": "b", "mass": 1.017 * M_earth, "orbital_radius": 0.01154 * au},
            {"name": "c", "mass": 1.156 * M_earth, "orbital_radius": 0.01580 * au},
            {"name": "d", "mass": 0.297 * M_earth, "orbital_radius": 0.02228 * au},
            {"name": "e", "mass": 0.772 * M_earth, "orbital_radius": 0.02925 * au},
            {"name": "f", "mass": 0.934 * M_earth, "orbital_radius": 0.03849 * au},
            {"name": "g", "mass": 1.148 * M_earth, "orbital_radius": 0.04683 * au},
            {"name": "h", "mass": 0.331 * M_earth, "orbital_radius": 0.06189 * au},
        ],
    },
    
    "kepler_90": {
        "name": "Kepler-90 system",
        "star_mass": 1.13 * M_sun,
        "star_radius": 1.2 * R_sun,
        "description": "8-planet system (like Solar System)",
        "planets": [
            {"name": "b", "mass": 2.0 * M_earth, "orbital_radius": 0.074 * au},
            {"name": "c", "mass": 3.0 * M_earth, "orbital_radius": 0.089 * au},
            {"name": "i", "mass": 2.5 * M_earth, "orbital_radius": 0.1234 * au},
            {"name": "d", "mass": 8.0 * M_earth, "orbital_radius": 0.32 * au},
            {"name": "e", "mass": 10.0 * M_earth, "orbital_radius": 0.42 * au},
            {"name": "f", "mass": 12.0 * M_earth, "orbital_radius": 0.48 * au},
            {"name": "g", "mass": 200.0 * M_earth, "orbital_radius": 0.71 * au},
            {"name": "h", "mass": 300.0 * M_earth, "orbital_radius": 1.01 * au},
        ],
    },
    
    "hd_10180": {
        "name": "HD 10180 system",
        "star_mass": 1.06 * M_sun,
        "star_radius": 1.2 * R_sun,
        "description": "7+ planet system",
        "planets": [
            {"name": "b", "mass": 1.35 * M_earth, "orbital_radius": 0.0222 * au},
            {"name": "c", "mass": 13.1 * M_earth, "orbital_radius": 0.0641 * au},
            {"name": "d", "mass": 11.75 * M_earth, "orbital_radius": 0.1286 * au},
            {"name": "e", "mass": 25.0 * M_earth, "orbital_radius": 0.2699 * au},
            {"name": "f", "mass": 23.9 * M_earth, "orbital_radius": 0.4929 * au},
            {"name": "g", "mass": 21.4 * M_earth, "orbital_radius": 1.422 * au},
            {"name": "h", "mass": 64.4 * M_earth, "orbital_radius": 3.4 * au},
        ],
    },
    
    "55_cancri": {
        "name": "55 Cancri system",
        "star_mass": 0.905 * M_sun,
        "star_radius": 0.943 * R_sun,
        "description": "5-planet system with hot Jupiter",
        "planets": [
            {"name": "e", "mass": 8.08 * M_earth, "orbital_radius": 0.01544 * au},
            {"name": "b", "mass": 263.98 * M_earth, "orbital_radius": 0.1148 * au},
            {"name": "c", "mass": 52.4 * M_earth, "orbital_radius": 0.2403 * au},
            {"name": "f", "mass": 45.8 * M_earth, "orbital_radius": 0.781 * au},
            {"name": "d", "mass": 1214.0 * M_earth, "orbital_radius": 5.74 * au},
        ],
    },
    
    "gj_667c": {
        "name": "Gliese 667 C system",
        "star_mass": 0.33 * M_sun,
        "star_radius": 0.42 * R_sun,
        "description": "Multiple habitable zone planets",
        "planets": [
            {"name": "b", "mass": 5.6 * M_earth, "orbital_radius": 0.0505 * au},
            {"name": "c", "mass": 3.8 * M_earth, "orbital_radius": 0.125 * au},
            {"name": "d", "mass": 5.1 * M_earth, "orbital_radius": 0.276 * au},
            {"name": "e", "mass": 2.7 * M_earth, "orbital_radius": 0.213 * au},
            {"name": "f", "mass": 2.7 * M_earth, "orbital_radius": 0.156 * au},
        ],
    },
    
    "hr_8799": {
        "name": "HR 8799 system",
        "star_mass": 1.56 * M_sun,
        "star_radius": 1.44 * R_sun,
        "description": "4 directly imaged giant planets",
        "planets": [
            {"name": "e", "mass": 2857.0 * M_earth, "orbital_radius": 15.0 * au},
            {"name": "d", "mass": 2857.0 * M_earth, "orbital_radius": 27.0 * au},
            {"name": "c", "mass": 2857.0 * M_earth, "orbital_radius": 43.0 * au},
            {"name": "b", "mass": 1905.0 * M_earth, "orbital_radius": 68.0 * au},
        ],
    },
    
    "proxima_centauri_system": {
        "name": "Proxima Centauri system",
        "star_mass": 0.1221 * M_sun,
        "star_radius": 0.1542 * R_sun,
        "description": "Nearest exoplanet system",
        "planets": [
            {"name": "b", "mass": 1.27 * M_earth, "orbital_radius": 0.0485 * au},
            {"name": "c", "mass": 7.0 * M_earth, "orbital_radius": 1.489 * au},
            {"name": "d", "mass": 0.26 * M_earth, "orbital_radius": 0.02885 * au},
        ],
    },
    
    "tau_ceti_system": {
        "name": "Tau Ceti system",
        "star_mass": 0.783 * M_sun,
        "star_radius": 0.793 * R_sun,
        "description": "Nearby Sun-like star with planets",
        "planets": [
            {"name": "g", "mass": 1.75 * M_earth, "orbital_radius": 0.133 * au},
            {"name": "h", "mass": 1.83 * M_earth, "orbital_radius": 0.243 * au},
            {"name": "e", "mass": 3.93 * M_earth, "orbital_radius": 0.538 * au},
            {"name": "f", "mass": 3.93 * M_earth, "orbital_radius": 1.334 * au},
        ],
    },
}


BINARY_SYSTEMS = {
    # ==========================================================================
    # BINARY SYSTEMS (10+ systems)
    # ==========================================================================
    
    "psr_j0737": {
        "name": "PSR J0737-3039 (Double Pulsar)",
        "type": "Binary neutron star",
        "mass_a": 1.3381 * M_sun,
        "mass_b": 1.2489 * M_sun,
        "separation": 8.0e8 * u.m,
        "orbital_period": 2.4 * u.hour,
        "description": "Closest known double pulsar",
    },
    
    "sirius_ab": {
        "name": "Sirius A/B",
        "type": "Binary (A1V + white dwarf)",
        "mass_a": 2.063 * M_sun,
        "mass_b": 1.018 * M_sun,
        "separation": 20.0 * au,
        "orbital_period": 50.1 * u.year,
        "description": "Brightest binary system",
    },
    
    "alpha_centauri_ab": {
        "name": "Alpha Centauri A/B",
        "type": "Binary (G2V + K1V)",
        "mass_a": 1.1 * M_sun,
        "mass_b": 0.907 * M_sun,
        "separation": 23.4 * au,
        "orbital_period": 79.91 * u.year,
        "description": "Nearest stellar binary",
    },
    
    "hulse_taylor": {
        "name": "PSR B1913+16 (Hulse-Taylor)",
        "type": "Binary neutron star",
        "mass_a": 1.4398 * M_sun,
        "mass_b": 1.3886 * M_sun,
        "separation": 1.95e9 * u.m,
        "orbital_period": 7.75 * u.hour,
        "description": "First binary pulsar (Nobel Prize)",
    },
    
    "algol": {
        "name": "Algol (Beta Persei)",
        "type": "Eclipsing binary",
        "mass_a": 3.17 * M_sun,
        "mass_b": 0.70 * M_sun,
        "separation": 0.062 * au,
        "orbital_period": 2.867 * u.day,
        "description": "Prototype eclipsing binary",
    },
    
    "cygnus_x1_binary": {
        "name": "Cygnus X-1 system",
        "type": "X-ray binary (BH + O-star)",
        "mass_a": 21.2 * M_sun,
        "mass_b": 40.6 * M_sun,
        "separation": 0.2 * au,
        "orbital_period": 5.6 * u.day,
        "description": "First confirmed BH binary",
    },
    
    "gw170817": {
        "name": "GW170817 progenitor",
        "type": "Binary neutron star (merged)",
        "mass_a": 1.46 * M_sun,
        "mass_b": 1.27 * M_sun,
        "separation": 300.0 * u.km,
        "orbital_period": 0.01 * u.s,
        "description": "First NS merger with EM counterpart",
    },
    
    "capella": {
        "name": "Capella (Alpha Aurigae)",
        "type": "Binary (G8III + G1III)",
        "mass_a": 2.69 * M_sun,
        "mass_b": 2.56 * M_sun,
        "separation": 0.72 * au,
        "orbital_period": 104.0 * u.day,
        "description": "Bright giant binary",
    },
}


# =============================================================================
# Helper Functions
# =============================================================================

def compute_schwarzschild_radius(M: u.Quantity) -> u.Quantity:
    """
    Berechne Schwarzschild-Radius: r_s = 2GM/c²
    
    Parameters
    ----------
    M : astropy.Quantity
        Masse
        
    Returns
    -------
    r_s : astropy.Quantity
        Schwarzschild-Radius
    """
    r_s = 2 * G * M / c**2
    return r_s.to(u.km)


def compute_photon_sphere(M: u.Quantity) -> u.Quantity:
    """
    Berechne Photon-Sphere Radius: r_ph = 3GM/c² = 1.5 r_s
    
    Parameters
    ----------
    M : astropy.Quantity
        Masse
        
    Returns
    -------
    r_ph : astropy.Quantity
        Photon-Sphere Radius
    """
    r_ph = 3 * G * M / c**2
    return r_ph.to(u.km)


def compute_isco(M: u.Quantity) -> u.Quantity:
    """
    Berechne ISCO (Innermost Stable Circular Orbit): r_ISCO = 6GM/c² = 3 r_s
    
    Parameters
    ----------
    M : astropy.Quantity
        Masse
        
    Returns
    -------
    r_ISCO : astropy.Quantity
        ISCO Radius
    """
    r_ISCO = 6 * G * M / c**2
    return r_ISCO.to(u.km)


def get_system_data(system_name: str, category: str = "stellar") -> dict:
    """
    Hole Systemdaten aus vordefinierten Katalogen.
    
    Parameters
    ----------
    system_name : str
        Name des Systems (z.B. "sun", "sirius_a")
    category : str
        Kategorie: "stellar", "exoplanet", "binary"
        
    Returns
    -------
    data : dict
        Systemdaten mit Masse, Radius, etc.
    """
    if category == "stellar":
        if system_name not in STELLAR_SYSTEMS:
            raise ValueError(f"Unknown stellar system: {system_name}")
        data = STELLAR_SYSTEMS[system_name].copy()
        
    elif category == "exoplanet":
        if system_name not in EXOPLANET_SYSTEMS:
            raise ValueError(f"Unknown exoplanet system: {system_name}")
        data = EXOPLANET_SYSTEMS[system_name].copy()
        
    elif category == "binary":
        if system_name not in BINARY_SYSTEMS:
            raise ValueError(f"Unknown binary system: {system_name}")
        data = BINARY_SYSTEMS[system_name].copy()
        
    else:
        raise ValueError(f"Unknown category: {category}")
    
    # Berechne zusätzliche Parameter für kompakte Objekte
    if "mass" in data and data.get("radius") is None:
        M = data["mass"]
        r_s = compute_schwarzschild_radius(M)
        r_ph = compute_photon_sphere(M)
        r_ISCO = compute_isco(M)
        
        data["schwarzschild_radius"] = r_s
        data["photon_sphere"] = r_ph
        data["ISCO"] = r_ISCO
        
        # Setze r_in auf ISCO wenn nicht definiert
        if data.get("r_in") is None:
            data["r_in"] = r_ISCO
    
    return data


def query_simbad(object_name: str) -> dict:
    """
    Hole Daten von SIMBAD (falls astroquery verfügbar).
    
    Parameters
    ----------
    object_name : str
        Name des Objekts (z.B. "Sirius", "Betelgeuse")
        
    Returns
    -------
    data : dict
        SIMBAD-Daten
    """
    if not HAS_ASTROQUERY:
        print("⚠️  astroquery nicht verfügbar. Verwende vordefinierte Daten.")
        return None
    
    try:
        # Konfiguriere SIMBAD
        customSimbad = Simbad()
        customSimbad.add_votable_fields('sptype', 'ra', 'dec', 'plx', 'pmra', 'pmdec')
        
        # Query
        result_table = customSimbad.query_object(object_name)
        
        if result_table is None:
            print(f"❌ Objekt '{object_name}' nicht in SIMBAD gefunden")
            return None
        
        # Extrahiere Daten
        result = result_table[0]
        
        data = {
            "name": object_name,
            "ra": result['RA'],
            "dec": result['DEC'],
            "spectral_type": result['SP_TYPE'] if 'SP_TYPE' in result.colnames else None,
            "parallax": result['PLX_VALUE'] if 'PLX_VALUE' in result.colnames else None,
            "source": "SIMBAD",
        }
        
        return data
        
    except Exception as e:
        print(f"❌ SIMBAD-Abfrage fehlgeschlagen: {e}")
        return None


def list_available_systems():
    """
    Zeige alle verfügbaren Systeme an.
    """
    print("\n" + "=" * 80)
    print("VERFÜGBARE ASTRONOMISCHE SYSTEME")
    print("=" * 80)
    
    print("\n[STARS] STERNE:")
    print("-" * 80)
    for key, data in STELLAR_SYSTEMS.items():
        print(f"  {key:20s} - {data['name']:30s} ({data['type']})")
        if data['mass']:
            print(f"                         M = {data['mass'].to(M_sun):.3f}")
    
    print("\n[EXOPLANETS] EXOPLANETEN-SYSTEME:")
    print("-" * 80)
    for key, data in EXOPLANET_SYSTEMS.items():
        print(f"  {key:20s} - {data['name']:30s}")
        print(f"                         {len(data['planets'])} Planeten, M_star = {data['star_mass'].to(M_sun):.3f}")
    
    print("\n[BINARY] BINAERSYSTEME:")
    print("-" * 80)
    for key, data in BINARY_SYSTEMS.items():
        print(f"  {key:20s} - {data['name']:30s} ({data['type']})")
        print(f"                         M_A = {data['mass_a'].to(M_sun):.3f}, M_B = {data['mass_b'].to(M_sun):.3f}")
    
    print("\n" + "=" * 80)


def save_system_data(filename: str = "astronomical_systems.json"):
    """
    Speichere alle Systemdaten als JSON.
    
    Parameters
    ----------
    filename : str
        Output-Dateiname
    """
    # Konvertiere zu serialisierbarem Format
    output = {
        "stellar_systems": {},
        "exoplanet_systems": {},
        "binary_systems": {},
    }
    
    for key, data in STELLAR_SYSTEMS.items():
        output["stellar_systems"][key] = {
            k: str(v) if hasattr(v, 'unit') else v
            for k, v in data.items()
        }
    
    for key, data in EXOPLANET_SYSTEMS.items():
        output["exoplanet_systems"][key] = {
            k: str(v) if hasattr(v, 'unit') else v
            for k, v in data.items()
        }
    
    for key, data in BINARY_SYSTEMS.items():
        output["binary_systems"][key] = {
            k: str(v) if hasattr(v, 'unit') else v
            for k, v in data.items()
        }
    
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    
    print(f"💾 Daten gespeichert: {filename}")


# =============================================================================
# Main / Tests
# =============================================================================

if __name__ == "__main__":
    
    print("=" * 80)
    print("FETCH REAL ASTRONOMICAL DATA")
    print("=" * 80)
    
    # Zeige verfügbare Systeme
    list_available_systems()
    
    # Test: Hole Daten für einige Systeme
    print("\n" + "=" * 80)
    print("BEISPIEL: DATEN FÜR EINZELNE SYSTEME")
    print("=" * 80)
    
    test_systems = ["sun", "sirius_a", "sgr_a_star", "neutron_star"]
    
    for sys_name in test_systems:
        print(f"\n{sys_name.upper()}:")
        print("-" * 40)
        data = get_system_data(sys_name, category="stellar")
        
        print(f"  Name:        {data['name']}")
        print(f"  Type:        {data['type']}")
        print(f"  Mass:        {data['mass'].to(M_sun) if data.get('mass') else 'N/A'}")
        
        if data.get('schwarzschild_radius'):
            print(f"  r_s:         {data['schwarzschild_radius']:.3e}")
            print(f"  r_ph:        {data['photon_sphere']:.3e}")
            print(f"  r_ISCO:      {data['ISCO']:.3e}")
        
        print(f"  r_in:        {data['r_in']}")
        print(f"  r_out:       {data['r_out']}")
    
    # Speichere alle Daten
    print("\n" + "=" * 80)
    save_system_data("astronomical_systems.json")
    
    # Optional: SIMBAD-Test
    if HAS_ASTROQUERY:
        print("\n" + "=" * 80)
        print("TEST: SIMBAD-Abfrage")
        print("=" * 80)
        
        simbad_data = query_simbad("Sirius")
        if simbad_data:
            print("\nSIMBAD-Daten für Sirius:")
            for key, val in simbad_data.items():
                print(f"  {key:20s}: {val}")
    
    print("\n" + "=" * 80)
    print("✅ Fetch-Script erfolgreich ausgeführt")
    print("=" * 80)
