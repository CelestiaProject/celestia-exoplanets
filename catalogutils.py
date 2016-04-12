#!/usr/bin/python3


from __future__ import division

import csv
import re
from numpy import log10, pi, sqrt


GREEKS = (
    'alf|bet|gam|del|eps|zet|eta|tet|iot|kap|lam|mu|nu|xi|omi|pi|rho|sig|tau|'
    'ups|phi|chi|psi|ome'
    )

CONSTELLATIONS = (
    'And|Ant|Aps|Aql|Aqr|Ara|Ari|Aur|Boo|Cae|Cam|Cap|Car|Cas|Cen|Cep|Cet|Cha|'
    'Cir|CMa|CMi|Cnc|Col|Com|CrA|CrB|Crt|Cru|Crv|CVn|Cyg|Del|Dor|Dra|Equ|Eri|'
    'For|Gem|Gru|Her|Hor|Hya|Hyi|Ind|Lac|Leo|Lep|Lib|LMi|Lup|Lyn|Lyr|Men|Mic|'
    'Mon|Mus|Nor|Oct|Oph|Ori|Pav|Peg|Per|Phe|Pic|PsA|Psc|Pup|Pyx|Ret|Scl|Sco|'
    'Sct|Ser|Sex|Sge|Sgr|Tau|Tel|TrA|Tri|Tuc|UMa|UMi|Vel|Vir|Vol|Vul'
    )


def normalise_name(name):
    """Normalises names to follow Celestia conventions."""
    if name[0:3] == 'GJ ':
        split_gliese = name.split(' ')
        if int(split_gliese[1]) < 1000:
            return 'Gliese ' + name[3:]
    
    m = re.match('^(' + GREEKS + ')(?: *([0-9]+)) +(' + CONSTELLATIONS + ')$',
                 name)
    if m:
        if m.groups(2) is None:
            return m.group(1).upper() + ' ' + m.group(3)
        else:
            return m.group(1).upper() + m.group(2) + ' ' + m.group(3)
    
    return name


def remove_primary_suffix(catalog_name):
    """Removes the A suffix for catalog names."""
    if catalog_name[-2:] == ' A':
        return catalog_name[:-2]
    return catalog_name


def filter_names(names):
    """Removes empty and duplicate names."""
    seen = set()
    result = []
    for name in names:
        if not name or name in seen:
            continue
        seen.add(name)
        result.append(name)
    return result


def teff_to_bc(teff):
    """Convert effective temperature to bolometric correction."""
    logT4 = log10(teff)-4
    return (-8.499 * (logT4**4) +
            13.421 * (logT4**3) -
            8.131 * (logT4**2) -
            3.901 * logT4 - 0.438)


def jhk_to_v(j, h, k):
    """Convert infrared magnitudes to visual."""
    jh = j - h
    hk = h - k
    bv = 1.622*jh + 0.912*hk + 0.044
    ri = 0.954*jh + 0.593*hk + 0.025
    vk = 1.896*bv + 1.131*ri - 0.004
    return round(k + vk, 1)


def mag_to_distance(v, radius, teff):
    """Calculate photometric distance."""
    bolo_mag = v + teff_to_bc(teff)
    bolo_abs = 4.75 - 5*log10(radius) - 10*log10(teff/5780.0)
    return round(10 ** ((bolo_mag-bolo_abs)/5 + 1), 0)


def mass_to_radius(mass):
    """Planetary mass-radius relationship."""
    radius = mass ** (1 / 2.06)
    if radius > 15:
        radius = 15
    return round(radius, 1)


def kepler3_semimajor(period, stmass):
    """Compute the semimajor axis using Kepler's 3rd law."""
    return round(((2.95995e-4*stmass*(period**2)) / (4*(pi**2))) ** (1/3), 3)


def kepler3_period(semimajor, stmass):
    """Compute the period using Kepler's 3rd law."""
    return round(sqrt(4 * (pi**2) * (semimajor**3) / (2.95995e-4*stmass)), 3)


class NasaCatalog(object):
    """Class to handle the NASA Exoplanet Database"""
    
    def __init__(self, stardb):
        self.stardb = stardb


    def parse_file(self, f):
        rows = csv.DictReader(f)
        systems = {}
        
        for row in rows:
            system_key = row["pl_hostname"]
            
            if system_key not in systems:
                systems[system_key] = {
                    'star': self.create_star(row),
                    'planets': [],
                    }
            
            system = systems[system_key]
            system['planets'].append(self.create_planet(row, system['star']))
            
        return systems


    def create_star(self, row):
        """Create a star data record."""
        host_names = filter_names([
            normalise_name(row["pl_hostname"]),
            row['hd_name'],
            row['hip_name'],
            remove_primary_suffix(row['hd_name']),
            remove_primary_suffix(row['hip_name']),
            ])
        
        hip = None
        exists = False
        preferred_name = host_names[0]
        
        for name in host_names:
            hip = self.stardb.lookup_name(name)
            if hip is not None:
                preferred_name = name
                exists = (hip == -1 or hip in self.stardb.hip_stars)
                break
        
        host_names[:] = [name for name in host_names
                        if not self.stardb.is_catalog_name(name)]
        
        star_data = {
            'hip': hip,
            'exists': exists,
            'names': host_names,
            'preferred_name': preferred_name,
            }
        
        if row['ra'] and row['dec']:
            star_data['ra'] = float(row['ra'])
            star_data['dec'] = float(row['dec'])
        
        if row['st_vj']:
            star_data['magV'] = float(row['st_vj'])
        elif row['st_j'] and row['st_h'] and row['st_k']:
            star_data['magV'] = jhk_to_v(float(row['st_j']),
                                        float(row['st_h']),
                                        float(row['st_k']),
                )
            star_data['#magV'] = 'estimated from 2MASS magnitudes'

        if row['st_rad']:
            star_data['rad'] = float(row['st_rad'])

        if row['st_mass']:
            star_data['mass'] = float(row['st_mass'])

        if row['st_teff']:
            star_data['teff'] = float(row['st_teff'])
            star_data['bc'] = teff_to_bc(float(row['st_teff']))

        if row['st_dist']:
            star_data['dist'] = float(row['st_dist'])
        elif row['st_plx'] and row['st_plxerr1']:
            plx = float(row['st_plx'])
            plxerr = float(row['st_plxerr2']) if row['st_plxerr2'] else float(row['st_plxerr1'])
            if plx > abs(plxerr):
                star_data['dist'] = round(1000/float(row['st_plx']), 0)
                star_data['#dist'] = 'from parallax value'
        elif 'magV' in star_data and 'rad' in star_data and 'teff' in star_data:
            star_data['dist'] = mag_to_distance(star_data['magV'],
                                                star_data['rad'],
                                                star_data['teff'])
            star_data['#dist'] = 'from photometry'

        if row['st_spstr']:
            if row['st_spstr'] == 'WD':
                star_data['spec'] = 'D'
            else:
                star_data['spec'] = row['st_spstr']

        return star_data


    def create_planet(self, row, star_data):
        """Create a planet data record."""
        planet_data = {'name': row['pl_letter']}
        
        if row['pl_masse']:
            planet_data['mass'] = float(row['pl_masse'])
        elif row['pl_msinie']:
            planet_data['mass'] = float(row['pl_msinie'])
            planet_data['#mass'] = 'minimum mass'
        
        if row['pl_rade']:
            planet_data['rad'] = float(row['pl_rade'])
        elif 'mass' in planet_data:
            planet_data['rad'] = mass_to_radius(planet_data['mass'])
            planet_data['#rad'] = 'from mass-radius relationship'

        if row['pl_orbsmax']:
            planet_data['a'] = float(row['pl_orbsmax'])
        elif row['pl_orbper'] and 'mass' in star_data:
            planet_data['a'] = kepler3_semimajor(float(row['pl_orbper']),
                                                star_data['mass'])
            planet_data['#a'] = "from Kepler's 3rd law"

        if row['pl_orbper']:
            planet_data['P'] = float(row['pl_orbper'])
        elif 'a' in planet_data and 'mass' in star_data:
            planet_data['P'] = kepler3_period(planet_data['a'], star_data['mass'])
            planet_data['#P'] = "from Kepler's 3rd law"

        if row['pl_orbeccen']:
            planet_data['e'] = float(row['pl_orbeccen'])

        if row['pl_orblper']:
            planet_data['arg_peri'] = float(row['pl_orblper'])
        
        if row['pl_tranmid']:
            planet_data['t_transit'] = float(row['pl_tranmid'])
        
        if row['pl_orbtper']:
            planet_data['t_peri'] = float(row['pl_orbtper'])

        if row['pl_orbincl']:
            planet_data['inc'] = float(row['pl_orbincl'])
        elif row['pl_tranflag'] == '1':
            planet_data['inc'] = 90
            planet_data['#inc'] = 'Assumed 90 degrees for transiting planet'

        return planet_data
