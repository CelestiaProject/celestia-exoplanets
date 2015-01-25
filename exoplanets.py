#!/usr/bin/python3
#
# celestia-exoplanets: Create an exoplanet catalogue for Celestia
# Copyright (C) 2014  Andrew Tribick
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301, USA.

# For description and acknowledgements, see the file headers contained in the
# constants STC_HEADER and SSC_HEADER below.

import csv
import os.path
import time
import re
from numpy import arctan2, cos, degrees, log10, pi, radians, sin, sqrt
from celestiautils import StarDatabase
from orbitconvert import OrbitElements, OrbitConverter
from optparse import OptionParser


STC_HEADER = """# Catalog of exoplanet host stars.
# Generated at {0}
#
# Data extracted from the NASA Exoplanet Archive
#     http://exoplanetarchive.ipac.caltech.edu/
#
# Stars which are not present in Celestia's default data files are included if
# there is sufficient information present in the NASA Exoplanet Archive data
# to model them in Celestia, and at least one of the planets meets the output
# requirements (see the extrasolar.ssc file for details).
#
# As the NASA Exoplanet Archive does not contain information on binary and
# multiple star systems, these are not currently included.
#
# Distance information is taken from the following sources, in decreasing
# order of preference: the distance field 'st_dist' in the CSV file, the
# parallax 'st_plx' if greater than the parallax error 'st_plxerr', the
# combination of V magnitude and bolometric magnitude/correction estimated
# from the stellar effective temperature 'st_teff' and radius 'st_rad', and
# the combination of 2MASS J, H, Ks magnitudes and bolometric magnitude and
# correction estimated from the stellar properties. The formulae for
# converting the 2MASS magnitudes are taken from Bilir et al. (2008). The
# formula for estimating bolometric correction from stellar effective
# temperature is from Reed (1998).
#
# References:
#
# Bilir et al. (2008) "Transformations between 2Mass, SDSS and BVRI
#     photometric systems: bridging the near-infrared and optical",
#     MNRAS v384, p1178-1188
#
# Reed, B. C. (1998) "The Composite Observational-Theoretical HR Diagram",
#     Journal of the Royal Astronomical Society of Canada, v92, p36
#
# Acknowledgements:
#
# Thanks goes to Grant Hutchison who previously maintained the extrasolar
# planet and host star catalogues. Without his efforts Celestia would be a far
# less rich environment. His orbit transformation spreadsheets were invaluable
# in verifying the code to process orbital elements, and several decisions in
# the output of the script were inspired by his extrasolar planets files.


"""

SSC_HEADER = """# Catalog of confirmed exoplanets.
# Generated at {0}
#
# Data extracted from the NASA Exoplanet Archive
#     http://exoplanetarchive.ipac.caltech.edu/
#
# Planets are included if the host star is already present in Celestia's
# catalogue files, or if there is sufficient information in the NASA Exoplanet
# archive to create the star (see the extrasolar.stc file for details). The 
# planets are then further filtered to ensure that the orbital period,
# semimajor axis and radius can be obtained. If one of the orbital period and
# semimajor axis is missing and the stellar mass value is available, the
# missing value is calculated using Kepler's third law. If the radius is
# missing, it will be estimated from the mass or minimum mass if either is
# present. The solar system mass-radius relationship M = R^2.06 (from Lissauer
# et al. 2011) is used.
#
# As the NASA Exoplanet Archive does not provide values for the position angle
# of the orbital node, none of the orbits can be considered to have fully-
# specified orientations. In case of missing node and inclination information,
# the values from the ecliptic plane (Celestia's default reference frame) are
# used to fill in the gaps. Some planets are noted in the NASA Exoplanet
# Archive as being transiting but lack inclination information, in this case
# an inclination of 90 degrees is assumed.
#
# References:
#
# Lissauer et al. (2011) "Architecture and Dynamics of Kepler's Candidate
#     Multiple Transiting Planet Systems", ApJ Supplement v197, article id. 8
#
# Acknowledgements:
#
# Thanks goes to Grant Hutchison who previously maintained the extrasolar
# planet and host star catalogues. Without his efforts Celestia would be a far
# less rich environment. His orbit transformation spreadsheets were invaluable
# in verifying the code to process orbital elements, and several decisions in
# the output of the script were inspired by his extrasolar planets files.


"""


parser = OptionParser()
parser.add_option("-c", "--celestia-dir", dest="celestia_dir",
                  help="path to Celestia directory", metavar="DIR")

(options, args) = parser.parse_args()

if (options.celestia_dir is None):
    parser.error('No Celestia directory was specified.')
else:
    celestia_path = os.path.expanduser(options.celestia_dir)


STC_FILES = [
    'revised.stc',
    'nearstars.stc',
    'visualbins.stc',
    'spectbins.stc',
    ]

stardb = StarDatabase()

data_dir = os.path.join(celestia_path, 'data')

stars_dat = os.path.join(data_dir, 'stars.dat')
starnames_dat = os.path.join(data_dir, 'starnames.dat')
hdxindex_dat = os.path.join(data_dir, 'hdxindex.dat')

stardb.parse_files((stars_dat, starnames_dat, hdxindex_dat))
stardb.parse_files((os.path.join(data_dir, fname) for fname in STC_FILES))

systems = {}

greeks = (
    'alf|bet|gam|del|eps|zet|eta|tet|iot|kap|lam|mu|nu|xi|omi|pi|rho|sig|tau|'
    'ups|phi|chi|psi|ome'
    )

constellations = (
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
            return 'Gliese ' + system_key[3:]
    
    m = re.match('^(' + greeks + ')(?: *([0-9]+)) +(' + constellations + ')$',
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


def create_star(row):
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
        hip = stardb.lookup_name(name)
        if hip is not None:
            preferred_name = name
            exists = (hip == -1 or hip in stardb.hip_stars)
            break
    
    host_names[:] = [name for name in host_names
                     if not stardb.is_catalog_name(name)]
    
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
    elif (row['st_plx'] and row['st_plxerr'] and
          float(row['st_plx']) > float(row['st_plxerr'])):
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


def mass_to_radius(mass):
    """Planetary mass-radius relationship."""
    radius = mass ** (1 / 2.06)
    if radius > 15:
        radius = 15
    return round(radius, 1)


def kepler3_semimajor(period, stmass):
    """Compute the semimajor axis using Kepler's 3rd law."""
    return round(((2.95995e-4*stmass*(float(row['pl_orbper'])**2)) /
                  (4*(pi**2))) ** (1/3), 3)


def kepler3_period(semimajor, stmass):
    """Compute the period using Kepler's 3rd law."""
    return round(sqrt(4 * (pi**2) * (semimajor**3) / (2.95995e-4*stmass)), 3)


def transit_anomaly(arg_peri, ecc):
    """Calculate the mean anomaly at transit time."""
    theta = radians(90-arg_peri)
    E = arctan2(sqrt(1-ecc**2) * sin(theta), ecc+cos(theta))
    M = degrees(E - ecc*sin(E))
    if M < 0:
        return M + 360
    else:
        return M


def create_planet(row, star_data):
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


def star_ok(system):
    """Checks if required information is present for the star."""
    star_data = system['star']
    if not system['planets']:
        print(star_data['preferred_name'], 'has no representable planets')
        return False
    if star_data['exists']:
        return True
    if 'ra' not in star_data:
        print(star_data['preferred_name'], 'missing right ascension')
        return False
    if 'dec' not in star_data:
        print(star_data['preferred_name'], 'missing declination')
        return False
    if 'dist' not in star_data:
        print(star_data['preferred_name'], 'missing distance')
        return False
    if 'magV' not in star_data:
        print(star_data['preferred_name'], 'missing apparent magnitude')
        return False
    return True


def planet_ok(planet_data, star_name):
    """Checks if required information is present for the planet."""
    if 'a' not in planet_data:
        print(star_name, planet_data['name'], 'missing semimajor axis')
        return False
    if 'P' not in planet_data:
        print(star_name, planet_data['name'], 'missing orbital period')
        return False
    if 'rad' not in planet_data:
        print(star_name, planet_data['name'], 'missing radius')
        return False
    return True


def output_star(f, star):
    """Outputs star data to an stc file."""
    if star['hip'] is not None:
        f.write(str(star['hip']))
    if star['names']:
        if star['hip'] is not None:
            f.write(' ')
        f.write('"' + ':'.join(star['names']) + '"\n')
    else:
        f.write('\n')
    f.write('{\n')
    f.write('\tRA<deg> {0:.6f}\n'.format(star['ra']))
    f.write('\tDec<deg> {0:.6f}\n'.format(star['dec']))
    
    if '#dist' in star:
        f.write('\tDistance<pc> {0:g} # {1}\n'.format(star['dist'],
                                                      star['#dist']))
    else:
        f.write('\tDistance<pc> {0:g}\n'.format(star['dist']))

    if '#magV' in star:
        f.write('\tAppMag {0:g} # {1}\n'.format(star['magV'],
                                                star['#magV']))
    else:
        f.write('\tAppMag {0:g}\n'.format(star['magV']))

    if 'spec' in star:
        f.write('\tSpectralType "{0}"\n'.format(star['spec']))
    else:
        f.write('\tSpectralType "?"\n')
    
    if 'teff' in star:
        f.write('\tTemperature {0:g}\n'.format(star['teff']))
    
    if 'bc' in star:
        f.write('\tBoloCorrection {0:.2f}\n'.format(star['bc']))

    if 'rad' in star:
        f.write('\tRadius<rS> {0:g}\n'.format(star['rad']))

    f.write('}\n\n')


def output_planet(f, planet, star):
    """Outputs planet data to an ssc file."""
    f.write('"{0}" "{1}"\n'.format(planet['name'], star['preferred_name']))
    f.write('{\n')
    
    if '#rad' in planet:
        f.write('\tRadius<rE> {0:g} # {1}\n'.format(planet['rad'],
                                                    planet['#rad']))
    else:
        f.write('\tRadius<rE> {0:g}\n'.format(planet['rad']))
    
    if 'mass' in planet:
        if '#mass' in planet:
            f.write('\tMass {0:g} # {1}\n'.format(planet['mass'],
                                                  planet['#mass']))
        else:
            f.write('\tMass {0:g}\n'.format(planet['mass']))
    
    f.write('\n\tEllipticalOrbit {\n')
    
    if '#a' in planet:
        f.write('\t\tSemiMajorAxis<AU> {0:g} # {1}\n'.format(planet['a'],
                                                             planet['#a']))
    else:
        f.write('\t\tSemiMajorAxis<AU> {0:g}\n'.format(planet['a']))
    
    if '#P' in planet:
        f.write('\t\tPeriod<d> {0:g} # {1}\n'.format(planet['P'],
                                                     planet['#P']))
    else:
        f.write('\t\tPeriod<d> {0:g}\n'.format(planet['P']))
    
    if 'e' in planet:
        f.write('\t\tEccentricity {0:g}\n'.format(planet['e']))
    
    if 't_transit' in planet:
        f.write('\n\t\tEpoch {0:.6f} # transit midpoint\n'.format(
            planet['t_transit'],
            ))
        arg_peri = planet['arg_peri'] if 'arg_peri' in planet else 0
        ecc = planet['e'] if 'e' in planet else 0
        f.write('\t\tMeanAnomaly {0:.3f}\n'.format(
            transit_anomaly(arg_peri, ecc)))
        pass
    elif 't_peri' in planet:
        f.write('\n\t\tEpoch {0:.6f} # periastron\n'.format(planet['t_peri']))
        f.write('\t\tMeanAnomaly<deg> 0\n')
        pass
    
    
    converter = OrbitConverter(star['ra'], star['dec'])
    elements = converter.ecliptic_plane()
    # orbital elements quoted are for the stellar reflex orbit, so rotate by
    # 180 degrees to get the planet's orbit.
    if 'arg_peri' in planet:
        elements.arg_peri = planet['arg_peri'] + 180
    else:
        elements.arg_peri = 180
        planet['#arg_peri'] = 'unknown, 0 assumed'
    
    if 'inc' in planet:
        elements.inclination = planet['inc']
    else:
        planet['#inc'] = 'unknown, using ecliptic'
        
    planet['#node'] = 'unknown, using ecliptic'
    
    elements = converter.transform(elements)
    if '#arg_peri' in planet:
        f.write('\n\t\tArgOfPericenter<deg> {0:0.3f} # {1}\n'.format(
            elements.arg_peri,
            planet['#arg_peri'],
            ))
    else:
        f.write('\n\t\tArgOfPericenter<deg> {0:0.3f}\n'.format(
            elements.arg_peri,
            ))
    
    if '#inc' in planet:
        f.write('\t\tInclination<deg> {0:0.3f} # {1}\n'.format(
            elements.inclination,
            planet['#inc'],
            ))
    else:
        f.write('\t\tInclination<deg> {0:0.3f}\n'.format(
            elements.inclination,
            ))
    
    if '#node' in planet:
        f.write('\t\tAscendingNode<deg> {0:0.3f} # {1}\n'.format(
            elements.node,
            planet['#node'],
            ))
    else:
        f.write('\t\tAscendingNode<deg> {0:0.3f} # {1}\n'.format(
            elements.node,
            ))
    
    f.write('\t}\n\n\tUniformRotation {\n')
    f.write('\t\tInclination<deg> {0:0.3f} # to match orbit\n'.format(
        elements.inclination,
        ))
    f.write('\t\tAscendingNode<deg> {0:0.3f} # to match orbit\n'.format(
        elements.node,
        ))
    
    f.write('\t}\n}\n\n')


with open('exoplanets.csv', 'r') as f:
    csv = csv.DictReader(f)
    for row in csv:
        system_key = row["pl_hostname"]
        
        if system_key not in systems:
            systems[system_key] = {
                'star': create_star(row),
                'planets': [],
                }
        
        system = systems[system_key]
        system['planets'].append(create_planet(row, system['star']))

for system_key, system in systems.items():
    preferred_name = system['star']['preferred_name']
    system['planets'][:] = [planet for planet in system['planets']
                            if planet_ok(planet, preferred_name)]
    system['planets'].sort(key=lambda planet: planet['a'])

systems = {system_key: system
           for system_key, system in systems.items()
           if star_ok(system)}

with open('extrasolar.stc', 'w') as out_stars, \
        open('extrasolar.ssc', 'w') as out_planets:

    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S%z', time.localtime())
    out_stars.write(STC_HEADER.format(timestamp))
    out_planets.write(SSC_HEADER.format(timestamp))
    
    for system in sorted(systems.values(), key=lambda s: s['star']['ra']):
        if not system['star']['exists']:
            output_star(out_stars, system['star'])
        for planet in system['planets']:
            output_planet(out_planets, planet, system['star'])