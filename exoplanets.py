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


from __future__ import print_function

import os.path
import time
from numpy import arctan2, cos, degrees, radians, sin, sqrt
from celestiautils import StarDatabase
from catalogutils import NasaCatalog
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


def transit_anomaly(arg_peri, ecc):
    """Calculate the mean anomaly at transit time."""
    theta = radians(90-arg_peri)
    E = arctan2(sqrt(1-ecc**2) * sin(theta), ecc+cos(theta))
    M = degrees(E - ecc*sin(E))
    if M < 0:
        return M + 360
    else:
        return M


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
        f.write('\t\tMeanAnomaly<deg> {0:.3f}\n'.format(
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



def clean_systems(systems):
    """Removes stars and planets which cannot be represented in Celestia."""
    for system_key, system in systems.items():
        preferred_name = system['star']['preferred_name']
        system['planets'][:] = [planet for planet in system['planets']
                                if planet_ok(planet, preferred_name)]
        system['planets'].sort(key=lambda planet: planet['a'])

    return {system_key: system
            for system_key, system in systems.items()
            if star_ok(system)}


def write_systems(systems,out_stars, out_planets):
    """Outputs systems to files."""
    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S%z', time.localtime())
    out_stars.write(STC_HEADER.format(timestamp))
    out_planets.write(SSC_HEADER.format(timestamp))
    
    for system in sorted(systems.values(), key=lambda s: s['star']['ra']):
        if not system['star']['exists']:
            output_star(out_stars, system['star'])
        for planet in system['planets']:
            output_planet(out_planets, planet, system['star'])   


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

catalog = NasaCatalog(stardb)

with open('exoplanets.csv', 'r') as f:
    systems = clean_systems(catalog.parse_file(f))

with open('extrasolar.stc', 'w') as out_stars, \
        open('extrasolar.ssc', 'w') as out_planets:
    write_systems(systems, out_stars, out_planets)
