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
import struct
import time
from numpy import arctan2, cos, degrees, log10, pi, radians, sin, sqrt
from stcparser import StcParser
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
    'spectbins.stc'
]


def parse_stars_dat(f):
    """Extracts the HIP numbers from the stars.dat file."""
    f.seek(14)
    data = f.read()
    hip_stars = {struct.unpack('<I', data[pos:pos+4])[0]
                 for pos in range(0, len(data), 20)}
    return hip_stars


def parse_starnames(f):
    """Extracts the star names from starnames.dat"""
    star_names = dict()
    for line in f:
        split_line = line.split(':')
        hip = int(split_line[0])
        star_names.update({name: hip for name in split_line[1:]})
    return star_names


def parse_hdxindex(f):
    """Extracts the HD names from hdxindex.dat"""
    f.seek(10)
    data = f.read()
    star_names = dict()
    for pos in range(0, len(data), 8):
        hd, hip = struct.unpack("<II", data[pos:pos+8])
        star_names["HD "+str(hd)] = hip
    return star_names


def parse_stc_file(f, parser=None):
    """Extracts the HIP numbers from an stc file."""
    if parser is None:
        parser = StcParser()
    stc_stars = parser.parse(f.read())
    hip_stars = set()
    star_names = dict()
    for star in stc_stars:
        if star["HIP"] is not None:
            hip_stars.add(star["HIP"])
        if star["Name"] is not None:
            split_names = star["Name"].split(':')
            star_names.update({name: -1 for name in split_names})
    return hip_stars, star_names


def star_ok(row):
    """Checks if all required information to create a star is present."""
    is_ok = True
    if row["ra"] == '' or row["dec"] == '':
        is_ok = False
    
    has_jhk = row["st_j"] != '' and row["st_h"] != '' and row["st_k"] != ''
    
    mag_estimatable = row["st_vj"] != '' or has_jhk    
    is_ok &= mag_estimatable
    
    dist_estimatable = True
    if row["st_dist"] == '':
        dist_estimatable = False
        if row["st_plx"] != '' and row["st_plxerr"] != '':
            if float(row["st_plx"]) > float(row["st_plxerr"]):
                dist_estimatable = True
        if row["st_teff"] != '' and row["st_rad"] != '':
            if row["st_vj"] != '' or has_jhk:
                dist_estimatable = True
    is_ok &= dist_estimatable
    
    return is_ok


def planet_ok(row):
    """Checks if all required information to create a planet is present."""
    is_ok = True
    if row["pl_orbper"] == '':
        if row["pl_orbsmax"] == '' or row["st_mass"] == '':
            is_ok = False
    elif row["pl_orbsmax"] == '':
        if row["pl_orbper"] == '' or row["st_mass"] == '':
            is_ok = False
    
    if (row['pl_rade'] == '' and
            row['pl_masse'] == '' and
            row['pl_msinie'] == ''):
        is_ok = False
    
    return is_ok


def out_star(f, row, hip, hostname, star_names):
    """Outputs a star definition to the stc file."""
    if hip is None:
        names = list()
        if (row["hip_name"] != '' and row["hip_name"] != hostname):
            names.append(row["hip_name"])
        names.append(hostname)
        if (row["hd_name"] != '' and row["hd_name"] != hostname):
            names.append(row["hd_name"])
        f.write('"' + ':'.join(names) + '"\n')
    else:
        names = list()
        if hostname not in star_names:
            names.append(hostname)
        if (row["hd_name"] != '' and
                row["hd_name"] != hostname and
                row["hd_name"] not in star_names):
            names.append(hostname)
        f.write(str(hip))
        if len(names) > 0:
            f.write(' "' + ':'.join(names) + '"')
        f.write('\n')

    f.write('{\n')
    f.write('\tRA<deg> ' + row['ra'] + '\n')
    f.write('\tDec<deg> ' + row['dec'] + '\n')
    bolo_correct = None
    t_eff = None
    if row['st_teff'] != '':
        t_eff = float(row['st_teff'])
        logT4 = log10(t_eff) - 4
        bolo_correct = (-8.499 * (logT4**4) +
                        13.421 * (logT4**3) -
                        8.131 * (logT4**2) -
                        3.901 * logT4 - 0.438)
    
    radius = None
    if row['st_rad'] != '':
        radius = float(row['st_rad'])
    
    v_mag = None
    v_actual = True
    if row['st_vj'] != '':
        v_mag = float(row['st_vj'])
    else:
        j = float(row['st_j'])
        h = float(row['st_h'])
        k = float(row['st_k'])
        jh = j - h
        hk = h - k
        bv = 1.622*jh + 0.912*hk + 0.044
        ri = 0.954*jh + 0.593*hk + 0.025
        vk = 1.896*bv + 1.131*ri - 0.004
        v_mag = k + vk
        v_actual = False
    
    if row['st_dist'] != '':
        f.write('\tDistance<pc> ' + row['st_dist'] + '\n')
    elif (row['st_plx'] != '' and
            row['st_plxerr'] != '' and
            float(row['st_plx']) > float(row['st_plxerr'])):
        plx_distance = 1000 / float(row['st_plx'])
        f.write('\tDistance<pc> {0:0.1f} # from parallax value\n'.format(
            plx_distance))
    else:
        bolo_mag = v_mag + bolo_correct
        bolo_abs = 4.75 - 5*log10(radius) - 10*log10(t_eff/5780)
        phot_distance = round(10 ** ((bolo_mag-bolo_abs)/5 + 1), 0)
        if v_actual:
            estimate_message = 'estimated from V magnitude'
        else:
            estimate_message = 'estimated from 2MASS magnitudes'
            
        f.write('\tDistance<pc> {0:0.0f} # {1}\n'.format(phot_distance,
                                                         estimate_message))
        
    f.write('\tAppMag {0:0.2f}\n'.format(v_mag))
    
    if row['st_spstr'] == '':
        f.write('\tSpectralType "?"\n')
    else:
        if row['st_spstr'] == 'WD':
            f.write('\tSpectralType "D"\n')
        else:
            f.write('\tSpectralType "' + row['st_spstr'] + '"\n')
    if row['st_teff'] != '':
        f.write('\tTemperature ' + row['st_teff'] + '\n')
    if bolo_correct is not None:
        f.write('\tBoloCorrection {0:0.2f}\n'.format(bolo_correct))

    if row['st_rad'] != '':
        f.write('\tRadius<rS> ' + row['st_rad'] + '\n')
    
    f.write('}\n\n')

def out_planet(f, row, use_host):
    """Outputs a planet definition to the ssc file."""
    f.write('"{0}" "{1}"\n'.format(row["pl_letter"], use_host))
    f.write('{\n')
    
    if row['pl_rade'] == '':
        if row['pl_masse'] == '':
            mass = float(row['pl_msinie'])
        else:
            mass = float(row['pl_masse'])
        radius = mass ** (1 / 2.06)
        if radius > 15:
            radius = 15
        f.write(
            '\tRadius<rE> {0:0.1f} # from mass-radius relationship\n'.format(
                radius))
    else:
        f.write('\tRadius<rE> ' + row['pl_rade'] + '\n')
    
    if row['pl_masse'] != '':
        f.write('\tMass ' + row['pl_masse'] + '\n')
    elif row['pl_msinie'] != '':
        f.write('\tMass ' + row['pl_msinie'] + ' # minimum mass\n')
    
    f.write('\n\tEllipticalOrbit {\n')
    stmass = None
    if row['st_mass'] != '':
        stmass = float(row['st_mass'])
    
    if row['pl_orbsmax'] == '':
        semimajor = (
            (2.95995e-4 * stmass * (float(row['pl_orbper'])**2)) /
            (4 * (pi**2))) ** (1/3)
        f.write(
            "\t\tSemiMajorAxis<AU> {0:0.2f} # from Kepler's 3rd law\n".format(
                semimajor))
    else:
        semimajor = float(row['pl_orbsmax'])
        f.write('\t\tSemiMajorAxis<AU> ' + row['pl_orbsmax'] + '\n')
    
    if row['pl_orbper'] == '':
        period = sqrt(4 * (pi**2) * (semimajor**3) / (2.95995e-4*stmass))
        f.write("\t\tPeriod<d> {0:0.3f} # from Kepler's 3rd law\n".format(period))
    else:
        f.write('\t\tPeriod<d> ' + row['pl_orbper'] + '\n')
    
    if row['pl_orbeccen'] != '':
        ecc = float(row['pl_orbeccen'])
        f.write('\t\tEccentricity ' + row['pl_orbeccen'] + '\n')
    else:
        ecc = 0
    
    if row['pl_orblper'] == '':
        arg_peri = 0
    else:
        arg_peri = float(row['pl_orblper'])

    if row['pl_tranmid'] != '':
        f.write('\n\t\tEpoch ' + row['pl_tranmid'] + ' # transit midpoint\n')
        
        theta = radians(90 - arg_peri)
        E = arctan2(sqrt(1-ecc**2) * sin(theta), ecc + cos(theta))
        M = degrees(E - ecc*sin(E))
        if M < 0:
            M += 360
        f.write('\t\tMeanAnomaly<deg> {0:0.3f}\n'.format(M))
    elif row['pl_orbtper'] != '':
        f.write(
            '\n\t\tEpoch ' + row['pl_orbtper'] + ' # time of periastron\n')
        f.write('\t\tMeanAnomaly<deg> 0\n')
    
    converter = OrbitConverter(float(row['ra']), float(row['dec']))
    elements = converter.ecliptic_plane()
    # orbital elements quoted are for the stellar reflex orbit, so rotate by
    # 180 degrees to get the planet's orbit.
    elements.arg_peri = arg_peri + 180
    
    if row['pl_orbincl'] == '' and row['pl_tranflag'] != '1':
        f.write('\n\t\t# Unknown orientation: using ecliptic plane\n')
        elements = converter.transform(elements)
        f.write('\t\tArgOfPericenter<deg> {0:0.3f}\n'.format(elements.arg_peri))
        f.write('\t\tInclination<deg> 0\n')
        f.write('\t\tAscendingNode<deg> 0\n')
        f.write('\t}\n\n')
        f.write('\tUniformRotation {\n')
        f.write('\t\tInclination<deg> 0 # to match orbit\n')
        f.write('\t\tAscendingNode<deg> 0 # to match orbit\n')
        f.write('\t}\n')
    else:
        if row['pl_orbincl'] == '':
            f.write(
                '\n\t\t# Transiting planet without inclination, '
                '90 degrees assumed.\n')
            elements.inclination = 90
        else:
            f.write(
                '\n\t\t# Known inclination, using node of ecliptic plane.\n')
            elements.inclination = float(row['pl_orbincl'])
        elements = converter.transform(elements)
        
        f.write('\t\tArgOfPericenter<deg> {0:0.3f}\n'.format(
            elements.arg_peri))
        f.write('\t\tInclination<deg> {0:0.3f}\n'.format(
            elements.inclination))
        f.write('\t\tAscendingNode<deg> {0:0.3f}\n'.format(
            elements.node))
        f.write('\t}\n\n')
        f.write('\tUniformRotation {\n')
        f.write('\t\tInclination<deg> {0:0.3f} # to match orbit\n'.format(
            elements.inclination))
        f.write('\t\tAscendingNode<deg> {0:0.3f} # to match orbit\n'.format(
            elements.node))
        f.write('\t}\n')

    f.write('}\n\n')


data_dir = os.path.join(celestia_path, 'data')

stars_dat = os.path.join(data_dir, 'stars.dat')
with open(stars_dat, 'rb') as f:
    hip_stars = parse_stars_dat(f)


starnames_dat = os.path.join(data_dir, 'starnames.dat')
with open(starnames_dat, 'r') as f:
    star_names = parse_starnames(f)
    

hdxindex_dat = os.path.join(data_dir, 'hdxindex.dat')
with open(hdxindex_dat, 'rb') as f:
    star_names.update(parse_hdxindex(f))


parser = StcParser()
for filename in STC_FILES:
    stc_file = os.path.join(data_dir, filename)
    with open(stc_file, 'r', encoding='latin-1') as f:
        stc_hip, stc_names = parse_stc_file(f, parser)
        hip_stars |= stc_hip
        star_names.update(stc_names)

output_stars = set()

with open('exoplanets.csv', 'r') as f, \
        open('extrasolar.stc', 'w') as out_stars, \
        open('extrasolar.ssc', 'w') as out_planets:
    
    timestamp = time.strftime('%Y-%m-%dT%H:%M:%S%z', time.localtime())
    out_stars.write(STC_HEADER.format(timestamp))
    out_planets.write(SSC_HEADER.format(timestamp))
    
    csv = csv.DictReader(f)
    for row in csv:
        host_name = row["pl_hostname"]
        # in Celestia, use GJ
        if host_name[0:3] == 'GJ ':
            split_gliese = host_name.split(' ')
            if int(split_gliese[1]) < 1000:
                host_name = "Gliese " + host_name[3:]

        hip = None
        use_host = None
        
        if host_name in star_names:
            hip = star_names[host_name]
            use_host = host_name
            needs_star = (star_names[host_name] != -1 and
                          star_names[host_name] not in hip_stars)
        elif row["hd_name"] in star_names:
            hip = star_names[row["hd_name"]]
            use_host = row["hd_name"]
            needs_star = (star_names[row["hd_name"]] != -1 and
                          star_names[row["hd_name"]] not in hip_stars)
        elif row["hip_name"] != '':
            split_hip = row["hip_name"].split(' ')
            if len(split_hip) == 2 or split_hip[2] == 'A':
                hip = int(split_hip[1])
                use_host = "HIP " + str(hip)
                needs_star = hip not in hip_stars
            else:
                use_host = row["hip_name"]
                needs_star = True
        else:
            use_host = host_name
            needs_star = True

        if needs_star and not star_ok(row):
            print("Missing parameters for exoplanet host", row["pl_hostname"])
        elif not planet_ok(row):
            print("Missing parameters for exoplanet", row["pl_name"])
        else:
            if needs_star and row["pl_hostname"] not in output_stars:
                out_star(out_stars, row, hip, host_name, star_names)
                output_stars.add(row["pl_hostname"])
            out_planet(out_planets, row, use_host)