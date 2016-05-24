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

# Acknowledgements:
# Thanks goes to Grant Hutchison, whose orbit transform spreadsheets were
# invaluable in checking the output of these converters.


import numpy as np
import numpy.linalg as linalg
from numpy import sin, cos, arccos, arctan2, radians, degrees

OBLIQUITY = 23.4392911

def rotation_x(angle):
    c = cos(radians(angle))
    s = sin(radians(angle))
    return np.matrix((
        (1, 0, 0),
        (0, c, -s),
        (0, s, c),
        ))


def rotation_y(angle):
    c = cos(radians(angle))
    s = sin(radians(angle))
    return np.matrix((
        (c, 0, s),
        (0, 1, 0),
        (-s, 0, c),
        ))


def rotation_z(angle):
    c = cos(radians(angle))
    s = sin(radians(angle))
    return np.matrix((
        (c, -s, 0),
        (s, c, 0),
        (0, 0, 1),
        ))


class OrbitElements(object):
    """Contains orbital elements"""

    def __init__(self, arg_peri, inclination, node):
        if arg_peri is not None:
            while arg_peri < 0:
                arg_peri += 360
            while arg_peri > 360:
                arg_peri -= 360
        self.arg_peri = arg_peri
        
        self.inclination = inclination
        
        if node is not None:
            while node < 0:
                node += 360
            while node > 360:
                node -= 360
        self.node = node


    def __repr__(self):
        return str.format(
            'OrbitElements(arg_peri={0!r}, inclination={1!r}, node={2!r})',
            self.arg_peri,
            self.inclination,
            self.node,
            )


class OrbitConverter(object):
    """Provides transformations for orbital elements between the sky plane and
    the ecliptic frames."""
    
    def __init__(self, ra, dec):
        self._ra = ra
        self._dec = dec
                
        self._transform = (
            rotation_x(-OBLIQUITY) *
            rotation_z(ra) *
            rotation_y(-dec - 90),
            )


    def transform(self, elements):
        """Converts a set of orbital elements from the sky plane to the
        ecliptic frame."""
        sky_matrix = (
            rotation_z(elements.node) *
            rotation_x(-elements.inclination) * 
            rotation_z(elements.arg_peri)
            )

        ecl_matrix = self._transform * sky_matrix
        
        m22 = ecl_matrix[2, 2]
        if m22 < -1:
            m22 = -1
        elif m22 > 1:
            m22 = 1
            
        inclination = degrees(arccos(m22))

        if abs(ecl_matrix[2, 2] == 1.0):
            arg_peri = degrees(arctan2(-ecl_matrix[0, 1], ecl_matrix[0, 0]))
            node = 0.0
        else:
            arg_peri = degrees(arctan2(ecl_matrix[2, 0], ecl_matrix[2, 1]))
            node = degrees(arctan2(ecl_matrix[0, 2], -ecl_matrix[1, 2]))
        return OrbitElements(arg_peri, inclination, node)


    def ecliptic_plane(self):
        """Gets the inclination and node of the ecliptic plane in the sky
        frame."""
        v = linalg.inv(self._transform).dot(np.array((0, 0, 1)))
        inclination = degrees(arccos(v[0, 2]))
        node = degrees(arctan2(v[0, 0], -v[0, 1])) + 180
        if node > 360:
            node -= 360
        return OrbitElements(None, inclination, node)


    @property
    def ra(self):
        return _ra


    @property
    def dec(self):
        return _dec


    def __repr__(self):
        return str.format(
            'OrbitConverter(ra={0!r}, dec={1!r})',
            self._ra,
            self._dec,
            )