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


from __future__ import print_function

import os.path
import struct
import re
import ply.lex as lex
import ply.yacc as yacc
from io import open

reserved_words = {
    'Barycenter': 'BARYCENTER',
    'Add': 'DISPOSITION',
    'Modify': 'DISPOSITION',
    'Replace': 'DISPOSITION',
    }

STC_TOKENS = (
    'DISPOSITION',
    'BARYCENTER',
    'PROPERTY',
    'LBRACE',
    'RBRACE',
    'LBRACKET',
    'RBRACKET',
    'LANGLE',
    'RANGLE',
    'NUMERIC',
    'QUOTE',
    'UNESCAPED',
    'ESCAPE',
    'BACKSLASH',
    'NEWLINE',
    'UNICODE',
    )

class StcLexer(object):
    """Lexer for STC files.
    
    This class should be used in conjuction with the StcParser.
    """
    
    def __init__(self, **kwargs):
        """Create a new StcLexer with the given keyword arguments."""
        self.lexer = lex.lex(module=self, **kwargs)

    # List of token names
    tokens = STC_TOKENS + ('COMMENT',)

    # List of states
    states = (
        ('string', 'exclusive'),
        ('escape', 'exclusive'),
        )
    
    t_LBRACE = r'\{'
    t_RBRACE = r'\}'
    t_LBRACKET = r'\['
    t_RBRACKET = r'\]'
    t_LANGLE = r'<'
    t_RANGLE = r'>'

    def t_PROPERTY(self, t):
        r'[A-Z][A-Za-z]+'
        t.type = reserved_words.get(t.value, 'PROPERTY')
        return t

    def t_NUMERIC(self, t):
        r'[-+]?[0-9]+(\.[0-9]+)?([Ee][-+]?[0-9]+)?'
        if '.' in t.value or 'E' in t.value or 'e' in t.value:
            t.value = float(t.value)
        else:
            t.value = int(t.value)
        return t
    
    def t_QUOTE(self, t):
        r'\"'
        t.lexer.push_state('string')
        return t
    
    t_ignore_COMMENT = r'\#.*'
    t_ignore = ' \t\r\n'
    
    t_string_UNESCAPED = r'[\x20-\x21,\x23-\x5B,\x5D-\xFF]+'
    
    def t_string_ESCAPE(self, t):
        r'\\'
        t.lexer.push_state('escape')
        return t
    
    def t_string_QUOTE(self, t):
        r'\"'
        t.lexer.pop_state()
        return t
    
    t_string_ignore = ''
    
    def t_escape_BACKSLASH(self, t):
        r'\\'
        t.value = '\\'
        t.lexer.pop_state()
        return t
    
    def t_escape_NEWLINE(self, t):
        r'n'
        t.value = '\n'
        t.lexer.pop_state()
        return t
    
    def t_escape_QUOTE(self, t):
        r'\"'
        t.value = '"'
        t.lexer.pop_state()
        return t
    
    def t_escape_UNICODE(self, t):
        r'u[0-9a-fA-F]{4}'
        t.value = chr(int(t.value[1:5], 16))
        t.lexer.pop_state()
        return t
    
    t_escape_ignore = ''
    
    def t_error(self, t):
        print("Illegal token '{0}'".format(t.value))
        t.lexer.skip(len(t.value))
        
    def t_string_error(self, t):
        print("Illegal character in string '{0}'".format(t.value[0]))
        t.lexer.skip(1)
    
    def t_escape_error(self, t):
        print("Unknown string escape '{0}'".format(t.value[0]))
        t.lexer.pop_state()
        t.lexer.skip(1)


class StcParser(object):
    """Parser for STC files.
    
    This class should be used in conjunction with the StcLexer.
    
    Example usage:
       parser = StcParser()
       starlist = parser.parse(data)
    """
    
    def __init__(self, lexer=None, **kwargs):
        """Create a new StcParser with option to supply an external lexer."""
        if lexer is None:
            self.lexer = StcLexer().lexer
        else:
            if isinstance(lexer, StcLexer):
                self.lexer = lexer.lexer
            else:
                self.lexer = lexer
        self.parser = yacc.yacc(module=self, **kwargs)

    tokens = STC_TOKENS
    
    def p_stars_blank(self, p):
        '''stars :'''
        p[0] = list()
    
    def p_stars_star(self, p):
        '''stars : stars star
                 | stars DISPOSITION star'''
        if len(p) == 3:
            p[2]['Disposition'] = None
            p[2]['IsBarycenter'] = False
            p[1].append(p[2])
        else:
            p[3]['Disposition'] = p[2]
            p[3]['IsBarycenter'] = False
            p[1].append(p[3])
        p[0] = p[1]
            
    def p_stars_barycenter(self, p):
        '''stars : stars BARYCENTER star
                 | stars DISPOSITION BARYCENTER star'''
        if len(p) == 4:
            p[3]['Disposition'] = None
            p[3]['IsBarycenter'] = True
            p[1].append(p[3])
        else:
            p[3]['Disposition'] = p[2]
            p[3]['IsBarycenter'] = True
            p[1].append(p[3])
        p[0] = p[1]

    def p_star_numeric(self, p):
        '''star : NUMERIC block'''
        p[2]['HIP'] = p[1]
        p[2]['Name'] = None
        p[0] = p[2]
    
    def p_star_named(self, p):
        '''star : qstring block'''
        p[2]['HIP'] = None
        p[2]['Name'] = p[1]
        p[0] = p[2]
    
    def p_star_namednumeric(self, p):
        '''star : NUMERIC qstring block'''
        p[3]['HIP'] = p[1]
        p[3]['Name'] = p[2]
        p[0] = p[3]
    
    def p_block(self, p):
        '''block : LBRACE rules RBRACE'''
        p[0] = p[2]
    
    def p_rules(self, p):
        '''rules :
                 | rules rule'''
        if len(p) == 1:
            p[0] = {}
        else:
            p[1].update(p[2])
            p[0] = p[1]
    
    def p_rule_simple(self, p):
        '''rule : PROPERTY value
                | PROPERTY block'''
        p[0] = { p[1]: p[2] }
    
    def p_rule_units(self, p):
        '''rule : PROPERTY units value'''
        p[0] = { p[1]: p[3], p[1]+"_units": p[2] }
    
    def p_units(self, p):
        '''units : LANGLE PROPERTY RANGLE'''
        p[0] = p[2]

    def p_value(self, p):
        '''value : NUMERIC
                 | qstring
                 | vector'''
        p[0] = p[1]
    
    def p_qstring(self, p):
        '''qstring : QUOTE string QUOTE'''
        p[0] = p[2]
    
    def p_string(self, p):
        '''string :
                  | string characters'''
        if len(p) == 1:
            p[0] = ''
        else:
            p[1] += p[2]
            p[0] = p[1]
    
    def p_characters_unescaped(self, p):
        '''characters : UNESCAPED'''
        p[0] = p[1]
        
    def p_characters_escaped(self, p):
        '''characters : ESCAPE BACKSLASH
                      | ESCAPE NEWLINE
                      | ESCAPE QUOTE
                      | ESCAPE UNICODE'''
        p[0] = p[2]
    
    def p_vector(self, p):
        '''vector : LBRACKET elements RBRACKET'''
        p[0] = p[2]
    
    def p_elements(self, p):
        '''elements :
                    | elements NUMERIC'''
        if len(p) == 1:
            p[0] = list()
        else:
            p[1].append(p[2])
            p[0] = p[1]
    
    def parse(self, data, lexer=None, *args, **kwargs):
        """Parses textual data into a list of dictionaries."""
        if lexer is None:
            lexer = self.lexer
        return self.parser.parse(data, lexer=lexer, *args, **kwargs)


class StarDatabase(object):
    """Manages the star database"""
    
    def __init__(self):
        """Creates a new star database"""
        self.hip_stars = set()
        self.star_names = {}
    
    def parse_files(self, file_list):
        """Parse the given file list"""
        for filename in file_list:
            basename = os.path.basename(filename)
            if basename == 'hdxindex.dat':
                self.parse_hdxindex(filename)
            elif basename == 'starnames.dat':
                self.parse_starnames(filename)
            elif basename == 'stars.dat':
                self.parse_stars_dat(filename)
            else:
                extension = os.path.splitext(filename)[-1]
                if extension == '.stc':
                    self.parse_stc_file(filename)
    
    def parse_stars_dat(self, filename):
        """Extracts the HIP numbers from the stars.dat file."""
        with open(filename, 'rb') as f:
            f.seek(14)
            data = f.read()
            self.hip_stars |= {
                struct.unpack('<I', data[pos:pos+4])[0]
                for pos in range(0, len(data), 20)
                }


    def parse_starnames(self, filename):
        """Extracts the star names from starnames.dat"""
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                split_line = line.split(':')
                hip = int(split_line[0])
                self.star_names.update({name: (hip, False)
                                        for name in split_line[1:]})


    def parse_hdxindex(self, filename):
        """Extracts the HD names from hdxindex.dat"""
        with open(filename, 'rb') as f:
            f.seek(10)
            data = f.read()
            for pos in range(0, len(data), 8):
                hd, hip = struct.unpack("<II", data[pos:pos+8])
                self.star_names["HD "+str(hd)] = hip, True


    def parse_stc_file(self, filename, parser=None):
        """Extracts the HIP numbers from an stc file."""
        if parser is None:
            parser = StcParser()

        with open(filename, 'r', encoding='latin-1') as f:
            stc_stars = parser.parse(f.read())
            for star in stc_stars:
                if star["HIP"] is not None:
                    self.hip_stars.add(star["HIP"])
                if star["Name"] is not None:
                    split_names = star["Name"].split(':')
                    self.star_names.update({name: (-1, False)
                                            for name in split_names})

    
    def lookup_name(self, name):
        """Looks up a name in the star database."""
        if name in self.star_names:
            return self.star_names[name][0]
        m = re.match('^HIP ([0-9]+)$', name)
        if m:
            return int(m.group(1))
        else:
            return None
    
    def is_catalog_name(self, name):
        """Checks if the name in question is from the cross-index"""
        if name in self.star_names:
            return self.star_names[name][1]
        m = re.match('^HIP ([0-9]+)$', name)
        if m:
            return True
        else:
            return False