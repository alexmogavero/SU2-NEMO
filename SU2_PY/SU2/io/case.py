#!/usr/bin/env python

## \file case.py
#  \brief python package for CFD test case handling
#  \author T. Lukaczyk, F. Palacios
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, re
import numpy as np
from .tools import *
from config_options import *
from config import Config

inf = 1.0e20


class Case(Config):
    """ config = SU2.io.Case(filename="")
        
        Starts a case class, an extension of 
        Config()
    """
    
    def __init__(self,*args,**kwarg):
        super(Case, self).__init__(*args, **kwarg)
        
        self.name = os.path.splitext(os.path.basename(self._filename))[0]
        self.root = os.path.dirname(self._filename) 
        
        self.commit = None
        if os.path.isdir(self.root):
            log_file = [f for f in os.listdir(self.root) if os.path.isfile(os.path.join(self.root, f))]
            log_file = [f for f in log_file if 'log' in f]
            
            log_file = open(os.path.join(self.root, log_file[0]), "r")
            for l in log_file:
                mt = re.match('.*Git commit: *([0-9a-z]*)', l)
                if mt:
                    self.commit = mt.group(1)
                    break
            log_file.close()
    
    def diff(self, konfig):
        out = super(Case, self).diff(konfig)
        if out.keys():
            out.name = [k.name for k in konfig]
            out.root = [k.root for k in konfig]
            out.commit = [k.commit for k in konfig]
            
        return out
        
        


