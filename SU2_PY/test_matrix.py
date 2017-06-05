#!/usr/bin/env python

## \file test_matrix.py
#  \brief show the test matrix of a group of test cases.
#  \author A. Mogavero
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

from SU2 import io
import shutil, os

def read_cases(root_dir):
    confs = []
    for (dirpath, dirnames, filenames) in os.walk(root_dir):
        conf_files = [f for f in filenames if '.cfg' == f[-4:]]
        for c in conf_files:
            confs.append(io.Config(filename=os.path.join(dirpath, c)))
                         
    return confs

if __name__=="__main__":
    
    confs = read_cases("/home/trb12187/Documents/kinetic/Tests/naca0012/")
    
    print(confs[0].diff(confs[1:]))
    