#!/usr/bin/env python

import re
import os
import sys
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from SU2.mesh import tools
import getopt
from math import sqrt
from scipy import spatial

typName = {3:'Line', 5:'Triangle', 9:'Quadrilateral', 10:'Tetrahedral', 12:'Hexahedral', 13:'Wedge', 14:'Pyramid'};
nodeNum = {3:2, 5:3, 9:4, 10:4, 12:8, 13:6, 14:5};

def repair(el):
    print 'attempting repair of element: ' + str(el[-1])
    
    seen = set()
    unique = []
    for nd in el[1:-1]:
        if nd not in seen:
            seen.add(nd)
            unique.append(nd)
            
    newType = [nodeNum[t] for t in nodeNum if nodeNum[t]==len(unique)]
    if len(newType) > 0:
        newType = newType[0]
        el[0:] = [newType] + unique + el[-1]
        
        return []
        
    else:
        print 'No suitable type found for the collapsed element, triangulating...'
        
        if el[0] == 12: #Hexahedral
            triang = [[1,2,3,6],
                      [1,3,4,8],
                      [1,3,6,8],
                      [3,6,7,8],
                      [1,5,6,8]]
        else:
            raise RuntimeError('Triangulation not implemented for element type ' + typName[el[0]])
        
        newEl = []
        for t in triang:
            nEl = [el[tt] for tt in t]
            if len(set(nEl)) == 4:
                newEl.append([10] + nEl + [-1])
                
        newEl[0][-1] = el[-1]
        el[0:] = newEl[0]
                
        return newEl[1:]
        

if __name__=='__main__':
    # inputs
    rep = False
    outputfile = None
    try:
        opts, args = getopt.getopt(sys.argv[1:-1],"hr:o:",["help", "repair", "ofile="])
    except getopt.GetoptError:
        print 'checkMesh.py [[-r | --repair] [-o | --ofile] <outputfile>] <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
           print 'checkMesh.py [[-r | --repair] [-o | --ofile] <outputfile>] <inputfile>'
           sys.exit()
        elif opt in ("-r", "--repair"):
           rep = True
        elif opt in ("-o", "--ofile"):
           outputfile = arg

    inputfile = sys.argv[-1]
    print 'Reading file...'
    data = tools.read(inputfile)
    
    typCnt = {3:0, 5:0, 9:0, 10:0, 12:0, 13:0, 14:0};
    for el in data['ELEM']:
        if el[0] not in typName:
            raise RuntimeError('Element type:' + str(el[0]) + ' unrecognized.')
        
        typCnt[el[0]] += 1
        if len(el) - 2 != nodeNum[el[0]]:
            print 'Wrong number of nodes for element ' + str(el[-1]) + ' of type ' + typName[el[0]]  
            
        unqNode = set(el[1:-1])
        if len(unqNode) != nodeNum[el[0]]:
            print 'Element ' + str(el[-1]) + ' has collapsed nodes.'
            if rep:
                newEl = repair(el)
                for nEl in newEl:
                    nEl[-1] = data['NELEM']
                    data['ELEM'].append(nEl)
                    data['NELEM'] += 1
                    
            print ' '
                            
    for tp in typCnt:
        print typName[tp] + ' ' + str(typCnt[tp])
    print ' '
     
    if outputfile:
        print 'Writing repaired mesh...'
        tools.write(outputfile, data)
        