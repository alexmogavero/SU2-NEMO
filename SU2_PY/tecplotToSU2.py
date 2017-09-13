#!/usr/bin/env python

import re
import os
import sys
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from SU2.mesh import tools
from copy import deepcopy
from math import sqrt
from scipy import spatial

def readZone(f, plotZone=False):
    param = {}
    zoneName = ''
    for l in f:
        if re.match('^ZONE', l):
            out = re.match('^.*T="(.*)"', l)
            zoneName = out.group(1)
            param['T'] = zoneName
            break
    if zoneName == '':
        raise StopIteration
        
    for l in f:
        out = l.split(',')
        for o in out:
            gr = re.match(' *(.*) *= *(.*) *', o)
            param[gr.group(1)] = gr.group(2)
        
        if 'DT' in param:
            break
    param['Nodes'] = long(param['Nodes'])
    param['Elements'] = long(param['Elements'])
    
    Coord = []
    i = 0
    prog = re.compile('^ *([0-9eE\-+.]*) *([0-9eE\-+.]*) *([0-9eE\-+.]*)')
    for l in f:
        out = prog.match(l)
        Coord.append([float(out.group(1)), float(out.group(2)), float(out.group(3))])

        i += 1
        if i==param['Nodes']:
            break
        
    i = 0
    elem = []
    for l in f:
        out = re.split(' *', l)
        out = out[1:]
        elem.append([long(o) for o in out])

        i += 1
        if i==param['Elements']:
            break
        
    
    if plotZone:
        fig = pl.figure()
        ax = fig.add_subplot(111, projection='3d')
        for e in elem:
            I = e + [e[0]]
            ax.plot([Coord[i-1][0] for i in I], [Coord[i-1][1] for i in I], [Coord[i-1][2] for i in I], '-b')
            
        pl.show()
        
    return param, Coord, elem

def febrickToVTK(elem):
    elvtk = []
    prevEl = -1
    for e in elem:
        if e != prevEl:
            elvtk.append(e)
        prevEl = e
    
    return elvtk
       
#     el1 = elem[0:4]
#     el2 = elem[4:]
#     n1 = len(set(el1))
#     n2 = len(set(el2))
#     
#     if n1 < n2:
#         eltmp = el1
#         n1 = n1tmp
#         el1 = el2
#         el2 = el1tmp
#         n1 = n2
#         n2 = n1tmp
#     
#     el1vtk = []
#     el2vtk = []
#     if n2 < 4:
#         if n1 < 4:
#             ePrev = -1
#             for e in el1:
#                 if e != ePrev:
#                     el1vtk.append(e)
#                 ePrev = e
#             
#         ePrev = -1
#         for e in el2:
#             if e != ePrev:
#                 el2vtk.append(e)
#             ePrev = e
#             
#         elvtk = el1vtk + el2vtk
#         if len(elvtk) != n1 + n2:
#             raise RuntimeError('Inconsistency in the element data.')
#     else:
#         elvtk = elem
#         
#     return elvtk

def convertTecplotSU2(tecData):
    data = {}
    data['NDIME'] = 3 #TODO Add handling of 2D meshes
    
    #internal mesh
    data['NPOIN'] = tecData[0][0]['Nodes']
    data['POIN'] = [tecData[0][1][i] + [i] for i in range(data['NPOIN'])]
    
    data['NELEM'] = tecData[0][0]['Elements']
    data['ELEM'] = []
    i = long(0)
    for e in tecData[0][2]:
        el = febrickToVTK([ee-1 for ee in e]) #list(set([ee-1 for ee in e]))
        if len(el) == 4:
            type = 10
        elif len(el) == 8:
            type = 12
        elif len(el) == 6:
            type = 13
        elif len(el) == 5:
            type = 14
        else:
            raise RuntimeError("Error: number of nodes not consistent with any element type.")
        data['ELEM'].append([type] + el + [i])
        i += 1
        
    # Buondary mesh
    print '\tCreating KD tree...'
    tree = spatial.KDTree([d[0:-1] for d in data['POIN']])
    
    data['NMARK'] = long(len(tecData) - 1)
    data['MARKS'] = {}
    points = [pInter[0:-1] for pInter in data['POIN']]
    values = [pInter[-1] for pInter in data['POIN']]
    for bDatTec in tecData[1:]:
        print '\t' + bDatTec[0]['T']
        
        data['MARKS'][bDatTec[0]['T']] = {'TAG':bDatTec[0]['T'], 'NELEM':bDatTec[0]['Elements'], 'ELEM':[]}
        
        # Find internal point id correspondent to the boundary point
        # In tecplot the boundary point are replicted and a new id is created
        id = [-1]*len(bDatTec[1])
        j = 0
        maxDis = 0
        iMaxDis = -1
        for p in bDatTec[1]:
            dis, i = tree.query(p)
            id[j] = data['POIN'][i][-1]
            if dis > maxDis:
                maxDis = dis
                iMaxDis = i

            j += 1
        if iMaxDis != -1:
            print '\t\tMax error=' + repr(maxDis) + ' on point ' + repr(iMaxDis)
            print '\t\tof coordinate ' + repr(data['POIN'][iMaxDis][0:-1])
                
        i = long(0)
        for e in bDatTec[2]:
            el = febrickToVTK([id[ee-1] for ee in e]) #list(set([id[ee-1] for ee in e]))
            if len(el) == 3:
                type = 5
            elif len(el) == 4:
                type = 9
            else:
                raise RuntimeError("Error: number of nodes not consistent with any element type.")
            data['MARKS'][bDatTec[0]['T']]['ELEM'].append([type] + el + [i])
            i += 1
    
    return data

if __name__=='__main__':

    if len(sys.argv) != 3:
        raise RuntimeError("2 arguments required\n" +
                           "Usage:\n\ttecplotToSU2 inputfile outputfile")
    
    print "Opening file " + sys.argv[1]
    f = open(sys.argv[1], 'r')
    
    print "Reading zones:"
    data = []
    for n in range(100):
        try:
            data.append(readZone(f))
            print "\t" + data[-1][0]['T']
        except StopIteration: 
            break
    
    print "Converting to SU2 format"
    dataSU2 = convertTecplotSU2(data)
    
    tools.write(sys.argv[2], dataSU2)
    print "Mesh written to file " + sys.argv[2]
    
        