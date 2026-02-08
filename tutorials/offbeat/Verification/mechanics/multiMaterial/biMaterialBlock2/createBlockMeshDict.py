#!/bin/env python3

from ofblockmeshdicthelper import BlockMeshDict, Vertex, ArcEdge, SimpleGrading
import numpy as np
from math import pi, sin, cos, radians, atan

class BlockNamingIter():
    def __init__(self):
        self.count = 0
    
    def __call__(self):
        self.count += 1
        return ('b%d' % (self.count-1))

class VertexNamingIter():
    def __init__(self):
        self.count = 0
    
    def __call__(self):
        self.count += 1
        return ('v%d' % (self.count-1))

class ArcNamingIter():
    def __init__(self):
        self.count = 0
    
    def __call__(self):
        self.count += 1
        return ('arc%d' % (self.count-1))

bNamer = BlockNamingIter()
vNamer = VertexNamingIter()
arcNamer = ArcNamingIter()

class NewVertex(Vertex):
    def __init__(self, x, y, z, name=None):
        if name is None:
            name = vNamer()
        Vertex.__init__(self, x, y, z, name)
        self.index = vNamer.count
    
    def rotZ(self, theta):
        q = radians(theta)
        c = cos(q)
        s = sin(q)
        return NewVertex(c*self.x - s*self.y, s*self.x + c*self.y, self.z)
    
    def __add__(self, v):
        if isinstance(v, Vertex):
            return NewVertex(self.x + v.x, self.y + v.y, self.z + v.z)
        else:
            return NewVertex(self.x + v[0], self.y + v[1], self.z + v[2])
    
    def __repr__(self):
        return ('Vertex %s (%g,%g,%g)' % (self.name, self.x, self.y, self.z))

class NewArcEdge(ArcEdge):
    def __init__(self, v1, v2, interVertex, name=None):
        if name is None:
            name = arcNamer()
        ArcEdge.__init__(self, [v1.name, v2.name], name, interVertex)
    
    def __repr__(self):
        return ('Arc %s (%d,%d) through (%g %g %g)' % (self.name,
                self.vnames[0], self.vnames[0],
                self.interVertex.x, self.interVertex.y, self.interVertex.z))

class UserBlockMeshDict(BlockMeshDict):
    def addVertex(self, v):
        return self.add_vertex(v.x, v.y, v.z, v.name)
    def addArc(self, arc):
        return self.add_arcedge(arc.vnames, arc.name, arc.interVertex)
    def addHexBlock(self, vertices, cells, name=None, grading=None):
        kwargs = {}
        if name is None:
            name = bNamer()
        if not grading is None:
            kwargs.update(dict(grading=grading))
        
        vnames = [v.name for v in vertices]
        return self.add_hexblock(vnames, cells, name, **kwargs)



bmd = UserBlockMeshDict()
bmd.set_metric('cm')

va = [NewVertex(x, 0, 0) for x in [-5, 0, 5]]
theta = 0
vb = [(v + (4*sin(theta), 4*cos(theta), 0)) for v in va]

vc = [v + (0,0,0.1) for v in va]
vd = [v + (0,0,0.1) for v in vb]

for v in va+vb+vc+vd:
    bmd.addVertex(v)

blockNames = ['left', 'right']
blocks = []
for name, vertices in zip(blockNames, zip(va[:-1], va[1:], vb[1:], vb[:-1], 
                                          vc[:-1], vc[1:], vd[1:], vd[:-1])):
    blocks.append(bmd.addHexBlock(vertices, (40, 20, 1), name=name))

frontAndBackFaces = [face for b in blocks for face in [b.face('b'), b.face('t')]]
topFaces = [b.face('n') for b in blocks]
bottomFaces = [b.face('s') for b in blocks]
bmd.add_boundary('empty', 'frontAndBack', frontAndBackFaces)
bmd.add_boundary('patch', 'left', [blocks[0].face('w')])
bmd.add_boundary('patch', 'right', [blocks[-1].face('e')])
bmd.add_boundary('patch', 'top', topFaces)
bmd.add_boundary('patch', 'bottom', bottomFaces)

bmd.assign_vertexid()
print(bmd.format())

with open('system/blockMeshDict', 'w') as f:
    f.write(bmd.format())


   