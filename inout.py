import numpy as np
import sys


def writeMesh(meshd, name):
        
    plex = meshd.mesh.topology._plex
    dim = meshd.mesh._topological_dimension
    cStart, cEnd = plex.getHeightStratum(0)
    eStart, eEnd = plex.getHeightStratum(1)
    vStart, vEnd = plex.getDepthStratum(0)
    nbrVer = vEnd - vStart
    
    coords = meshd.mesh.coordinates.dat.data
    
    f = open(name+".mesh", 'w')
    f.write("MeshVersionFormatted 1\n\nDimension 2\n\n")
    f.write("Vertices\n%d\n" % nbrVer)
    for iVer in range(nbrVer):
        off = meshd.section.getOffset(iVer+vStart)/dim
        f.write("%1.15e %1.15e 1\n" % (coords[off][0], coords[off][1]))
    f.write("Triangles\n%d\n" % (cEnd-cStart))
    for c in range(cStart,cEnd):
        closure = plex.getTransitiveClosure(c)[0]
        tri = [0,0,0]
        i = 0
        for cl in closure:
            if (cl < vStart or cl >= vEnd): continue
            tri[i] = cl-vStart+1
            i += 1
        f.write("%d %d %d 1\n" % (tri[0], tri[1], tri[2]))
    nbrEdg = plex.getStratumSize("exterior_facets", 1)
    f.write("Edges\n%d\n" % nbrEdg)
    for e in range(eStart, eEnd):
        l = plex.getLabelValue("boundary_ids", e)
        if l != -1:
            ver = plex.getCone(e)
            f.write("%d %d %d\n" % (ver[0]-vStart+1, ver[1]-vStart+1, l))    
    f.close()
    
def writeMetric(meshd, metric, name):
    
    plex = meshd.mesh.topology._plex
    dim = meshd.mesh._topological_dimension
    vStart, vEnd = plex.getDepthStratum(0)
    NbrVer = vEnd - vStart
    
    f = open(name+".sol", 'w')
    f.write("MeshVersionFormatted 2\n\nDimension 2\n\n")
    f.write("SolAtVertices\n%d\n1 3\n\n" % NbrVer)
    for iVer in range(NbrVer):
        off = meshd.section.getOffset(iVer+vStart)/dim
        f.write("%f %f %f\n" % (metric.dat.data[off][0][0], metric.dat.data[off][0][1], metric.dat.data[off][1][1]))
    f.write("\nEnd")
    f.close()
    
def writeMetric2(meshd, metric, name):

    plex = meshd.mesh.topology._plex
    dim = meshd.mesh._topological_dimension
    vStart, vEnd = plex.getDepthStratum(0)
    NbrVer = vEnd - vStart

    f = open(name+".sol", 'w')
    f.write("MeshVersionFormatted 2\n\nDimension 2\n\n")
    f.write("SolAtVertices\n%d\n1 3\n\n" % NbrVer)
    for iVer in range(NbrVer):
        f.write("%f %f %f\n" % (metric[iVer][0], metric[iVer][1], metric[iVer][2]))
    f.write("\nEnd")
    f.close()
    
def writeSol(meshd, u, name):

    plex = meshd.mesh.topology._plex
    dim = meshd.mesh._topological_dimension
    vStart, vEnd = plex.getDepthStratum(0)
    NbrVer = vEnd - vStart

    f = open(name+".sol", 'w')
    f.write("MeshVersionFormatted 2\n\nDimension 2\n\n")
    f.write("SolAtVertices\n%d\n1 1\n\n" % NbrVer)
    for iVer in range(NbrVer):
        off = meshd.section.getOffset(iVer+vStart)/dim
        f.write("%f\n" % (u.dat.data[off]))
    f.write("\nEnd")
    f.close()
    
    
def readMetric(meshd, metric, coordSection, name):
    
    plex = meshd.mesh.topology._plex
    dim = meshd.mesh._topological_dimension
    vStart, vEnd = plex.getDepthStratum(0)
    NbrVer = vEnd - vStart
    
    f = open(name+".sol", 'r')
    count = 0
    iVer = 0
    for line in f:
        count += 1
        if count < 9:
            continue
        h = line.split(" ")
        h1 = float(h[0])
        h2 = float(h[1])
        h3 = float(h[2])
        off = meshd.section.getOffset(iVer+vStart)/dim
        metric.dat.data[off] = [[h1, h2], [h2, h3]]
        iVer += 1
        if iVer == NbrVer:
            break
    f.close()