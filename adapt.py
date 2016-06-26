from firedrake import *
import numpy as np
from inout import *
from mesh import *
from numpy import linalg as LA



def adaptInternal(meshd, metric) :
    
    newmesh = adapt(meshd.mesh, metric)
    newmeshd = Meshd(newmesh)
    
    return newmeshd


def normalizeMetrics(hessianMetrics, meshes, options) :
    
    p = options.p
    N = options.N
    a = options.a
    usa2 = 1/(a*a)
    hmin = options.hmin
    ushmin2 = 1/(hmin*hmin)
    hmax = options.hmax
    ushmax2 = 1/(hmax*hmax)
    
    cofGlob = 0
    for H, meshd in zip(hessianMetrics, meshes) :
        
        mesh = meshd.mesh
        V = FunctionSpace(mesh, 'CG', 1)    
        det = Function(V)
        # compute determinant
        for iVer in range(mesh.topology.num_vertices()):
            hesLoc = H.dat.data[iVer]
            detLoc = hesLoc[0,0]*hesLoc[1,1] - hesLoc[0,1]*hesLoc[1,0]
            H.dat.data[iVer] *= pow(detLoc, -1./(2*p+2))
            det.dat.data[iVer] = pow(detLoc, float(p)/(2*p+2))
        
        # intergrate determinant over space and assemble gloabl normalization term
        cofGlob += assemble(det*dx)
        
    cofGlob = N/cofGlob
    
    j = 0
    for H, meshd in zip(hessianMetrics, meshes) :
        
        mesh = meshd.mesh
        j += 1
        H.dat.data[...] *= cofGlob
        
        for iVer in range(mesh.topology.num_vertices()):
            lbd, v = LA.eig(H.dat.data[iVer])
            # truncation of eigenvalues
            lbd1 = min(ushmin2, max(ushmax2,lbd[0]))
            lbd2 = min(ushmin2, max(ushmax2,lbd[1]))
            maxLbd = max(lbd1, lbd2)
            lbd1 = max(lbd1, usa2*maxLbd)
            lbd2 = max(lbd2, usa2*maxLbd)
            v1, v2 = v[0], v[1]
            # reconstruction of |Hu|
            H.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
            H.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
            H.dat.data[iVer][1,0] = H.dat.data[iVer][0,1]
            H.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
        
#        writeMesh(meshd, "metric.%d" % j)
#        writeMetric(meshd, H, "metric.%d" % j)
        
    return hessianMetrics