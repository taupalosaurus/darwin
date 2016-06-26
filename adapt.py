from firedrake import *
import numpy as np
from inout import *
from numpy import linalg as LA



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
    for H, mesh in zip(hessianMetrics, meshes) :
        
        V = FunctionSpace(mesh, 'CG', 1)    
        det = Function(V)
        # compute determinant
        for iVer in range(mesh.topology.num_vertices()):
            hesLoc = H.dat.data[iVer]
            detLoc = hesLoc[0,0]*hesLoc[1,1] - hesLoc[0,1]*hesLoc[1,0]
#            lbd, v = LA.eig(hesLoc)
#            lbd1 = max(abs(lbd[0]), 1e-10)
#            lbd2 = max(abs(lbd[1]), 1e-10)
#            detLoc = lbd1*lbd2
            H.dat.data[iVer] *= pow(detLoc, -1./(2*p+2))
            det.dat.data[iVer] = pow(detLoc, float(p)/(2*p+2))
        
        # intergrate determinant over space and assemble gloabl normalization term
        cofGlob += assemble(det*dx)
        
    cofGlob = N/cofGlob
    
    j = 0
    for H, mesh in zip(hessianMetrics, meshes) :
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
        
#        entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
#        entity_dofs[0] = mesh.geometric_dimension()
#        coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
#        writeMesh(mesh, coordSection, "metric.%d" % j)
#        writeMetric(mesh, H, coordSection, "metric.%d" % j)
        
    return hessianMetrics