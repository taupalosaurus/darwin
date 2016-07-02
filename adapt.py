from firedrake import *
import numpy as np
from inout import *
from mesh import *
from numpy import linalg as LA



def adaptInternal(meshd, metric) :
    
    newmesh = adapt(meshd.mesh, metric)
    newmeshd = Meshd(newmesh)
    
    return newmeshd


def normalizeUnsteadyMetrics(hessianMetrics, meshes, options) :
    
    p = options.p
    N = options.N*options.nbrAdap
    a = options.a
    usa2 = 1./(a*a)
    hmin = options.hmin
    ushmin2 = 1./(hmin*hmin)
    hmax = options.hmax
    ushmax2 = 1./(hmax*hmax)
    
    cofGlob = 0
    for H, meshd in zip(hessianMetrics, meshes) :
        
        mesh = meshd.mesh
        V = FunctionSpace(mesh, 'CG', 1)    
        detH = Function(V)
        # compute determinant
        detH.interpolate(det(H))
        detH_pow1 = np.power(detH.dat.data, -1./(2*p+2))
        H.dat.data[...] *= detH_pow1[:, np.newaxis, np.newaxis]
        detH.dat.data[...] = np.power(detH.dat.data, float(p)/(2*p+2))

#        for iVer in range(mesh.topology.num_vertices()):
#            hesLoc = H.dat.data[iVer]
#            detLoc = hesLoc[0,0]*hesLoc[1,1] - hesLoc[0,1]*hesLoc[1,0]
#            H.dat.data[iVer] *= pow(detLoc, -1./(2*p+2))
#            det.dat.data[iVer] = pow(detLoc, float(p)/(2*p+2))
        
        # intergrate determinant over space and assemble gloabl normalization term
        cofGlob += assemble(detH*dx)
        
    #cofGlob *= options.Tend/options.nbrSpl
    cofGlob = float(N)/cofGlob
    
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
    
    
    
    
def computeSteadyMetric(meshd, hessian, sol, options) :
    
    p = options.p
    N = options.N*options.nbrAdap
    use = N
    umin = 0.01
    a = options.a
    usa2 = 1./(a*a)
    hmin = options.hmin
    ushmin2 = 1./(hmin*hmin)
    hmax = options.hmax
    ushmax2 = 1./(hmax*hmax)
    
    mesh = meshd.mesh
    metric = hessian
    
    
    if options.steadyMetric == 1:
      for iVer in range(mesh.topology.num_vertices()):
          hesLoc = hessian.dat.data[iVer] * 1/max(abs(sol.dat.data[iVer]), umin) * use
          meanDiag = 0.5*(hesLoc[0][1] + hesLoc[1][0]);
          hesLoc[0][1] = meanDiag
          hesLoc[1][0] = meanDiag
          lbd, v = LA.eig(hesLoc)
          # truncation of eigenvalues
          lbd1 = min(ushmin2, max(ushmax2,abs(lbd[0])))
          lbd2 = min(ushmin2, max(ushmax2,abs(lbd[1])))
          maxLbd = max(lbd1, lbd2)
          lbd1 = max(lbd1, usa2*maxLbd)
          lbd2 = max(lbd2, usa2*maxLbd)
          v1, v2 = v[0], v[1]
          # reconstruction of |Hu|
          metric.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
          metric.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
          metric.dat.data[iVer][1,0] = metric.dat.data[iVer][0,1]
          metric.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
    elif options.steadyMetric == 2:
        print "####  ERROR option 2 for steady metric not implemented yet"
        exit(1)
    elif metOpt == 'Loseille2011' :
        # computation of the metric of Loseille 2011
        detHes = Function(meshd.V)
        for iVer in range(mesh.topology.num_vertices()):
            hesLoc = hessian.dat.data[iVer]
            meanDiag = 0.5*(hesLoc[0][1] + hesLoc[1][0]);
            hesLoc[0][1] = meanDiag
            hesLoc[1][0] = meanDiag
            lbd, v = LA.eig(hesLoc)
            lbd1 = max(abs(lbd[0]), 1e-10)
            lbd2 = max(abs(lbd[1]), 1e-10)
            det = lbd1*lbd2
            v1 = v[0]
            v2 = v[1]
            metric.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
            metric.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
            metric.dat.data[iVer][1,0] = metric.dat.data[iVer][0,1]
            metric.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
            metric.dat.data[iVer] *= pow(det, -1./(2*p+2))
            detHes.dat.data[iVer] = pow(det, p/(2.*p+2))
    
        detHesInt = assemble(detHes*dx)
        glbNrmCof = (N/detHesInt)
        metric *= glbNrmCof
    
        for iVer in range(mesh.topology.num_vertices()):
            lbd, v = LA.eig(metric.dat.data[iVer])
            # truncation of eigenvalues
            lbd1 = min(ushmin2, max(ushmax2,lbd[0]))
            lbd2 = min(ushmin2, max(ushmax2,lbd[1]))
            maxLbd = max(lbd1, lbd2)
            lbd1 = max(lbd1, usa2*maxLbd)
            lbd2 = max(lbd2, usa2*maxLbd)
            # reconstruction of |Hu|
            v1 = v[0]
            v2 = v[1]
            metric.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
            metric.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
            metric.dat.data[iVer][1,0] = metric.dat.data[iVer][0,1]
            metric.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
    else :
        print "####  ERROR ivalid option %d for steady metric" % options.steadyMetric
        exit(1)
    return metric
