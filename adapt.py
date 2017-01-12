from firedrake import *
import numpy as np
from inout import *
from mesh import *
from numpy import linalg as LA

import time, gc

from hessian import *



def adaptInternal(j,dim) :
    print "DEBUG  read mesh and metric %d" %j; sys.stdout.flush()
    chrono1 = time.clock()
    mesh = Meshd(readGmfMesh("bubble.%d" %j, dim, "boundary_ids"), computeAltMin=False)
    metric = Function(TensorFunctionSpace(mesh.mesh, 'CG', 1))
    readGmfSol(mesh.mesh, metric, "metric.%d" %j, 5, mesh.section)
    chrono2 = time.clock()
    print "DEBUG  end read mesh and metric. Elapsed time: %1.2e" %(chrono2-chrono1);  sys.stdout.flush()
    print "##### Adap procedure started %d" %j ; sys.stdout.flush()
    chrono1 = time.clock()
    newmesh = adapt(mesh.mesh, metric)
    gc.collect()
    newmeshd = Meshd(newmesh)
    chrono2 = time.clock()
    print "##### Adap procedure completed %d. Elapsed time: %1.2e" %(j, chrono2-chrono1) ; sys.stdout.flush()
    print "DEBUG  write mesh"; sys.stdout.flush()
    chrono1 = time.clock()
    writeGmf(newmeshd.mesh, 1, "boundary_ids", "newmesh.%d" % j, None, None, None, newmeshd.section)
    chrono2 = time.clock()
    print "DEBUG  end write mesh. Elapsed time: %1.2e" %(chrono2-chrono1);  sys.stdout.flush()


def adaptInternal_star(j_dim):
    return adaptInternal(*j_dim)


def adaptInternal_kernel(meshd, metric) :
    
    newmesh = adapt(meshd.mesh, metric)
    newmeshd = Meshd(newmesh)
    
    return newmeshd



def computeGlobalNormalizationCoef(options, meshd, H, coef) :

    p = options.p
    if options.dim == 2 :
        lpPow1 = -1./(2*p+2)
        lpPow2 = float(p)/(2*p+2)
    else :
        lpPow1 = -1./(2*p+3)
        lpPow2 = float(p)/(2*p+3)

    mesh = meshd.mesh
    V = FunctionSpace(mesh, 'CG', 1)    
    detH = Function(V)
    # compute determinant
    detH.interpolate(det(H))
    detH_pow1 = np.power(detH.dat.data, lpPow1)
    H.dat.data[...] *= detH_pow1[:, np.newaxis, np.newaxis]
    detH.dat.data[...] = np.power(detH.dat.data, lpPow2)
    
    # intergrate determinant over space and assemble gloabl normalization term
    coef += assemble(detH*dx)

    return coef



def normalizeUnsteadyMetrics(options, coef) :

    N = options.N*options.nbrAdap
    a = options.a
    usa2 = 1./(a*a)
    hmin = options.hmin
    ushmin2 = 1./(hmin*hmin)
    hmax = options.hmax
    ushmax2 = 1./(hmax*hmax)
    if options.dim == 2 :
        lpPow = 1.
    else :
        lpPow = 2./3

    cofGlob = pow(float(N)/coef, lpPow)

    lbdMin = op2.Global(1, ushmax2, dtype=float);
    lbdMax = op2.Global(1, ushmin2, dtype=float);
    rat = op2.Global(1, usa2, dtype=float);
    
    for j in range(1,options.nbrAdap+1) :

        meshd = Meshd(readGmfMesh("bubble.%d" %j, options.dim, "boundary_ids"), reorderPlex=False, computeAltMin=False)
        H = Function(TensorFunctionSpace(meshd.mesh, 'CG', 1))
        readGmfSol(meshd.mesh, H, "hesmet.%d" %j, 5, meshd.section)
        H.dat.data[...] *= cofGlob
        op2.par_loop(options.absTruncMetric_kernel, H.node_set().superset, H.dat(op2.RW), lbdMin(op2.READ), lbdMax(op2.READ), rat(op2.READ))
        writeGmf(meshd.mesh, 0, "", "", H, 5, "metric.%d"%j, meshd.section)
        gc.collect()




def normalizeUnsteadyMetrics_old(options) :
    
    p = options.p
    N = options.N*options.nbrAdap
    a = options.a
    usa2 = 1./(a*a)
    hmin = options.hmin
    ushmin2 = 1./(hmin*hmin)
    hmax = options.hmax
    ushmax2 = 1./(hmax*hmax)
    if options.dim == 2 :
        lpPow1 = -1./(2*p+2)
        lpPow2 = float(p)/(2*p+2)
        lpPow3 = 1.
    else :
        lpPow1 = -1./(2*p+3)
        lpPow2 = float(p)/(2*p+3)
        lpPow3 = 2./3

    meshes = []
    hessianMetrics = []

    for j in range(1, options.nbrAdap+1) :
        meshd = Meshd(readGmfMesh("bubble.%d" %j, options.dim, "boundary_ids"), reorderPlex=False, computeAltMin=False)
        H = Function(TensorFunctionSpace(meshd.mesh, 'CG', 1))
        readGmfSol(meshd.mesh, H, "hesmet.%d" %j, 5, meshd.section)
        meshes.append(meshd)
        hessianMetrics.append(H)
    
    cofGlob = 0
    for H, meshd in zip(hessianMetrics, meshes) :
        
        mesh = meshd.mesh
        V = FunctionSpace(mesh, 'CG', 1)    
        detH = Function(V)
        # compute determinant
        detH.interpolate(det(H))
        detH_pow1 = np.power(detH.dat.data, lpPow1)
        H.dat.data[...] *= detH_pow1[:, np.newaxis, np.newaxis]
        detH.dat.data[...] = np.power(detH.dat.data, lpPow2)
        
        # intergrate determinant over space and assemble gloabl normalization term
        cofGlob += assemble(detH*dx)
        
    cofGlob = pow(float(N)/cofGlob, lpPow3)

    lbdMin = op2.Global(1, ushmax2, dtype=float);
    lbdMax = op2.Global(1, ushmin2, dtype=float);
    rat = op2.Global(1, usa2, dtype=float);
    
    j = 0
    for H, meshd in zip(hessianMetrics, meshes) :
        
        mesh = meshd.mesh
        j += 1
        H.dat.data[...] *= cofGlob

        op2.par_loop(options.absTruncMetric_kernel, H.node_set().superset, H.dat(op2.RW), lbdMin(op2.READ), lbdMax(op2.READ), rat(op2.READ))

        writeGmf(mesh, 0, "", "", H, 5, "metric.%d"%j, meshd.section)
        

    
    
    
    
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
#        print "####  ERROR option 2 for steady metric not implemented yet"
#        exit(1)
#    elif metOpt == 'Loseille2011' :
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


if __name__ == '__main__':

    parameters["pyop2_options"]["log_level"] = "WARNING"
    parameters["assembly_cache"]["enabled"] = False

    adaptInternal(4,3)
    adaptInternal(6,3)
