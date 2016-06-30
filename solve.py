from firedrake import *
import sys, os, os.path
from numpy import linalg as LA
import numpy as np

from inout import *





def computeDtAdvecOld(meshd, cn, cxExpr, cyExpr, options) :
    
    nrmc = "sqrt((%s)*(%s)+(%s)*(%s))" % (cxExpr, cxExpr, cyExpr, cyExpr)
    cn.interpolate(Expression(nrmc))
    cmax = cn.dat.data.max()
    hmin = 1./options.n #### TODO don't leave that !!!
    dt = 0.5*hmin/cmax

    return dt
    

    
def computeDtAdvec(meshd, cn, cxExpr, cyExpr, options) :
    
#    dt = 1e10
#    for iVer in range(meshd.mesh.topology.num_vertices()):
#        dtloc = meshd.altMin.dat.data[iVer] / (cn.dat.data[iVer]+1e-10)
#        dt = min(dt, dtloc)   

    dt = (meshd.altMin.dat.data / (cn.dat.data+1e-10)).min()
    dt *= options.cfl

    return dt

    
    
def computeAvgHessian(meshd, sol, t, tIni, tEnd, nbrSpl, H, hessian) :
    
    print "DEBUG  hessian-metric assembly"
    
    mesh = meshd.mesh
    
#    # compute the hessian
#    # computation of the hessian
#    sigma = TestFunction(M)
#    H = Function(M)
#    n = FacetNormal(mesh)
#    L = inner(sigma, H)*dx
#    L += inner(div(sigma), grad(sol))*dx - (sigma[0, 1]*n[1]*sol.dx(0) + sigma[1, 0]*n[0]*sol.dx(1))*ds
#    H_prob = NonlinearVariationalProblem(L, H)
#    H_solv = NonlinearVariationalSolver(H_prob)
#    H_solv.solve()

    detMax = Function(meshd.V).interpolate(abs(det(H))).dat.data.max()   
#    detMax = 0
#    for iVer in range(mesh.topology.num_vertices()):
#        detLoc = H.dat.data[iVer][0,0]*H.dat.data[iVer][1,1] + H.dat.data[iVer][0,1]*H.dat.data[iVer][1,0]
#        detMax = max(detMax, abs(detLoc))

    lbdmax = sqrt(detMax)
    lbdmin = 1.e-20 * lbdmax
    
    # take the absolute value of the hessian
    for iVer in range(mesh.topology.num_vertices()):
        hesLoc = H.dat.data[iVer]
        meanDiag = 0.5*(hesLoc[0][1] + hesLoc[1][0]);
        hesLoc[0][1] = meanDiag
        hesLoc[1][0] = meanDiag
        lbd, v = LA.eig(hesLoc)
        lbd1 = max(abs(lbd[0]), lbdmin)
        lbd2 = max(abs(lbd[1]), lbdmin)
#        det = lbd1*lbd2
        v1, v2 = v[0], v[1]
        H.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
        H.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
        H.dat.data[iVer][1,0] = H.dat.data[iVer][0,1]
        H.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];

    # sum with previous hessians
    # 3 cases if t=tIni or tEnd or between
    cof = float(tEnd-tIni)/(nbrSpl-1)
    if (t == tIni) or (t == tEnd) : cof *= 0.5
    hessian.dat.data[...] += cof*H.dat.data





def solveAdvec(meshd, solIni, tIni, tEnd, options):
    
    mesh = meshd.mesh
    
    T = options.T
    
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 2)
    Vn = FunctionSpace(mesh, "CG", 1)
    M = TensorFunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(Q), TestFunction(Q)
    u0 = Function(Q)     
    c0, cn = Function(V), Function(Vn) # velocity base and its norm
    hessian = Function(M)


            
    zero = Constant(0)
    bc = DirichletBC(Q, zero, 1)

    time = Constant(0)
    delta_t = Constant(0.1)
    cxExpr = "-sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])" 
    cyExpr = "sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])"
    c0.interpolate(Expression([cxExpr, cyExpr]))
    c = c0*cos(2*pi*time/T)

    a = v*u*dx + 0.5*delta_t*(v*dot(c, grad(u))*dx)
    L = v*u0*dx - 0.5*delta_t*(v*dot(c, grad(u0))*dx)
    # Add SUPG stabilization terms
    ra = u + 0.5*delta_t*(dot(c, grad(u)))
    rL = u0 - 0.5*delta_t*(dot(c, grad(u0)))
    cnorm = sqrt(dot(c, c))
    h = meshd.altMin
    a += (h/(2.*cnorm))*dot(c, grad(v))*ra*dx
    L += (h/(2.*cnorm))*dot(c, grad(v))*rL*dx


    nbrSpl = options.nbrSpl 
    dtSpl = float(tEnd-tIni)/(nbrSpl-1)

    if options.nbrSav > 0 :
        dtSav = float(tEnd-tIni)/(options.nbrSav)
   
    step = 0
    stepSav = 0
    stepSpl = 0
    print "####  step %d " % step ; sys.stdout.flush()
    
    t = tIni
    dt = 0
    u0.assign(solIni)

    sigma = TestFunction(M)
    H = Function(M)
    n = FacetNormal(mesh)
    Lh = inner(sigma, H)*dx
    Lh += inner(div(sigma), grad(u0))*dx - (sigma[0, 1]*n[1]*u0.dx(0) + sigma[1, 0]*n[0]*u0.dx(1))*ds
    H_prob = NonlinearVariationalProblem(Lh, H)
    H_solv = NonlinearVariationalSolver(H_prob)
    H_solv.solve()
    computeAvgHessian(meshd, u0, t, tIni, tEnd, nbrSpl, H, hessian) 

    if options.nbrSav > 0 :
        stepSpl = 0
        writeMesh(meshd, "film_tmp.%d" % stepSpl)
        writeSol(meshd, u0, "film_tmp.%d" % stepSpl)
    
    while t < tEnd:
    
        step += 1
        print "####  step %d " % step ; sys.stdout.flush()
    
        time.assign(t)
#        c.interpolate(c0*cos(2*pi*t/T))
#        cnorm = sqrt(dot(c, c))
        cn.interpolate(cnorm)
    
        dt = computeDtAdvec(meshd, cn, cxExpr, cyExpr, options)
        if (options.nbrSav > 0) and (t < tIni+(stepSav+1)*dtSav) and (t+dt >= tIni+(stepSav+1)*dtSav) :
            print "DEBUG  Trunc dt for solution saving"
            # truncation of the dt for solution saving
            dt = (tIni+(stepSav+1)*dtSav) - t + 1.e-5*dt
        if (t < tIni+(stepSpl+1)*dtSpl) and (t+dt >= tIni+(stepSpl+1)*dtSpl) :
            print "DEBUG  Trunc dt for hessian sampling"
            # truncation of the dt for hessian sampling
            dt = (tIni+(stepSpl+1)*dtSpl) - t + 1.e-5*dt
        endSol = 0
        if (t+dt > tEnd) : 
            print "DEBUG  Trunc dt to final time"
            endSol = 1
            dt = tEnd - t + + 1.e-5*dt
        delta_t.assign(dt)
        doSav = 0
        if (options.nbrSav > 0) and ((t < tIni+(stepSav+1)*dtSav) and (t+dt >= tIni+(stepSav+1)*dtSav) or endSol):
            stepSav += 1
            doSav = 1
        doSpl = 0
        if ((t < tIni+(stepSpl+1)*dtSpl) and (t+dt >= tIni+(stepSpl+1)*dtSpl) or endSol) :
            stepSpl += 1
            doSpl = 1
        
        print "DEBUG  t:  %1.3e -> %1.3e with dt: %1.7e" %(t, t+dt, dt)
        
#        a = v*u*dx + 0.5*delta_t*(v*dot(c, grad(u))*dx)
#        L = v*u0*dx - 0.5*delta_t*(v*dot(c, grad(u0))*dx)
#        # Add SUPG stabilization terms
#        ra = u + 0.5*delta_t*(dot(c, grad(u)))
#        rL = u0 - 0.5*delta_t*(dot(c, grad(u0)))
##        cnorm = sqrt(dot(c, c))
#        h = meshd.altMin
#        a += (h/(2.*cnorm))*dot(c, grad(v))*ra*dx
#        L += (h/(2.*cnorm))*dot(c, grad(v))*rL*dx
        solve(a==L, u0, bcs=bc)
        
        t += dt
        
        if doSpl :
            H_solv.solve()
            computeAvgHessian(meshd, u0, t, tIni, tEnd, nbrSpl, H, hessian) 

        if doSav :
            #writeMesh(meshd, "film.%d" % stepSav)
            if os.path.exists("film_tmp.%d.mesh" % stepSav) : os.remove("film_tmp.%d.mesh" % stepSav)
            os.symlink("film_tmp.0.mesh", "film_tmp.%d.mesh" % stepSav)
            writeSol(meshd, u0, "film_tmp.%d" % stepSav)
            
            
    return [u0, hessian]
    
    
    
def solIniAdvec(meshd):
    
    Q = FunctionSpace(meshd.mesh, "CG", 1)
    icExpr = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.75)*(x[1]-0.75) < 0.15*0.15)")
    u0 = Function(Q).interpolate(icExpr)
    
    return u0
    
    
    