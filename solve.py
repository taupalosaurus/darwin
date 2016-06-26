from firedrake import *
import sys
from numpy import linalg as LA
import numpy as np





def computeDtAdvec(mesh, cn, cxExpr, cyExpr):
    
    n = 50
    
    nrmc = "sqrt((%s)*(%s)+(%s)*(%s))" % (cxExpr, cxExpr, cyExpr, cyExpr)
    cn.interpolate(Expression(nrmc))
    cmax = cn.dat.data.max()
    hmin = 1./n #### TODO don't leave that !!!
    dt = 0.5*hmin/cmax

    return dt
    
    
    
def computeAvgHessian(mesh, sol, t, tIni, tEnd, nbrSpl, M, hessian) :
    
    print "DEBUG  computing hessian sample"
    
    # compute the hessian
    # computation of the hessian
    sigma = TestFunction(M)
    H = Function(M)
    n = FacetNormal(mesh)
    L = inner(sigma, H)*dx
    L += inner(div(sigma), grad(sol))*dx - (sigma[0, 1]*n[1]*sol.dx(0) + sigma[1, 0]*n[0]*sol.dx(1))*ds
    H_prob = NonlinearVariationalProblem(L, H)
    H_solv = NonlinearVariationalSolver(H_prob)
    H_solv.solve()
    
    # take the absolute value of the hessian
    for iVer in range(mesh.topology.num_vertices()):
        hesLoc = H.dat.data[iVer]
        meanDiag = 0.5*(hesLoc[0][1] + hesLoc[1][0]);
        hesLoc[0][1] = meanDiag
        hesLoc[1][0] = meanDiag
        lbd, v = LA.eig(hesLoc)
        lbd1 = max(abs(lbd[0]), 1e-10)
        lbd2 = max(abs(lbd[1]), 1e-10)
        det = lbd1*lbd2
        v1, v2 = v[0], v[1]
        H.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
        H.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
        H.dat.data[iVer][1,0] = H.dat.data[iVer][0,1]
        H.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];

    # sum with previous hessians
    # 3 cases if t=tIni or tEnd or between
    cof = (tEnd-tIni)/(nbrSpl-1)
    if (t == tIni) or (t == tEnd) : cof *= 0.5
    hessian.dat.data[...] += cof*H.dat.data





def solveAdvec(mesh, solIni, tIni, tEnd):
    
    T = 6
    
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 2)
    Vn = FunctionSpace(mesh, "CG", 2)
    M = TensorFunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(Q), TestFunction(Q)
    u0 = Function(Q)     
    c, cn = Function(V), Function(Vn) # velocity and its norm
    hessian = Function(M)
            
    zero = Constant(0)
    bc = DirichletBC(Q, zero, 1)

    nbrSpl = 5 
    dtSpl = (tEnd-tIni)/(nbrSpl-1)
    
            
    step = 0
    print "####  step %d " % step ; sys.stdout.flush()
    
    t = tIni
    dt = 0
    u0.assign(solIni)
    
    computeAvgHessian(mesh, u0, t, tIni, tEnd, nbrSpl, M, hessian) 
    
    while t < tEnd:
    
        step += 1
        print "####  step %d " % step ; sys.stdout.flush()
    
        cxExpr = "-sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*cos(2*pi*%f)" % (t/T)
        cyExpr = "sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(2*pi*%f)" % (t/T)
        c.interpolate(Expression([cxExpr, cyExpr]))
    
        dt = computeDtAdvec(mesh, cn, cxExpr, cyExpr)
        endSol = 0
        if (t+dt > tEnd) : 
            print "Trunc end run"
            endSol = 1
            dt = tEnd - t + + 1.e-5*dt
        doSpl = 0
        if (not endSol) and (int(t/dtSpl) != int((t+dt)/dtSpl)) :
            print "trunc sample"
            # truncation of the dt for hessian sampling
            stepSpl = int(t/dtSpl) + 1
            dt = stepSpl*dtSpl - t + 1.e-5*dt
            doSpl = 1
        
        print "DEBUG  t:  %1.3e -> %1.3e with dt: %1.7e" %(t, t+dt, dt)
        a = v*u*dx + 0.5*dt*(v*dot(c, grad(u))*dx)
        L = v*u0*dx - 0.5*dt*(v*dot(c, grad(u0))*dx)
        solve(a==L,u0)
        
        t += dt
        
        if doSpl or endSol :
            computeAvgHessian(mesh, u0, t, tIni, tEnd, nbrSpl, M, hessian) 
            
            
    return [u0, hessian]
    
    
    
def solIniAdvec(mesh):
    
    Q = FunctionSpace(mesh, "CG", 1)
    icExpr = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.75)*(x[1]-0.75) < 0.15*0.15)")
    u0 = Function(Q).interpolate(icExpr)
    
    return u0