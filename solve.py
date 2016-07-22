from firedrake import *
import pyop2 as op2
import sys, os, os.path
from numpy import linalg as LA
import numpy as np
from time import clock

from inout import *
from hessian import *





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

    print "DEBUG  start dt computation"; sys.stdout.flush()
    chrono1 = clock()

    dt = (meshd.altMin.dat.data / (cn.dat.data+1e-10)).min()
    dt *= options.cfl

    chrono2 = clock()
    print "DEBUG  end dt computation. Elapsed time : %1.2es" % (chrono2-chrono1); sys.stdout.flush()

    return dt

    
    
def computeAvgHessian(meshd, sol, t, tIni, tEnd, nbrSpl, H, hessian, options) :
    
    print "DEBUG  hessian-metric assembly"; sys.stdout.flush()
    chrono1 = clock()

    mesh = meshd.mesh

    detMax = Function(meshd.V).interpolate(abs(det(H))).dat.data.max()   
    lbdmax = sqrt(detMax)
    lbdmin = 1.e-20 * lbdmax

    lbdMin = op2.Global(1, lbdmin, dtype=float);
    op2.par_loop(options.absValHessian_kernel, H.node_set().superset, H.dat(op2.RW), lbdMin(op2.READ))

    # sum with previous hessians
    # 3 cases: if t=tIni or tEnd or in between
    cof = float(tEnd-tIni)/(nbrSpl-1)
    if (t == tIni) or (t == tEnd) : cof *= 0.5
    hessian.dat.data[...] += cof*H.dat.data

    chrono2 = clock()
    print "DEBUG  end Hessian-metric assembly. Elapsed time: %1.2e" %(chrono2-chrono1); sys.stdout.flush()




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
    if options.dim == 2:
        cxExpr = "-sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])" 
        cyExpr = "sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])"
        c0.interpolate(Expression([cxExpr, cyExpr]))
    else:
        cxExpr = "2*sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*sin(2*pi*x[2])"
        cyExpr = "-sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*sin(2*pi*x[2])"
        czExpr = "-sin(2*pi*x[0])*sin(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2])"
        c0.interpolate(Expression([cxExpr, cyExpr, czExpr]))
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
    
    sigma = TestFunction(M)
    H = Function(M)
    n = FacetNormal(mesh)
    Lh = inner(sigma, H)*dx + inner(div(sigma), grad(u0))*dx
    if options.dim == 2:
        Lh -= (sigma[0, 1]*n[1]*u0.dx(0) + sigma[1, 0]*n[0]*u0.dx(1))*ds
    else :
        Lh -= (sigma[0,0]*u0.dx(0)*n[0] + sigma[1,0]*u0.dx(1)*n[0] + sigma[2,0]*u0.dx(2)*n[0] \
            +  sigma[0,1]*u0.dx(0)*n[1] + sigma[1,1]*u0.dx(1)*n[1] + sigma[2,1]*u0.dx(2)*n[1] \
            +  sigma[0,2]*u0.dx(0)*n[2] + sigma[1,2]*u0.dx(1)*n[2] + sigma[2,2]*u0.dx(2)*n[2])*ds
        print "## Warning: Variational form for 3D hessian computation still under heavy work"    
    H_prob = NonlinearVariationalProblem(Lh, H)
    H_solv = NonlinearVariationalSolver(H_prob, solver_parameters={'snes_rtol': options.snes_rtol,
                                                                   'ksp_rtol': options.ksp_rtol,
                                                                   'ksp_gmres_restart': 20,
                                                                   'pc_type': 'sor',
                                                                   'snes_monitor': True,
                                                                   'snes_view': False,
                                                                   'ksp_monitor_true_residual': True,
                                                                   'snes_converged_reason': True,
                                                                   'ksp_converged_reason': True})

    nbrSpl = options.nbrSpl 
    dtSpl = float(tEnd-tIni)/(nbrSpl-1)

    if options.nbrSav > 0 :
        dtSav = float(tEnd-tIni)/(options.nbrSav)
   
    step = 0
    stepSav = 0
    stepSpl = 0
    print "\n####  step %d " % step ; sys.stdout.flush()
    chronostep1 = clock()
    
    t = tIni
    dt = 0
    u0.assign(solIni)

    if options.algo == 1:
        print "DEBUG  start hessian solve"; sys.stdout.flush()
        chrono1 = clock()
        H_solv.solve()
        chrono2 = clock()
        print "DEBUG  end hessian solve. Elapsed time: %1.2e" %(chrono2-chrono1); sys.stdout.flush()
        computeAvgHessian(meshd, u0, t, tIni, tEnd, nbrSpl, H, hessian, options) 

    if options.nbrSav > 0 :
        stepSpl = 0
        print "DEBUG  write mesh and solution"; sys.stdout.flush()
        chrono1 = clock()
        writeGmf(meshd.mesh, 1, "boundary_ids", "film_tmp.%d" % stepSav, u0, 1, "film_tmp.%d" % stepSav, meshd.section)
        chrono2 = clock()
        print "DEBUG  end write mesh and solution. Elapsed time: %1.2e" %(chrono2-chrono1);  sys.stdout.flush()

    
    while t < tEnd:
    
        chronostep1 = clock()
        step += 1
        print "\n####  step %d " % step ; sys.stdout.flush()
    
        time.assign(t)
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
        if (t < options.nbrGlobSav*options.dtSav) and (t+dt >= options.nbrGlobSav*options.dtSav) :
            print "DEBUG  Trunc dt for solution saving"
            # truncation of the dt for hessian sampling
            dt = options.nbrGlobSav*options.dtSav - t + 1.e-5*dt
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
        doSav2 = 0
        if (t < options.nbrGlobSav*options.dtSav) and (t+dt >= options.nbrGlobSav*options.dtSav) :
            doSav2 = 1
        doSpl = 0
        if ((t < tIni+(stepSpl+1)*dtSpl) and (t+dt >= tIni+(stepSpl+1)*dtSpl) or endSol) :
            stepSpl += 1
            doSpl = 1
        
        print "DEBUG  t:  %1.3e -> %1.3e with dt: %1.7e" %(t, t+dt, dt);  sys.stdout.flush()
        
        print "DEBUG  start solve"; sys.stdout.flush()
        chrono1 = clock()
        solve(a==L, u0, bcs=bc)
        chrono2 = clock()
        print "DEBUG  end solve. Elapsed time: %1.2e" %(chrono2-chrono1); sys.stdout.flush()

        t += dt
        
        if doSpl :
            print "DEBUG  start hessian solve"; sys.stdout.flush()
            chrono1 = clock()
            H_solv.solve()
            chrono2 = clock()
            print "DEBUG  end hessian solve. Elapsed time: %1.2e" %(chrono2-chrono1); sys.stdout.flush()
            computeAvgHessian(meshd, u0, t, tIni, tEnd, nbrSpl, H, hessian, options) 

        if doSav :
            if os.path.exists("film_tmp.%d.meshb" % stepSav) : os.remove("film_tmp.%d.meshb" % stepSav)
            os.symlink("film_tmp.0.meshb", "film_tmp.%d.meshb" % stepSav)
            print "DEBUG  write solution"; sys.stdout.flush()
            chrono1 = clock()
            writeGmf(meshd.mesh, 0, "", "", u0, 1, "film_tmp.%d" % stepSav, meshd.section)
            chrono2 = clock()
            print "DEBUG  end write solution. Elapsed time: %1.2e" %(chrono2-chrono1);  sys.stdout.flush()
        if doSav2 :
            writeGmf(meshd.mesh, 1, "boundary_ids", "bubble.%d" % options.nbrGlobSav, u0, 1, "bubble.%d" % options.nbrGlobSav, meshd.section)
            options.nbrGlobSav += 1
        if options.algo == 2 and (endSol or step == options.adaptStepFreq) :
            H_solv.solve()
            hessian = H
            break

        chronostep2 = clock()
        print "DEBUG  .Step elapsed time: %1.2e" %(chronostep2-chronostep1);  sys.stdout.flush()
            
    return [u0, hessian, t]



def solIniAdvec(meshd):
    
    V = FunctionSpace(meshd.mesh, "CG", 1)
    if meshd.mesh._topological_dimension == 2:
        icExpr = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.75)*(x[1]-0.75) < 0.15*0.15)")
    else: 
        icExpr = Expression("((x[0]-0.35)*(x[0]-0.35) + (x[1]-0.35)*(x[1]-0.35) + (x[2]-0.35)*(x[2]-0.35) < 0.15*0.15)")
    u0 = Function(V).interpolate(icExpr)
    
    return u0



def hessIniAdvec(meshd, sol):
    
    mesh = meshd.mesh

    M = TensorFunctionSpace(mesh, "CG", 1)

    sigma = TestFunction(M)
    H = Function(M)
    n = FacetNormal(mesh)
    Lh = inner(sigma, H)*dx
    Lh += inner(div(sigma), grad(sol))*dx - (sigma[0, 1]*n[1]*sol.dx(0) + sigma[1, 0]*n[0]*sol.dx(1))*ds
    H_prob = NonlinearVariationalProblem(Lh, H)
    H_solv = NonlinearVariationalSolver(H_prob)
    
    H_solv.solve()
    
    return H

    
    
    
