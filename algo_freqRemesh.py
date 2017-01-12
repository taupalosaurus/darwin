from firedrake import *
import numpy as np
from numpy import linalg as LA
import sys, os, time

from options import *
from mesh import *
from inout import *
from solve import *
import lswe
from adapt import *
from interpol import *







def freqRemesh(options) :
                    

    mesh = Meshd(UnitSquareMesh(options.n, options.n))
    
    tEndSolve = 0
    nbrAdap = 0
    tIni = 0;
        
    while (True) :
        
        print "\n\n\n##########  i: %d\n" % nbrAdap ; sys.stdout.flush()        
        
        if nbrAdap == 0 :
            sol = solIniAdvec(mesh)
            hessian = hessIniAdvec(mesh, sol)
            tEndSolve = 0
        else :
            sol, hessian, tEndSolve = solveAdvec(mesh, sol, tIni, options.Tend, options)
            
        if tEndSolve >= options.Tend :
            break
        
        metric = computeSteadyMetric(mesh, hessian, sol, options)

        print "##### Adap procedure started "; sys.stdout.flush()
        newmesh = adaptInternal_kernel(mesh, metric)
        print "##### Adap procedure finished"; sys.stdout.flush()
        writeGmf(newmesh.mesh, 1, "boundary_ids", "newmesh", None, None, None, newmesh.section)

        newsol = Function(newmesh.V)
        interpol(sol, mesh, newsol, newmesh)
        sol = newsol
        mesh = newmesh
        tIni = tEndSolve
        
        if nbrAdap == 0 :
            writeGmf(mesh.mesh, 1, "boundary_ids", "bubble.0", sol, 1, "bubble.0", mesh.section)
            options.nbrGlobSav += 1
        
        nbrAdap += 1
        
            

def freqRemesh_swe(options) :


    # 1000 km x 500km domain with 20x20 resolution
    mesh = Meshd(RectangleMesh(options.n, options.n/2, 1e6, 5e5))
#    mesh = Meshd(RectangleMesh(2, 2, 1e6, 5e5))
    delta_t = 15.0
    dt = Constant(delta_t) # time step

    xy = SpatialCoordinate(mesh.mesh)
    depth = conditional(xy[0]<5e5, 1e3, Max(1e3 * (1.0-(xy[0]-5e5)/5e5), 10.))

    model = lswe.LinearShallowWaterSolver(mesh.mesh, depth, dt)

    u0, eta0 = model.get_solution()

    x0, y0 = 5e5, 2.5e5 # centre of gaussian bump
    L = 5e4 # width of gaussian bump

    eta0.interpolate(exp(-(pow((xy[0]-x0)/L,2)+pow((xy[1]-y0)/L,2))))

    tEndSolve = 0
    nbrAdap = 0
    tIni = 0;

        
    while (True) :
        
        print "\n\n\n##########  i: %d\n" % nbrAdap ; sys.stdout.flush()        
        
#        writeGmf(mesh.mesh, 1, "boundary_ids", "swe_iniite.%d" % nbrAdap, eta0, 1, "swe_iniite.%d" % nbrAdap, mesh.section)
#        writeGmf(mesh.mesh, 1, "boundary_ids", "u_iniite.%d" % nbrAdap, u0, 2, "u_iniite.%d" % nbrAdap, mesh.section)

        if nbrAdap == 0 :
            tEndSolve = 0
        else :
            for i in range(options.adaptStepFreq):
                print "solver iteration ", i
                model.iterate()
                #outFile.write(eta0, u0)
            tEndSolve += options.adaptStepFreq*delta_t    
            
        File("out.%d.pvd" % nbrAdap).write(eta0, u0)

        if tEndSolve >= options.Tend :
            break

#        writeGmf(mesh.mesh, 1, "boundary_ids", "swe_before.%d" % nbrAdap, eta0, 1, "swe_before.%d" % nbrAdap, mesh.section)
#        writeGmf(mesh.mesh, 1, "boundary_ids", "u_before.%d" % nbrAdap, u0, 2, "u_before.%d" % nbrAdap, mesh.section)

        
        model.compute_hessian()
        hessian = model.get_hessian()

        metric = computeSteadyMetric(mesh, hessian, eta0, options)

        print "##### Adap procedure started "; sys.stdout.flush()
        newmesh = adaptInternal_kernel(mesh, metric)
        #newmesh = Meshd(RectangleMesh(options.n, options.n/2, 1e6, 5e5))
        print "##### Adap procedure finished"; sys.stdout.flush()
        writeGmf(newmesh.mesh, 1, "boundary_ids", "newmesh", None, None, None, newmesh.section)

        xy = SpatialCoordinate(newmesh.mesh)
        depth = conditional(xy[0]<5e5, 1e3, Max(1e3 * (1.0-(xy[0]-5e5)/5e5), 10.))
        newmodel = lswe.LinearShallowWaterSolver(newmesh.mesh, depth, dt)

        newu0, neweta0 = newmodel.get_solution()

        print "DEBUG  here"
        u0_proj = Function(VectorFunctionSpace(mesh.mesh, 'CG', 1)).interpolate(u0)
        newu0_proj = Function(VectorFunctionSpace(newmesh.mesh, 'CG', 1))
        interpol(u0_proj, mesh, newu0_proj, newmesh)
        newu0.interpolate(newu0_proj)
        print "DEBUG  there"
#        eta0_proj = Function(FunctionSpace(mesh.mesh, 'CG', 1)).interpolate(eta0)
#        neweta0_proj = Function(FunctionSpace(newmesh.mesh, 'CG', 1))
#        interpol(eta0_proj, mesh, neweta0_proj, newmesh)
#        neweta0.interpolate(neweta0_proj)
        interpol(eta0, mesh, neweta0, newmesh)
        model = newmodel
        u0, eta0 = model.get_solution()
        mesh = newmesh
        
        writeGmf(mesh.mesh, 1, "boundary_ids", "swe_after.%d" % nbrAdap, interpolate(eta0,FunctionSpace(mesh.mesh, 'CG', 1)), 1, "swe_after.%d" % nbrAdap, mesh.section)
        writeGmf(mesh.mesh, 1, "boundary_ids", "u_after.%d" % nbrAdap, interpolate(u0,VectorFunctionSpace(mesh.mesh, 'CG', 1)), 2, "u_after.%d" % nbrAdap, mesh.section)
        
        nbrAdap += 1
    
    
            
if __name__ == '__main__':


    parameters["pyop2_options"]["log_level"] = "WARNING"

    tini = time.clock()

    main()

    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print "#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5])
