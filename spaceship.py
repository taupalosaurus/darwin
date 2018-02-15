from firedrake import *
import pyop2 as op2

import sys, os, os.path
from numpy import linalg as LA
import numpy as np
from time import clock

from inout import *
from hessian import *
from options import *
from mesh import *



    
def computeDtSpaceship(meshd, uv, options) :
    
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




def solveSpaceship(meshd, uv_init, elev_init, tIni, tEnd, iteration, options):
    
    mesh2d = meshd.mesh
        
    # Solver for hessian computation
    M = TensorFunctionSpace(mesh2d, 'CG', 1)
    hessian = Function(M)  
    sigma = TestFunction(M)
    H = Function(M)
    n = FacetNormal(mesh2d)
    Lh = inner(sigma, H)*dx + inner(div(sigma), grad(u0))*dx
    Lh -= (sigma[0, 1]*n[1]*u0.dx(0) + sigma[1, 0]*n[0]*u0.dx(1))*ds
    H_prob = NonlinearVariationalProblem(Lh, H)
    H_solv = NonlinearVariationalSolver(H_prob, solver_parameters={'snes_rtol': options.snes_rtol,
                                                                   'ksp_rtol': options.ksp_rtol,
                                                                   'ksp_gmres_restart': 20,
                                                                   'pc_type': 'sor',
                                                                   'snes_monitor': True,
                                                                   'snes_view': False,
                                                                   'ksp_monitor_true_residual': False,
                                                                   'snes_converged_reason': True,
                                                                   'ksp_converged_reason': True})

    x = mesh2d.coordinates
    # bathymetry
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    # bathymetry2d = Function(P1_2d, name="bathymetry")
    bathymetry2d = project(Min(25.0 - 20.0*(x[0]-20000.)/(31500.-20000.), 25.0), P1_2d)
    
    bottom_drag = Constant(0.0025)
    hViscosity = Constant(100.0)

    # TODO  COMPUTE DT SOMEWHERE, and make sure it is compatible with nbrSav and nbrSpl
    dt = 20

    # --- create solver ---
    solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
    options_solv = solverObj.options
    options_solv.use_nonlinear_equations = True
    options_solv.simulation_export_time =  (tEnd-tIni)/options.nbrSav
    options_solv.simulation_end_time = tEnd - 1e-3
    options_solv.output_directory = outputDir+".%d" % i_adapt
    options_solv.check_volume_conservation_2d = True
    options_solv.fields_to_export = ['uv_2d', 'elev_2d']
    options_solv.timestepper_type = 'PressureProjectionPicard'
    options_solv.timestepper_options.implicitness_theta = 1.0
    options_solv.use_lax_friedrichs_velocity = False
    options_solv.quadratic_drag_coefficient = bottom_drag
    options_solv.horizontal_viscosity = hViscosity
    options_solv.element_family = 'dg-cg'
    options_solv.timestep = dt 


    # ---- Boundary conditions ------------
    # boundary conditions
    inflow_tag = 2
    elev_func=Function(P1_2d)
    elev_bc = {'elev': elev_func}
    solverObj.bnd_functions['shallow_water'] = {inflow_tag: elev_bc}


    nbrSpl = options.nbrSpl 
    dtSpl = float(tEnd-tIni)/(nbrSpl-1)

    if options.nbrSav > 0 :
        dtSav = float(tEnd-tIni)/(options.nbrSav)
   
    step = 0
    stepSpl = 0
    
    deltaT = 100. # of forcing data
    alpha = Constant(0)
    forcing_cur = Constant(0)
    forcing_next = Constant(0)

    def updateForcings(t):
        t -= dt/2.
        i = int(t/deltaT)
        alpha.assign(t/deltaT-i)
        forcing_cur.assign(forcing[i])
        forcing_next.assign(forcing[i+1])
        elev_func.assign((1.0-alpha)*forcing_cur + alpha*forcing_next)
        print_output("FORCING: {}, {}".format(t, (1.0-alpha.dat.data[0])*forcing[i] + alpha.dat.data[0]*forcing[i+1]))

        t += dt/2
        j = (t-tIni)/dtSpl
        if  abs(j-int(j))/j < 0.0001: #TODO does that really make sense ? trying to guess is (t-tIni) is divided by dtSpl      
            H_solv.solve()
            computeAvgHessian(meshd, u0, t, tIni, tEnd, nbrSpl, H, hessian, options)
            stepSpl += 1
    
#    if (i_adapt == 0) :   #TODO  do at the right place !!!
#        updateForcings(0.)
    
    solverObj.assign_initial_conditions(uv=uv_init, elev=elev_init)

    # Update solver time parameters to start at the right time
    if (i_adapt > 0):
        solverObj.i_export = i_adapt*options.nbrSav 
        solverObj.iteration = iteration_prev
        solverObj.next_export_t = tIni + (tEnd-tIni)/nbrSav
        solverObj.simulation_time = tIni
        for e in solverObj.exporters.values():
            e.set_next_export_ix(solverObj.i_export)

#    if options.algo == 1:
#        H_solv.solve()
#        computeAvgHessian(meshd, u0, t, tIni, tEnd, nbrSpl, H, hessian, options) 

    solverObj.iterate(update_forcings=updateForcings)
    uv_sol, elev_sol = solverObj.timestepper.solution.split()

    iteration = solverObj.iteration
    if abs(solverObj.simulation_time - tEnd)/tEnd > 0.0001 :
        print("ERROR  tend != tend: %1.3e %1.3e\n", solverObj.simulation_time, tEnd)
        sys.exit(1)

    print "DEBUG  End solve"
            
    return [uv_sol, elev_sol, hessian, t, iteration]



def solIniSpaceship(meshd):
    
    uv_init = Expression(("1e-8","0.0"))
    elev_init = Expression("0.")
    
    return uv_init, elev_init



def hessIniSpaceship(meshd, sol):
    
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

