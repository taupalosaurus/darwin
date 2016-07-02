import numpy as np
from numpy import linalg as LA
import sys, os, time

from options import *
from mesh import *
from inout import *
from solve import *
from adapt import *
from interpol import *

t0 = time.clock()
from firedrake import *
tloadfdk = time.clock()
print "## module firedrake loaded (%ds)" % (tloadfdk - t0); sys.stdout.flush()

INF = float("inf")






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
        newmesh = adaptInternal(mesh, metric)
        print "##### Adap procedure finished"; sys.stdout.flush()
        writeMesh(newmesh, "newmesh")

        newsol = Function(newmesh.V)
        interpol(sol, mesh, newsol, newmesh)
        sol = newsol
        mesh = newmesh
        tIni = tEndSolve
        
        if nbrAdap == 0 :
            writeMesh(mesh, "bubble.0")
            writeSol(mesh, sol, "bubble.0")
            options.nbrGlobSav += 1
        
        nbrAdap += 1
        
            
    
    
    
            
if __name__ == '__main__':


    parameters["pyop2_options"]["log_level"] = "WARNING"

    tini = time.clock()

    main()

    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print "#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5])
