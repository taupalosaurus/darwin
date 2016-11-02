from firedrake import *
#import numpy as np
#from numpy import linalg as LA
import sys, os, time, gc
from multiprocessing import Pool

from options import *
from mesh import *
from inout import *
from solve import *
from adapt import *
from interpol import *



def ptfx (options) :
   
    print "DEBUG  generating initial unit mesh"
    chrono1 = time.clock()
    if options.dim == 2 :
        mesh = Meshd(UnitSquareMesh(options.n, options.n))
    else :
        mesh = Meshd(UnitCubeMesh(options.n, options.n, options.n))
    chrono2 = time.clock()
    print "DEBUG  done generating initial mesh. Elapsed time: %1.2e" %(chrono2-chrono1); sys.stdout.flush()
    
    dtAdap = float(options.Tend)/options.nbrAdap


    for i in range(options.nbrPtfxIte) :

        print "\n\n\n##########  i: %d\n" % i ; sys.stdout.flush()

        j = 0
        globNormCoef = 0.

        for j in range(1, options.nbrAdap+1) :

            print "\n#####  j: %d\n" % j ; sys.stdout.flush()

            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            if i == 0:
                if j == 1 :
                    ####  Initial solution    
                    solIni = solIniAdvec(mesh)
                    writeGmf(mesh.mesh, 1, "boundary_ids", "bubble.0", solIni, 1, "bubble.0", mesh.section)
                else :
                    solIni = sol
            else :
                mesh = Meshd(readGmfMesh("newmesh.%d" %j, options.dim, "boundary_ids"))
                if j == 1 :
                    ####  Initial solution    
                    solIni = solIniAdvec(mesh)
                    writeGmf(mesh.mesh, 1, "boundary_ids", "bubble.0", solIni, 1, "bubble.0", mesh.section)
                else :
                    solIni = Function(FunctionSpace(mesh.mesh, 'CG', 1))
                    print "DEBUG  interpolation"; sys.stdout.flush()
                    chrono1 = time.clock()
                    interpol(sol, meshOld, solIni, mesh)
                    chrono2 = clock()
                    print "DEBUG  end interpolation. Elapsed time: %1.2e" %(chrono2-chrono1);  sys.stdout.flush()

            # solve 
            tIni, tEnd  = (j-1)*dtAdap, j*dtAdap
            sol, hesMet, t = solveAdvec(mesh, solIni, tIni, tEnd, options)
            globNormCoef = computeGlobalNormalizationCoef(options, mesh, hesMet, globNormCoef)
            print "DEBUG  globNormCoef: %1.18e" % globNormCoef

            writeGmf(mesh.mesh, 1, "boundary_ids", "bubble.%d" % j, sol, 1, "bubble.%d" % j, mesh.section)
            writeGmf(mesh.mesh, 0, "boundary_ids", "", hesMet, 5, "hesmet.%d" %j, mesh.section)
            if options.nbrSav > 0 :
                kIni = (j-1)*options.nbrSav
                for k in range(options.nbrSav+1) :
                    kGlob = kIni + k
                    print "DEBUG   film.%d -> film.%d" %(k, kGlob)
                    if k == 0 :
                        os.rename("film_tmp.%d.meshb" % k, "film.%d.meshb" % kGlob)
                    else :
                        if os.path.exists("film.%d.meshb" % kGlob) : os.remove("film.%d.meshb" % kGlob)
                        os.symlink("film.%d.meshb" % kIni, "film.%d.meshb" % kGlob)
                    os.rename("film_tmp.%d.solb" % k, "film.%d.solb" % kGlob)

            meshOld = mesh
            gc.collect()


        ######## End of the loop over sub-intervals
        
        if i < (options.nbrPtfxIte-1) :
        
            # normalizeMetric
            print "########## Metrics computation" ; sys.stdout.flush()
            chrono1 = time.clock()
            globNormCoef = 2.257682503752497780e+02
            metrics = normalizeUnsteadyMetrics(options, globNormCoef)
            print "########## End metrics computation. Elapsed time: %1.2e" %(chrono2-chrono1) ; sys.stdout.flush()
            gc.collect()
        
            # generate meshes
            print "########## Meshes generation" ; sys.stdout.flush()
            
            pool = Pool(15)
            pool.map(adaptInternal_star, zip(range(1,options.nbrAdap+1), [options.dim]*options.nbrAdap))
            gc.collect()



