from firedrake import *
#import numpy as np
#from numpy import linalg as LA
import sys, os, time, gc
from multiprocessing import Pool

from options import *
from mesh import *
from solve import *
from adapt import *
from interpol import *
from spaceship import *



def ptfx (options) :
   
    print("DEBUG  reading initial unit mesh")
    mesh = Meshd(Mesh('symmetric_b.msh'))
    
    dtAdap = float(options.Tend)/options.nbrAdap

    iteration = 0


    for i in range(options.nbrPtfxIte) :

        print("\n\n\n##########  i: %d\n" % i ); sys.stdout.flush()

        j = 0
        globNormCoef = 0.

        for j in range(1, options.nbrAdap+1) :

            print("\n#####  j: %d\n" % j ); sys.stdout.flush()

            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            if i == 0:
                if j == 1 :
                    ####  Initial solution    
                    uv_init, elev_init = solIniSpaceship(mesh)
#                    writeGmf(mesh.mesh, 1, "boundary_ids", "uv.0", uv_init, 2, "uv.0", mesh.section)
                else :
                    uv_init, elev_init = uv, elev
            else :
                mesh = Meshd(readGmfMesh("newmesh.%d" %j, options.dim, "boundary_ids"))
                if j == 1 :
                    ####  Initial solution    
                    uv_init, elev_init = solIniSpaceship(mesh)
#                    writeGmf(mesh.mesh, 1, "boundary_ids", "uv.0", uv_init, 2, "uv.0", mesh.section)
                else :
                    print("DEBUG  interpolation"); sys.stdout.flush()
                    uv_init, elev_init = transfer_solution(mesh.mesh, uv, elev)
                    print("DEBUG  end interpolation.");  sys.stdout.flush()

            # solve 
            tIni, tEnd  = (j-1)*dtAdap, j*dtAdap
            uv, elev, hesMet, t, iteration = solveSpaceship(mesh, uv_init, elev_init, tIni, tEnd, iteration, j-1, options)
            globNormCoef = computeGlobalNormalizationCoef(options, mesh, hesMet, globNormCoef)
            print("DEBUG  globNormCoef: %1.18e" % globNormCoef)

#            writeGmf(mesh.mesh, 1, "boundary_ids", "uv.%d" % j, uv, 2, "uv.%d" % j, mesh.section)
            writeGmf(mesh.mesh, 1, b"boundary_ids", b"hesmet.%d" %j, hesMet, 5, b"hesmet.%d" %j, mesh.section)
#            if options.nbrSav > 0 :
#                kIni = (j-1)*options.nbrSav
#                for k in range(options.nbrSav+1) :
#                    kGlob = kIni + k
#                    print("DEBUG   film.%d -> film.%d" %(k, kGlob))
#                    if k == 0 :
#                        os.rename("film_tmp.%d.meshb" % k, "film.%d.meshb" % kGlob)
#                    else :
#                        if os.path.exists("film.%d.meshb" % kGlob) : os.remove("film.%d.meshb" % kGlob)
#                        os.symlink("film.%d.meshb" % kIni, "film.%d.meshb" % kGlob)
#                    os.rename("film_tmp.%d.solb" % k, "film.%d.solb" % kGlob)

            meshOld = mesh
            gc.collect()


        ######## End of the loop over sub-intervals
        
        if i < (options.nbrPtfxIte-1) :
        
            # normalizeMetric
            print("########## Metrics computation") ; sys.stdout.flush()
            metrics = normalizeUnsteadyMetrics(options, globNormCoef)
            print("########## End metrics computation."); sys.stdout.flush()
            gc.collect()
        
            # generate meshes
            print("########## Meshes generation") ; sys.stdout.flush()
            
            pool = Pool(15)
            pool.map(adaptInternal_star, zip(range(1,options.nbrAdap+1), [options.dim]*options.nbrAdap))
            gc.collect()



