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






def main() :
    
    
#    options = Options(nbrPtfxIte=4,
#                      nbrAdap = 60,
#                      nbrSpl = 20,
#                      p = 2,
#                      N = 60000,
#                      hmin = 0.0005,
#                      hmax = 0.3,
#                      T = 6,
#                      Tend = 3,
#                      n = 200,
#                      nbrSav = 4)
    options = Options(nbrPtfxIte=3,
                      nbrAdap = 40,
                      nbrSpl = 20,
                      p = 2,
                      N = 15000,
                      hmin = 0.001,
                      hmax = 0.3,
                      T = 6,
                      Tend = 3,
                      n = 75,
                      nbrSav = 0)
#    options.setSmallTest()
                
    mesh = Meshd(UnitSquareMesh(options.n, options.n))

    
    for i in range(options.nbrPtfxIte) :
        
        print "\n\n\n##########  i: %d\n" % i ; sys.stdout.flush()
        
        j = 0
        if i == 0 :
            meshes = options.nbrAdap*[mesh]
        else : 
            meshes = newmeshes
            mesh = meshes[j]
                    
        ####  Initial solution    
        solIni = solIniAdvec(mesh)
        writeMesh(mesh, "bubble.%d" % j)
        writeSol(mesh, solIni, "bubble.%d" % j)
        
        hessianMetrics = []
        newmeshes = []
        
        for j in range(1, options.nbrAdap+1) :
            
            print "\n#####  j: %d\n" % j ; sys.stdout.flush()
                                        
            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            if i == 0:
                if j > 1 :
                    solIni = sol
            else :
                if j > 1 :
                    meshOld = mesh
                    mesh = meshes[j-1]
                    solIni = Function(FunctionSpace(mesh.mesh, 'CG', 1))
                    interpol(sol, meshOld, solIni, mesh)
            
            # solve 
            dtAdap = float(options.Tend)/options.nbrAdap
            tIni, tEnd  = (j-1)*dtAdap, j*dtAdap
            sol, hesMet = solveAdvec(mesh, solIni, tIni, tEnd, options)
            hessianMetrics.append(hesMet)                
            writeMesh(mesh, "bubble.%d" % j)
            writeSol(mesh, sol, "bubble.%d" % j)
            if options.nbrSav > 0 :
                kIni = (j-1)*options.nbrSav
                for k in range(options.nbrSav+1) :
                    kGlob = kIni + k
                    print "DEBUG   film.%d -> film.%d" %(k, kGlob)
                    if k == 0 :
                        os.rename("film_tmp.%d.mesh" % k, "film.%d.mesh" % kGlob)
                    else :
                        if os.path.exists("film.%d.mesh" % kGlob) : os.remove("film.%d.mesh" % kGlob)
                        os.symlink("film.%d.mesh" % kIni, "film.%d.mesh" % kGlob)
                    os.rename("film_tmp.%d.sol" % k, "film.%d.sol" % kGlob)



        ######## End of the loop over sub-intervals
        
        if i < (options.nbrPtfxIte-1) :
        
            # normalizeMetric
            print "########## Metrics computation" ; sys.stdout.flush()
            metrics = normalizeMetrics(hessianMetrics, meshes, options)
        
            # generate meshes
            print "########## Meshes generation" ; sys.stdout.flush()
            for j in range(options.nbrAdap) :
                print "##### Adap procedure started %d" %(j+1) ; sys.stdout.flush()
                newmesh = adaptInternal(meshes[j], metrics[j])
                print "##### Adap procedure completed %d" %(j+1) ; sys.stdout.flush()
                newmeshes.append(newmesh)
                writeMesh(newmesh, "newmesh.%d" % (j+1))
        
        
        
if __name__ == '__main__':


    parameters["pyop2_options"]["log_level"] = "WARNING"

    tini = time.clock()

    main()

    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print "#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5])
