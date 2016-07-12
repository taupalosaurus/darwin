from firedrake import *
#import numpy as np
#from numpy import linalg as LA
import sys, os, time

from options import *
from mesh import *
from inout import *
from solve import *
from adapt import *
from interpol import *



def ptfx (options) :
    
    if options.dim == 2:
        mesh = Meshd(UnitSquareMesh(options.n, options.n))
    else:
        mesh = Meshd(UnitCubeMesh(options.n, options.n, options.n))
    
    dtAdap = float(options.Tend)/options.nbrAdap

    
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
        writeGmf(mesh.mesh, 1, "boundary_ids", "bubble.%d" % j, solIni, 1, "bubble.%d" % j, mesh.section)
        
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
            tIni, tEnd  = (j-1)*dtAdap, j*dtAdap
            sol, hesMet, t = solveAdvec(mesh, solIni, tIni, tEnd, options)
            hessianMetrics.append(hesMet)
            writeGmf(mesh.mesh, 1, "boundary_ids", "bubble.%d" % j, sol, 1, "bubble.%d" % j, mesh.section)
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



        ######## End of the loop over sub-intervals
        
        if i < (options.nbrPtfxIte-1) :
        
            # normalizeMetric
            print "########## Metrics computation" ; sys.stdout.flush()
            metrics = normalizeUnsteadyMetrics(hessianMetrics, meshes, options)

        
            # generate meshes
            print "########## Meshes generation" ; sys.stdout.flush()
            for j in range(options.nbrAdap) :
                print "##### Adap procedure started %d" %(j+1) ; sys.stdout.flush()
                newmesh = adaptInternal(meshes[j], metrics[j])
                print "##### Adap procedure completed %d" %(j+1) ; sys.stdout.flush()
                newmeshes.append(newmesh)
                writeGmf(newmesh.mesh, 1, "boundary_ids", "newmesh.%d" % (j+1), None, None, None, newmesh.section)
