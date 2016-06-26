import numpy as np
from numpy import linalg as LA
import sys, os, time

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
    
    Tend = 0.1
    nbrPtfxIte = 3
    nbrAdap = 2
    
    dtAdap = Tend/nbrAdap
    
    mesh = UnitSquareMesh(50,50)
    V_tmp = FunctionSpace(mesh, 'DG', 0)
    
    entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
    entity_dofs[0] = mesh.geometric_dimension()
    coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
    
    
    
    for i in range(nbrPtfxIte) :
        
        print "\n\n\n##########  i: %d\n" % i
        
        j = 0
        if i == 0 :
            meshes = nbrAdap*[mesh]
        else : 
            meshes = newmeshes
            mesh = meshes[j]
            entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
            entity_dofs[0] = mesh.geometric_dimension()
            coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
            
        solIni = solIniAdvec(mesh)
        writeMesh(mesh, coordSection, "bubble.%d" % j)
        writeSol(mesh, solIni, coordSection, "bubble.%d" % j)
        
        hessianMetrics = []
        newmeshes = []
        
        for j in range(1, nbrAdap+1) :
            
            print "\n#####  j: %d\n" % j
                                        
            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            if i == 0:
                if j > 1 :
                    solIni = sol
            else :
                if j > 1 :
                    meshOld = mesh
                    mesh = meshes[j-1]
                    entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
                    entity_dofs[0] = mesh.geometric_dimension()
                    coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
                    
                    solIni = Function(FunctionSpace(mesh, 'CG', 1))
                    interpol(sol, meshOld, solIni, mesh)
            
            # solve 
            tIni = (j-1)*dtAdap
            tEnd = j*dtAdap
            sol, hesMet = solveAdvec(mesh, solIni, tIni, tEnd)
            hessianMetrics.append(hesMet)                
                        
            writeMesh(mesh, coordSection, "bubble.%d" % j)
            writeSol(mesh, sol, coordSection, "bubble.%d" % j)
        
        ######## End of the loop over sub-intervals
        
        if i < (nbrPtfxIte-1) :
        
            # normalizeMetric
            print "########## Metrics computation"
            metrics = normalizeMetrics(hessianMetrics, meshes)
        
            # generate meshes
            print "########## Meshes generation"
            for j in range(nbrAdap) :
                print "##### Adap procedure started %d" %(j+1)
                newmesh = adapt(meshes[j],metrics[j])
                print "##### Adap procedure completed %d" %(j+1)
                newmeshes.append(newmesh)
                FS_tmp = FunctionSpace(newmesh, 'CG', 1)
                
                entity_dofs = np.zeros(newmesh._topological_dimension+1, dtype=np.int32)
                entity_dofs[0] = newmesh.geometric_dimension()
                coordSection = newmesh._plex.createSection([1], entity_dofs, perm=newmesh.topology._plex_renumbering)
                writeMesh(newmesh, coordSection, "newmesh.%d" % (j+1))
        
        
        
        
if __name__ == '__main__':

    parameters["pyop2_options"]["log_level"] = "WARNING"

    main()