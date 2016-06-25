import numpy as np
from numpy import linalg as LA
import sys, os, time

from solve import *
from inout import *

t0 = time.clock()
from firedrake import *
tloadfdk = time.clock()
print "## module firedrake loaded (%ds)" % (tloadfdk - t0); sys.stdout.flush()

INF = float("inf")






def main() :
    
    Tend = 1.
    Nite = 1
    Nadap = 20
    
    dtAdap = Tend/Nadap
    
    mesh = UnitSquareMesh(50,50)
    V_tmp = FunctionSpace(mesh, 'DG', 0)
    
    entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
    entity_dofs[0] = mesh.geometric_dimension()
    coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
    
    
    
    for i in range(Nite):
        
        print "DEBUG  i: %d" % i
        
        solIni = solIniAdvec(mesh)
        j = 0
        writeMesh(mesh, coordSection, "bubble.%d" % j)
        writeSol(mesh, solIni, coordSection, "bubble.%d" % j)
        
        for j in range(1,Nadap):
            
            print "DEBUG  j: %d" % j
            
            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            #interpol(solPrev, meshPrev, solIni, meshIni)
            if j > 1 :
                solIni = sol
            
            # solve 
            tIni = j*dtAdap
            tEnd = (j+1)*dtAdap
            sol, hessianMetric = solveAdvec(mesh, solIni, tIni, tEnd)
                        
            writeMesh(mesh, coordSection, "bubble.%d" % j)
            writeSol(mesh, sol, coordSection, "bubble.%d" % j)
        
        # normalizeMetric
        
        
        # generate meshes
        
        
        
        
        
        
        
        
if __name__ == '__main__':

    parameters["pyop2_options"]["log_level"] = "WARNING"

    main()