import numpy as np
from numpy import linalg as LA
import sys, os, time
from eigVal import *
from metric import *
from inout import *

t0 = time.clock()
from firedrake import *
tloadfdk = time.clock()
print "## module firedrake loaded (%ds)" % (tloadfdk - t0); sys.stdout.flush()

INF = float("inf")






def main() :
    
    Tend = 6
    Nite = 3
    Nadap = 20
    
    dtAdap = Tend/Nadap
    
    
    for i in range(Nite):
        
        for j in range(Nadap):
            
            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            interpol(solPrev, meshPrev, solIni, meshIni)
            
            # solve 
            tIni = j*dtAdap
            tend = (j+1)*dtAdap
            sol, hessianMetric = solve(mesh, solIni, tIni, tEnd)
                        
            
        
        # normalizeMetric
        
        
        # generate meshes