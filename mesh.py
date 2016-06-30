from firedrake import *
import numpy as np


class Meshd() :
    
    
    def __init__(self, mesh) :
        
        self.mesh = mesh 
        
        self.altMin = Function(FunctionSpace(self.mesh, 'CG', 1))
        
        entity_dofs = np.zeros(self.mesh._topological_dimension+1, dtype=np.int32)
        entity_dofs[0] = self.mesh.geometric_dimension()
        self.section = self.mesh._plex.createSection([1], entity_dofs, perm=self.mesh.topology._plex_renumbering)
        
        #self.computeVerMinAlt()
        self.altMin.interpolate(2*CellVolume(self.mesh)/MaxCellEdgeLength(self.mesh))

        
    def computeVerMinAlt(self) :  # altMin is a Function(FunctionSpace(mesh, 'CG', 1))  
    
        plex = self.mesh._plex
        vStart, vEnd = plex.getDepthStratum(0)
        cStart, cEnd = plex.getHeightStratum(0)
    
        for iVer in range(0, vEnd-vStart) :
            self.altMin.dat.data[iVer] = 1e10
    
        for c in range(cStart, cEnd) :
    
            ver = []
            ne = []
            closure = plex.getTransitiveClosure(c)[0]
            for cl in closure:
                if cl >= vStart and cl < vEnd : 
                    off = self.section.getOffset(cl)/2
                    ver.append(self.mesh.coordinates.dat.data[off])
            # normals to edges
            ne.append([ver[0][1]-ver[1][1], ver[1][0]-ver[0][0]])
            ne.append([ver[1][1]-ver[2][1], ver[2][0]-ver[1][0]])
            ne.append([ver[2][1]-ver[0][1], ver[0][0]-ver[2][0]])
            are = 0.5*abs(ne[0][0]*ne[1][1]-ne[0][1]*ne[1][0])
            edgLenMax = ne[0][0]*ne[0][0] + ne[0][1]*ne[0][1]
            edgLenMax = max(edgLenMax, ne[1][0]*ne[1][0] + ne[1][1]*ne[1][1])
            edgLenMax = max(edgLenMax, ne[2][0]*ne[2][0] + ne[2][1]*ne[2][1])
            altMinTri = 2*are/sqrt(edgLenMax)
            for cl in closure:
                if cl >= vStart and cl < vEnd : 
                    off = self.section.getOffset(cl)/2
                    self.altMin.dat.data[off] = min(self.altMin.dat.data[off], altMinTri)
    
    
        
        
        