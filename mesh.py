from firedrake import *
import numpy as np


class Meshd() :
    
    
    def __init__(self, mesh, reorderPlex=True, computeAltMin=True) :
        
        self.mesh = mesh 

        self.V = FunctionSpace(self.mesh, 'CG', 1)
        
        self.altMin = Function(self.V)
        
        entity_dofs = np.zeros(self.mesh._topological_dimension+1, dtype=np.int32)
        entity_dofs[0] = self.mesh.geometric_dimension()
        self.section = self.mesh._plex.createSection([1], entity_dofs, perm=self.mesh.topology._plex_renumbering)

        if reorderPlex :
            with self.mesh.coordinates.dat.vec_ro as coords:
                self.mesh.topology._plex.setCoordinatesLocal(coords)

        if computeAltMin:
            if self.mesh._topological_dimension == 2:
                self.altMin.interpolate(2*CellVolume(self.mesh)/MaxCellEdgeLength(self.mesh))
            else :
                self.altMin.interpolate(Circumradius(self.mesh)*MinCellEdgeLength(self.mesh)/MaxCellEdgeLength(self.mesh))
                print "#### WARNING minimal altitude computed very approximately in 3D"
#                self.altMin.interpolate(6*CellVolume(self.mesh)/MaxCellFacetArea(self.mesh))
        
        
        