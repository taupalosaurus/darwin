# inspire de https://github.com/FEniCS/dolfin/blob/master/demo/undocumented/advection-diffusion/python/demo_advection-diffusion.py

import numpy as np
from numpy import linalg as LA
import sys, os, time
from inout import *

t0 = time.clock()
from firedrake import *
tloadfdk = time.clock()
print "## module firedrake loaded (%ds)" % (tloadfdk - t0); sys.stdout.flush()

INF = float("inf")


def main() :

    tini = time.clock()
    
    n = 100
    T = 6.
    Tend = T
    tsav = 0.01

    noAdap = 0    
    tadap = 0.01
    iteAdp = 3
    
    mesh = UnitSquareMesh(n,n)
    tloadmesh = time.clock()
    print "## mesh created (%ds, tot: %dm)  " % (tloadmesh - tini, (tloadmesh-tini)/60.); sys.stdout.flush()
    
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 2)
    Vn = FunctionSpace(mesh, "CG", 2)
    u, v = TrialFunction(Q), TestFunction(Q)
    u0 = Function(Q)     
    c, cn = Function(V), Function(Vn) # velocity and its norm
    
    # initial condition
    icExpr = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.75)*(x[1]-0.75) < 0.15*0.15)")
    u0.interpolate(icExpr)
        
    # boundary condition
    zero = Constant(0)
    bc = DirichletBC(Q, zero, 1)
        
    tsetup = time.clock()
    print "## setup finished (%ds, tot: %dm)  " % (tsetup - tloadmesh, (tsetup-tini)/60.); sys.stdout.flush()
    
    t = 0
    step = 0
    stepAdp = 0
    stepSav = 0
    
    print "\n\n###########  STEP %d  ###############" % step ; sys.stdout.flush()
    tinistep = time.clock() 

    if not noAdap :
        for i in range(iteAdp) :
            newmesh = adapt_int(mesh, u0, step)    
            mesh = newmesh    
            Q = FunctionSpace(mesh, "CG", 1)
            u0 = Function(Q)
            icExpr = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.75)*(x[1]-0.75) < 0.15*0.15)")
            u0.interpolate(icExpr)
        V = VectorFunctionSpace(mesh, "CG", 2)
        Vn = FunctionSpace(mesh, "CG", 2)
        u, v = TrialFunction(Q), TestFunction(Q)
        c = Function(V)
        cn = Function(Vn)
        zero = Constant(0)
        bc = DirichletBC(Q, zero, 1)
            


    print "### saving sol %d" % stepSav
    entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
    entity_dofs[0] = mesh.geometric_dimension()
    coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
    writeMesh(mesh, coordSection, "bubble.%d" % stepSav)
    writeSol(mesh, u0, coordSection, "bubble.%d" % stepSav)
    
    tendstep = time.clock()
    print "## end step (%ds, tot: %dm)  " % (tendstep-tinistep, (tendstep-tini)/60.); sys.stdout.flush()

    
    while t < Tend:
    
        step += 1
        print "\n\n###########  STEP %d  ###############" % step ; sys.stdout.flush()
        tinistep = time.clock() 
    
        cx = "-sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*cos(2*pi*%f)" % (t/T)
        cy = "sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(2*pi*%f)" % (t/T)
        c.interpolate(Expression([cx, cy]))
    
        nrmc = "sqrt((%s)*(%s)+(%s)*(%s))" % (cx, cx, cy, cy)
        cn.interpolate(Expression(nrmc))
        cmax = cn.dat.data.max()
        hmin = 1./n #### TODO don't leave that !!!
        dt = 0.5*hmin/cmax
        doSav = 0
        if int(t/tsav) != int((t+dt)/tsav) :
            doSav = 1
            stepSav += 1
            dt = stepSav*tsav - t
        doAdp = 0
        if (not noAdap) and (int(t/tadap) != int((t+dt)/tadap)) :
            doAdp = 1
            stepAdp += 1
            dt = stepAdp*tadap - t
        print "### t: %1.3f -> %1.3f      dt: %1.3e" % (t, t+dt, dt); sys.stdout.flush()
        
        a = v*u*dx + 0.5*dt*(v*dot(c, grad(u))*dx)
        L = v*u0*dx - 0.5*dt*(v*dot(c, grad(u0))*dx)
        solve(a==L,u0)
        
        t += dt
        
        if doAdp :
            for i in range(iteAdp):
                tbegadp = time.clock()
                newmesh = adapt_int(mesh, u0, stepAdp)        
                tendadp = time.clock()
                print "## metric + adaptation (%ds (%dm))  " % (tendadp-tbegadp, (tendadp-tbegadp)/60.); sys.stdout.flush()
                Q = FunctionSpace(newmesh, "CG", 1)
                u0_new = Function(Q)
                tbegitp = time.clock()
                interpol(u0, mesh, u0_new, newmesh)
                tenditp = time.clock()
                print "## interpolation (%ds (%dm))  " % (tenditp-tbegitp, (tenditp-tbegitp)/60.); sys.stdout.flush()
                mesh = newmesh 
                u0 = u0_new   
            V = VectorFunctionSpace(newmesh, "CG", 2)
            Vn = FunctionSpace(newmesh, "CG", 2)
            u, v = TrialFunction(Q), TestFunction(Q)
            c, cn = Function(V), Function(Vn)
            zero = Constant(0)
            bc = DirichletBC(Q, zero, 1)
        
        
        if doSav :
            print "### saving sol %d" % stepSav
            entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
            entity_dofs[0] = mesh.geometric_dimension()
            coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
            writeMesh(mesh, coordSection, "bubble.%d" % stepSav)
            writeSol(mesh, u0, coordSection, "bubble.%d" % stepSav)
                
        tendstep = time.clock()
        print "## end step (%ds, tot: %dm)  " % (tendstep-tinistep, (tendstep-tini)/60.); sys.stdout.flush()
    
    
    
    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print "#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5])






def adapt_int(mesh, u, i):

    p = 2.
    N = 2200
    use = 100
    umin = 8e-3
    a = 1000
    usa2 = 1/(a*a)
    hmin = 0.005
    ushmin2 = 1/(hmin*hmin)
    hmax = 0.3
    ushmax2 = 1/(hmax*hmax)

    M = TensorFunctionSpace(mesh, "CG", 1)

    # computation of the hessian
    sigma = TestFunction(M)
    H = Function(M)
    n = FacetNormal(mesh)
    L = inner(sigma, H)*dx
    L += (inner(div(sigma), grad(u))*dx - (sigma[0, 1]*n[1]*u.dx(0) + sigma[1, 0]*n[0]*u.dx(1))*ds)
    H_prob = NonlinearVariationalProblem(L, H)
    H_solv = NonlinearVariationalSolver(H_prob)
    H_solv.solve()
    
    # computation of the metric as in Pain 2001
    metric = Function(M)
    metOpt = 'Loseille2011'
    if metOpt == 'Pain2001':
      for iVer in range(mesh.topology.num_vertices()):
          hesLoc = H.dat.data[iVer] * 1/max(abs(u.sub(0).dat.data[iVer]), umin) * use
          meanDiag = 0.5*(hesLoc[0][1] + hesLoc[1][0]);
          hesLoc[0][1] = meanDiag
          hesLoc[1][0] = meanDiag
          lbd, v = LA.eig(hesLoc)
          # truncation of eigenvalues
          lbd1 = min(ushmin2, max(ushmax2,abs(lbd[0])))
          lbd2 = min(ushmin2, max(ushmax2,abs(lbd[1])))
          maxLbd = max(lbd1, lbd2)
          lbd1 = max(lbd1, usa2*maxLbd)
          lbd2 = max(lbd2, usa2*maxLbd)
          v1, v2 = v[0], v[1]
          # reconstruction of |Hu|
          metric.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
          metric.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
          metric.dat.data[iVer][1,0] = metric.dat.data[iVer][0,1]
          metric.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
    
    elif metOpt == 'Loseille2011' :
        # computation of the metric of Loseille 2011
        V = FunctionSpace(mesh, 'CG', 1)    
        detHes = Function(V)
        for iVer in range(mesh.topology.num_vertices()):
            hesLoc = H.dat.data[iVer]
            meanDiag = 0.5*(hesLoc[0][1] + hesLoc[1][0]);
            hesLoc[0][1] = meanDiag
            hesLoc[1][0] = meanDiag
            lbd, v = LA.eig(hesLoc)
            lbd1 = max(abs(lbd[0]), 1e-10)
            lbd2 = max(abs(lbd[1]), 1e-10)
            det = lbd1*lbd2
            v1, v2 = v[0], v[1]
            metric.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
            metric.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
            metric.dat.data[iVer][1,0] = metric.dat.data[iVer][0,1]
            metric.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
            metric.dat.data[iVer] *= pow(det, -1./(2*p+2))
            detHes.dat.data[iVer] = pow(det, p/(2.*p+2))
    
        detHesInt = assemble(detHes*dx)
        glbNrmCof = (N/detHesInt)
        metric *= glbNrmCof
    
        for iVer in range(mesh.topology.num_vertices()):
            lbd, v = LA.eig(metric.dat.data[iVer])
            # truncation of eigenvalues
            lbd1 = min(ushmin2, max(ushmax2,lbd[0]))
            lbd2 = min(ushmin2, max(ushmax2,lbd[1]))
            maxLbd = max(lbd1, lbd2)
            lbd1 = max(lbd1, usa2*maxLbd)
            lbd2 = max(lbd2, usa2*maxLbd)
            v1, v2 = v[0], v[1]
            # reconstruction of |Hu|
            metric.dat.data[iVer][0,0] = lbd1*v1[0]*v1[0] + lbd2*v2[0]*v2[0];
            metric.dat.data[iVer][0,1] = lbd1*v1[0]*v1[1] + lbd2*v2[0]*v2[1];
            metric.dat.data[iVer][1,0] = metric.dat.data[iVer][0,1]
            metric.dat.data[iVer][1,1] = lbd1*v1[1]*v1[1] + lbd2*v2[1]*v2[1];
    else :
        print "## ERROR : bad metop"
        exit(1)


    print "####  calling adapt procedure"
    mesh = adapt(mesh,metric)
    print "####  end call adapt procedure"

    return mesh


def interpol(u, mesh, unew, meshnew):
    
    plexnew = meshnew._plex
    vStart, vEnd = plexnew.getDepthStratum(0)
    
    entity_dofs = np.zeros(meshnew._topological_dimension+1, dtype=np.int32)
    entity_dofs[0] = meshnew.geometric_dimension()
    coordSectionnew = meshnew._plex.createSection([1], entity_dofs, perm=meshnew.topology._plex_renumbering)
        
    notInDomain = []
    for v in range(vStart, vEnd):
        offnew = coordSectionnew.getOffset(v)/2
        newCrd = meshnew.coordinates.dat.data[offnew]
        try :
            val = u.at(newCrd)
        except PointNotInDomainError :
            print "####  New vertex not in domain: %f %f" % (newCrd[0], newCrd[1])
            val = 0.
            notInDomain.append([v,INF,-1])
        finally :
            unew.dat.data[offnew] = val
    print "####  Number of points not in domain: %d / %d" % (len(notInDomain), meshnew.topology.num_vertices())
    
    if len(notInDomain) > 0 :
        plex = mesh._plex
        fStart, fEnd = plex.getHeightStratum(1)  # edges/facets
        vStart, vEnd = plex.getDepthStratum(0)
        
        entity_dofs = np.zeros(mesh._topological_dimension+1, dtype=np.int32)
        entity_dofs[0] = mesh.geometric_dimension()
        coordSection = mesh._plex.createSection([1], entity_dofs, perm=mesh.topology._plex_renumbering)
                
        for f in range(fStart, fEnd):
            if plex.getLabelValue("boundary_ids", f) == -1 : continue
            closure = plex.getTransitiveClosure(f)[0]
            crdE = [] # coordinates of the two vertices of the edge
            for cl in closure:
                if cl >= vStart and cl < vEnd : 
                    off = coordSection.getOffset(cl)/2
                    crdE.append(mesh.coordinates.dat.data[off])
            if len(crdE) != 2 : exit(16)
            vn =  [crdE[0][1]-crdE[1][1], crdE[0][0]-crdE[1][0]]# normal vector of the edge
            nrm = sqrt(vn[0]*vn[0]+vn[1]*vn[1])
            vn = [vn[0]/nrm, vn[1]/nrm]
            for nid in notInDomain:
                v = nid[0]
                offnew = coordSectionnew.getOffset(v)/2    
                crdP = meshnew.coordinates.dat.data[offnew]
                dst = abs(vn[0] * (crdE[0][0] - crdP[0]) + vn[1] * (crdE[0][1] - crdP[1]))
                if dst < nid[1]:
                    # control if the vertex is between the two edge vertices :  big assumption that corners are preserved
                    # e1P.e1e2
                    sca1 = (crdP[0] - crdE[0][0])*(crdE[1][0] - crdE[0][0]) + (crdP[1] - crdE[0][1])*(crdE[1][1] - crdE[0][1])
                    # e2P.e2e1
                    sca2 = (crdP[0] - crdE[1][0])*(crdE[0][0] - crdE[1][0]) + (crdP[1] - crdE[1][1])*(crdE[0][1] - crdE[1][1])
                    if sca1 >= 0 and sca2 > 0:
                        nid[1] = dst
                        nid[2] = f
                        
        for nid in notInDomain :
            v = nid[0]
            offnew = coordSectionnew.getOffset(v)/2    
            crdP = meshnew.coordinates.dat.data[offnew]
            val = -1
            if nid[1] > 0.01:
                cStart, cEnd = plex.getHeightStratum(0)
                barCrd = [0,0,0]
                inCell = 0
                for c in range(cStart, cEnd):
                    closure = plex.getTransitiveClosure(c)[0]
                    crdC = [] # coordinates of the three vertices of the triangle
                    val = [] # value of the function at the vertices of the traingle
                    for cl in closure:
                        if cl >= vStart and cl < vEnd :
                            off = coordSection.getOffset(cl)/2
                            crdC.append(mesh.coordinates.dat.data[off])
                            val.append(u.dat.data[off])
                    # barycentric coordinates of v
                    barCrd[0] = barCoord(crdP, crdC, 0)
                    barCrd[1] = barCoord(crdP, crdC, 1)
                    barCrd[2] = barCoord(crdP, crdC, 2)
                    if barCrd[0] >= 0 and barCrd[1] >= 0 and barCrd[2] >= 0 :
                        print "DEBUG  Cell : %1.4f %1.4f   %1.4f %1.4f   %1.4f %1.4f   bary:  %e %e %e" % (crdC[0][0], crdC[0][1], crdC[1][0], crdC[1][1], crdC[2][0], crdC[2][1], barCrd[0], barCrd[1], barCrd[2])
                        val = barCrd[0]*val[0] + barCrd[1]*val[1] + barCrd[2]*val[2]
                        inCell = 1
                        break
                if not inCell :
                    print "ERROR  vertex too far from the boundary but no enclosing cell found. Crd: %f %f" % (crdP[0], crdP[1])
                    exit(16)
            else :
                f = nid[2]
                if (f < fStart or f > fEnd):
                    print "## ERROR   f: %d,   fStart: %d,  fEnd: %d" % (f, fStart, fEnd)
                    exit(14)
                closure = plex.getTransitiveClosure(f)[0]
                crdE = [] # coordinates of the two vertices of the edge
                val = [] # value of the function at the vertices of the edge
                for cl in closure:
                    if cl >= vStart and cl < vEnd : 
                        off = coordSection.getOffset(cl)/2
                        crdE.append(mesh.coordinates.dat.data[off])
                        val.append(u.dat.data[off])
                if len(crdE) != 2 : 
                    print "## ERROR  number of points in crdE: %d" % len(crdE)
                    exit(16)
                edg =  [crdE[0][0]-crdE[1][0], crdE[0][1]-crdE[1][1]]   # normed vector e2e1
                nrm = sqrt(edg[0]*edg[0] + edg[1]*edg[1])
                edg = [edg[0]/nrm, edg[1]/nrm]
                # H = alpha e1 + (1-alpha) e2    and   alpha = e2P.e2e1/||e2e1||
                alpha = (crdP[0] - crdE[1][0])*edg[0] + (crdP[1] - crdE[1][1])*edg[1]
                if alpha > 1 : exit(23)
                val = alpha*val[0] + (1-alpha)*val[1]
            unew.dat.data[offnew] = val    
            
            
def barCoord(crdM, crdTri, i):   
    # computes the barycentric coordinate of M in triangle Tri = P0,P1,P2 w/r the ith vertex
    # crd = det(MPj,MPk)/det(PiPj,PiPk)
    j = (i+1) % 3
    k = (i+2) % 3
    res1 = (crdTri[j][0]-crdM[0])*(crdTri[k][1]-crdM[1]) - (crdTri[k][0]-crdM[0])*(crdTri[j][1]-crdM[1])
    res2 = (crdTri[j][0]-crdTri[i][0])*(crdTri[k][1]-crdTri[i][1]) - (crdTri[k][0]-crdTri[i][0])*(crdTri[j][1]-crdTri[i][1])
    res = res1 / res2
    return res

def isApprox(x, y):
    if x < y*1.001 and x > y*0.999: return 1
    else : return 1 
    
if __name__ == '__main__':
    
    parameters["pyop2_options"]["log_level"] = "WARNING"
    
    main()
    
    
    


# compute altMin
#lenMax = vnx[0]*vnx[0] + vny[0]*vny[0];
#lenMax = DMAX(lenMax, vnx[1]*vnx[1] + vny[1]*vny[1]);
#lenMax = DMAX(lenMax, vnx[2]*vnx[2] + vny[2]*vny[2]);
#altMin = twoAir*twoAir / lenMax;
#if ( minAlt > altMin ) {
#  minAlt = altMin;
#  iAlt   = iTri;
#}