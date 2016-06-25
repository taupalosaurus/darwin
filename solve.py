from firedrake import *
import sys





def computeDtAdvec(mesh, cn, cxExpr, cyExpr):
    
    n = 50
    
    nrmc = "sqrt((%s)*(%s)+(%s)*(%s))" % (cxExpr, cxExpr, cyExpr, cyExpr)
    cn.interpolate(Expression(nrmc))
    cmax = cn.dat.data.max()
    hmin = 1./n #### TODO don't leave that !!!
    dt = 0.5*hmin/cmax

    return dt




def solveAdvec(mesh, solIni, tIni, tEnd):
    
    T = 6
    
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 2)
    Vn = FunctionSpace(mesh, "CG", 2)
    u, v = TrialFunction(Q), TestFunction(Q)
    u0 = Function(Q)     
    c, cn = Function(V), Function(Vn) # velocity and its norm
            
    zero = Constant(0)
    bc = DirichletBC(Q, zero, 1)
            
    step = 0
    print "####  step %d " % step ; sys.stdout.flush()
    
    t = tIni
    u0.assign(solIni)
     
    
    while t < tEnd:
    
        step += 1
        print "####  step %d " % step ; sys.stdout.flush()
    
        cxExpr = "-sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*cos(2*pi*%f)" % (t/T)
        cyExpr = "sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(2*pi*%f)" % (t/T)
        c.interpolate(Expression([cxExpr, cyExpr]))
    
        dt = computeDtAdvec(mesh, cn, cxExpr, cyExpr)
        if (t+dt > tEnd) : dt = tEnd - t
        
        a = v*u*dx + 0.5*dt*(v*dot(c, grad(u))*dx)
        L = v*u0*dx - 0.5*dt*(v*dot(c, grad(u0))*dx)
        solve(a==L,u0)
        
        t += dt
        
            
            
    return [u0, 0]
    
    
    
def solIniAdvec(mesh):
    
    Q = FunctionSpace(mesh, "CG", 1)
    icExpr = Expression("((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.75)*(x[1]-0.75) < 0.15*0.15)")
    u0 = Function(Q).interpolate(icExpr)
    
    return u0