import pyop2 as op2
from hessian import *

class Options :
    
    
    def __init__(self, 
                 algo = 1,
                 nbrPtfxIte=1,
                 nbrAdap = 1,
                 nbrSpl = 0,
                 adaptStepFreq = 0,
                 p = 2,
                 N = 1000,
                 a = 1000,
                 hmin = 0.005,
                 hmax = 0.3,
                 steadyMetric = 1,
                 T = 6,
                 Tend = 1,
                 n = 50,
                 cfl = 0.95,
                 nbrSav = 0,
                 nbrSavTot = 0,
                 dtSav = 0) :
    
        # adaptation parameters
        self.algo = algo  # 1 for ptfx, 2 for frequent remeshes    
        self.nbrPtfxIte = nbrPtfxIte
        self.nbrAdap = nbrAdap
        self.nbrSpl = nbrSpl  # including initial and final samples, so this is the actual number
        self.adaptStepFreq = adaptStepFreq
        
        # metrics computation parameters
        self.p = p
        self.N = N
        self.a = a
        self.hmin = hmin
        self.hmax = hmax
        self.steadyMetric = steadyMetric # if algo = 2, 1 for pain 2001, 2 for fluidity, 3 for loseille 2011
        
        # solver parameters
        self.T = T
        self.Tend = Tend
        self.n = n
        self.cfl = cfl
        self.nbrSav = nbrSav  # number of sol saved per sub-interval in algo 1 initial save does not count (so I actually save nbrSav+1)
        self.nbrSavTot = nbrSavTot # total number of solutions saved for algo 2
        self.dtSav = dtSav
        self.nbrGlobSav  = 0 # keep track of the number of saved solutions
        
        if nbrSavTot > 0:
            self.dtSav = self.Tend/self.nbrSavTot
            self.nbrSav = 0 



        op2.init()
        self.absValHessian_kernel = op2.Kernel("""
void absValHessian(double * hess, double *lbdMin) {
  
  double lmin = *lbdMin;
  
%s
%s
%s
%s

}
""" % (computeEigVal_str, absValueHessian_str, truncLowHessian_str, rebuildHessian_str), "absValHessian")



#        usa2 = 1./(self.a*self.a)
#        ushmin2 = 1./(self.hmin*self.hmin)
#        ushmax2 = 1./(self.hmax*self.hmax)
#        op2.Const(1, ushmax2, dtype=float, name="lbdMin");
#        op2.Const(1, ushmin2, dtype=float, name="lbdMax");
#        op2.Const(1, usa2, dtype=float, name="usa2");
        self.absTruncMetric_kernel = op2.Kernel("""
void absTruncMetric_kernel(double * hess, double *lbdmin, double *lbdmax, double *usa2_p) {
    
    double lmin = *lbdmin;
    double lmax = *lbdmax;
    double usa2 = *usa2_p;
%s
%s
%s
%s
%s

}
""" % (computeEigVal_str, truncLowHessian_str, truncHighHessian_str, truncRatioHessian_str, rebuildHessian_str), "absTruncMetric_kernel")

                
    def setSmallTest(self):
        
        # adaptation parameters    
        self.nbrPtfxIte = 3
        self.nbrAdap = 2
        self.nbrSpl = 5 
        
        # metrics computation parameters
        self.p = 2
        self.N = 1000
        self.a = 1000
        self.hmin = 0.001
        self.hmax = 0.3
        
        # solver parameters
        self.T = 6
        self.Tend = 0.1
        self.n = 50
        self.cfl = 0.95
        self.nbrSav = 3
        