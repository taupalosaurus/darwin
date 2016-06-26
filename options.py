class Options :
    
    
    def setParam(self, 
                 nbrPtfx=1,
                 nbrAdap = 1,
                 nbrSpl = 0,
                 p = 2,
                 N = 1000,
                 a = 1000,
                 hmin = 0.005,
                 hmax = 0.3,
                 T = 6,
                 Tend = 1,
                 n = 50,
                 tSav = 0) :
    
        # adaptation parameters    
        self.nbrPtfxIte = nbrPtfxIte
        self.nbrAdap = nbrAdap
        self.nbrSpl = nbrSpl 
        
        # metrics computation parameters
        self.p = p
        self.N = N
        self.a = a
        self.hmin = hmin
        self.hmax = hmax
        
        # solver parameters
        self.T = T
        self.Tend = Tend
        self.n = n
        self.tSav = tSav
        
                
    def setSmallTest(self):
        
        # adaptation parameters    
        self.nbrPtfxIte = 3
        self.nbrAdap = 2
        self.nbrSpl = 5 
        
        # metrics computation parameters
        self.p = 2
        self.N = 10000
        self.a = 1000
        self.hmin = 0.005
        self.hmax = 0.3
        
        # solver parameters
        self.T = 6
        self.Tend = 0.1
        self.n = 50
        self.tSav = 0
        
    def setParam(self, 
                 nbrPtfx=1,
                 nbrAdap = 1,
                 nbrSpl = 0,
                 p = 2,
                 N = 1000,
                 a = 1000,
                 hmin = 0.005,
                 hmax = 0.3,
                 T = 6,
                 Tend = 1,
                 n = 50,
                 tSav = 0) :
    
        # adaptation parameters    
        self.nbrPtfxIte = nbrPtfxIte
        self.nbrAdap = nbrAdap
        self.nbrSpl = nbrSpl 
        
        # metrics computation parameters
        self.p = p
        self.N = N
        self.a = a
        self.hmin = hmin
        self.hmax = hmax
        
        # solver parameters
        self.T = T
        self.Tend = Tend
        self.n = n
        self.tSav = tSav