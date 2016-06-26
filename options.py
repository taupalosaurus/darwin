class Options :
    
    
    def __init__(self) :
        
        # adaptation parameters    
        self.nbrPtfxIte = 1
        self.nbrAdap = 1
        self.nbrSpl = 0 
        
        # metrics computation parameters
        self.p = 2
        self.N = 1000
        self.a = 1000
        self.hmin = 0.005
        self.hmax = 0.3
        
        # solver parameters
        self.T = 6
        self.Tend = 1
        self.n = 50
        self.tSav = 0
        
                
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
        