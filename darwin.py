import numpy as np
from numpy import linalg as LA
import sys, os, time, datetime

from options import *
from mesh import *
from inout import *
from algo_ptfx import *
from algo_freqRemesh import *
from solve import *
from adapt import *
from interpol import *

t0 = time.clock()
from firedrake import *
tloadfdk = time.clock()
print "## module firedrake loaded (%ds)" % (tloadfdk - t0); sys.stdout.flush()

INF = float("inf")






def main() :


    today = datetime.datetime.today()
    print "#############################"
    print "#           Darwin          #"
    print "#  Imperial College London  #"
    print "#    %4d-%02d-%02d %02d:%02d:%02d    #" %(today.year, today.month, today.day, today.hour, today.minute, today.second) 
    print "#############################\n"

    

    options = Options(algo=1,
                      dim = 3,
                      nbrPtfxIte=3,
                      nbrAdap = 15,
                      nbrSpl = 15,
                      p = 2,
                      N = 50000,
                      hmin = 0.003,
                      hmax = 0.3,
                      a = 2,
                      T = 6,
                      Tend = 1.5,
                      n = 70,
                      nbrSav = 3,
                      snes_rtol = 1e-2,
                      ksp_rtol = 1e-5)
#    options = Options( algo = 2,
#                       adaptStepFreq = 3,
#                       p = 2,
#                       N = 100,
#                       a = 1000,
#                       hmin = 0.001,
#                       hmax = 0.3,
#                       steadyMetric = 1,
#                       T = 6,
#                       Tend = 1,
#                       n = 150,
#                       cfl = 0.95,
#                       nbrSavTot = 150)

#    options.setSmallTest3d()
    options.printOptions()


    if options.algo == 1 :
        ptfx(options)
    elif options.algo == 2 :
        freqRemesh(options)
    else :
        print "####  ERROR  only two adaptation algrithms"


        
        
if __name__ == '__main__':


    parameters["pyop2_options"]["log_level"] = "WARNING"

    tini = time.clock()

    main()

    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print "#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5])
