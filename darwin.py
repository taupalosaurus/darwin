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
#print "## module firedrake loaded (%ds)" % (tloadfdk - t0); sys.stdout.flush()

INF = float("inf")






def main() :


    today = datetime.datetime.today()
    print "#############################"
    print "#           Darwin          #"
    print "#  Imperial College London  #"
    print "#    %4d-%02d-%02d %02d:%02d:%02d    #" %(today.year, today.month, today.day, today.hour, today.minute, today.second) 
    print "#############################\n"

    

#    options = Options(algo=1,
#                      dim = 3,
#                      nbrPtfxIte=5,
#                      nbrAdap = 100,
#                      nbrSpl = 15,
#                      p = 2,
#                      N = 150000,
#                      hmin = 0.0008,
#                      hmax = 0.5,
#                      a = 100,
#                      T = 6,
#                      Tend = 3,
#                      cfl = 2.,
#                      n = 80,
#                      nbrSav = 3,
#                      snes_rtol = 1e-2,
#                      ksp_rtol = 1e-5)
    options = Options(algo=1,
                      dim = 3,
                      nbrPtfxIte=5,
                      nbrAdap = 30,
                      nbrSpl = 15,
                      p = 2,
                      N = 300000,
                      hmin = 0.0008,
                      hmax = 0.5,
                      a = 100,
                      T = 6,
                      Tend = 1.5,
                      cfl = 5.,
                      n = 100,
                      nbrSav = 5,
                      snes_rtol = 1e-2,
                      ksp_rtol = 1e-5)
#    options = Options(algo=1,
#                      dim = 2,
#                      nbrPtfxIte=6,
#                      nbrAdap = 90,
#                      nbrSpl = 20,
#                      p = 2,
#                      N = 80000,
#                      hmin = 0.0007,
#                      hmax = 0.3,
#                      a = 25,
#                      T = 6,
#                      Tend = 6,
#                      cfl = 1.5,
#                      n = 300,
#                      nbrSav = 4,
#                      snes_rtol = 1e-2,
#                      ksp_rtol = 1e-5)
#    
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

#    options = Options(dim = 2)
    options.setSmallTest3d()
    options.printOptions()
    sys.stdout.flush()


    if options.algo == 1 :
        ptfx(options)
    elif options.algo == 2 :
        freqRemesh(options)
    else :
        print "####  ERROR  only two adaptation algrithms"


        
        
if __name__ == '__main__':


    parameters["pyop2_options"]["log_level"] = "WARNING"
  #  parameters["assembly_cache"]["enabled"] = False

    tini = time.clock()

    main()

    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print "#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5])
