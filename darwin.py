import numpy as np
from numpy import linalg as LA
import sys, os, time, datetime

from options import *
from mesh import *
from algo_ptfx import *
from solve import *
from adapt import *
from interpol import *

t0 = time.clock()
from firedrake import *
tloadfdk = time.clock()
#print("## module firedrake loaded (%ds)" % (tloadfdk - t0)); sys.stdout.flush()

INF = float("inf")






def main() :


    today = datetime.datetime.today()
    print("#############################")
    print("#           Darwin          #")
    print("#  Imperial College London  #")
    print("#    %4d-%02d-%02d %02d:%02d:%02d    #" %(today.year, today.month, today.day, today.hour, today.minute, today.second) )
    print("#############################\n")

    

    options = Options(algo=1,
                      dim = 2,
                      nbrPtfxIte=2,
                      nbrAdap = 2,
                      nbrSpl = 5,
                      p = 2,
                      N = 2000,
                      hmin = 100,
                      hmax = 1000,
                      a = 10,
                      Tend = 400,
                      cfl = 20,
                      n = 5,
                      T = 1.,
                      nbrSav = 5,
                      snes_rtol = 1e-2,
                      ksp_rtol = 1e-5)


#    options = Options(dim = 2)
#    options.setSmallTest3d()
    options.printOptions()
    sys.stdout.flush()


    ptfx(options)



        
        
if __name__ == '__main__':


    parameters["pyop2_options"]["log_level"] = "WARNING"
  #  parameters["assembly_cache"]["enabled"] = False

    tini = time.clock()

    main()

    tend = time.clock()
    ttotal = time.gmtime(tend-tini)
    print("#####  TOTAL TIME : %ddays %dhours %dmin %dsec" %(ttotal[2]-1, ttotal[3], ttotal[4], ttotal[5]))
