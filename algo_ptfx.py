

def ptfx (options) :
    
    mesh = Meshd(UnitSquareMesh(options.n, options.n))
    
    dtAdap = float(options.Tend)/options.nbrAdap

    
    for i in range(options.nbrPtfxIte) :
        
        print "\n\n\n##########  i: %d\n" % i ; sys.stdout.flush()
        
        j = 0
        if i == 0 :
            meshes = options.nbrAdap*[mesh]
        else : 
            meshes = newmeshes
            mesh = meshes[j]
                    
        ####  Initial solution    
        solIni = solIniAdvec(mesh)
        writeMesh(mesh, "bubble.%d" % j)
        writeSol(mesh, solIni, "bubble.%d" % j)
        
        hessianMetrics = []
        newmeshes = []
        
        for j in range(1, options.nbrAdap+1) :
            
            print "\n#####  j: %d\n" % j ; sys.stdout.flush()
                                        
            # interpolate previous solution on new mesh if necessary (ie i > 0 and j > 0)
            if i == 0:
                if j > 1 :
                    solIni = sol
            else :
                if j > 1 :
                    meshOld = mesh
                    mesh = meshes[j-1]
                    solIni = Function(FunctionSpace(mesh.mesh, 'CG', 1))
                    interpol(sol, meshOld, solIni, mesh)
            
            # solve 
            tIni, tEnd  = (j-1)*dtAdap, j*dtAdap
            sol, hesMet, t = solveAdvec(mesh, solIni, tIni, tEnd, options)
            hessianMetrics.append(hesMet)                
            writeMesh(mesh, "bubble.%d" % j)
            writeSol(mesh, sol, "bubble.%d" % j)
            if options.nbrSav > 0 :
                kIni = (j-1)*options.nbrSav
                for k in range(options.nbrSav+1) :
                    kGlob = kIni + k
                    print "DEBUG   film.%d -> film.%d" %(k, kGlob)
                    if k == 0 :
                        os.rename("film_tmp.%d.mesh" % k, "film.%d.mesh" % kGlob)
                    else :
                        if os.path.exists("film.%d.mesh" % kGlob) : os.remove("film.%d.mesh" % kGlob)
                        os.symlink("film.%d.mesh" % kIni, "film.%d.mesh" % kGlob)
                    os.rename("film_tmp.%d.sol" % k, "film.%d.sol" % kGlob)



        ######## End of the loop over sub-intervals
        
        if i < (options.nbrPtfxIte-1) :
        
            # normalizeMetric
            print "########## Metrics computation" ; sys.stdout.flush()
            metrics = normalizeUnsteadyMetrics(hessianMetrics, meshes, options)
        
            # generate meshes
            print "########## Meshes generation" ; sys.stdout.flush()
            for j in range(options.nbrAdap) :
                print "##### Adap procedure started %d" %(j+1) ; sys.stdout.flush()
                newmesh = adaptInternal(meshes[j], metrics[j])
                print "##### Adap procedure completed %d" %(j+1) ; sys.stdout.flush()
                newmeshes.append(newmesh)
                writeMesh(newmesh, "newmesh.%d" % (j+1))