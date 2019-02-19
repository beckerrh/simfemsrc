from __future__ import print_function
import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("@libpythonpath@")

import mesh.geometry
import tools.plot
import numpy as np
import simfempy
import simfemcdr
import time
import math

#---------------------------------------------------------#
def main():
    if len(sys.argv)>1:
        if sys.argv[1] == "1D":
            meshtypes = ["LineMesh"]
        elif sys.argv[1] == "2D":
            meshtypes = ["TriangleMesh"]
        elif  sys.argv[1] == "3D":
            meshtypes = ["TetrahedralMesh"]
        else:
            raise KeyError("unknwon argument", sys.argv[1])
    meshtypes = ["LineMesh", "TriangleMesh", "TetrahedralMesh"]
    meshtypes = ["TriangleMesh"]
    testname = "poisson"
    # testname = "cdrexp"
    # testname = "rdcosh"
    fem = "P1"
    # fem = "CR1"
    niter = 3
    methods=["traditionalintegration", "nitscheintegration", "newnitscheintegration","traditional", "nitsche", "newnitsche"]
    methods=["traditionalintegration", "nitscheintegration", "newnitscheintegration"]
    # methods=["traditionalintegration"]
    methods=["traditional", "nitsche", "newnitsche"]
    # methods=["traditional", "traditionalintegration"]
    # methods=["nitsche", "nitscheintegration"]
    # methods=["newnitsche", "newnitscheintegration"]
    for meshtype in meshtypes:
        test(meshtype, testname, methods, niter, fem)


#---------------------------------------------------------#
def test(meshtype, testname, methods, niter, fem):
    if  meshtype == "LineMesh":
        from mesh.geometries.unitline import GeometryClass
    elif  meshtype == "TriangleMesh":
        from mesh.geometries.unitsquare import GeometryClass
    elif  meshtype == "TetrahedralMesh":
        from mesh.geometries.unitcube import GeometryClass

    if testname == "cdrexp":
        diff = 0.001
        beta = "east"
        alpha = 0.0
        application="cdexplayer"
    elif testname == "rdcosh":
        diff = 0.0001
        beta = "zero"
        alpha = 1.0
        application="rdcosh"
    elif testname == "poisson":
        diff = 1.0
        beta = "zero"
        alpha = 0.0
        application="quadratic"
        # application="cosinus"
        # application="linear"
    errL1 = {}
    errL2 = {}
    errH1 = {}
    times = {}
    for method in methods:
        times[method] = 0.0
        errL1[method] =  np.zeros(niter)
        errL2[method] =  np.zeros(niter)
        errH1[method] =  np.zeros(niter)
    geom=GeometryClass(hmean=0.5)
    name = geom.runGmsh(outdir="Data", number=0)
    hmeans = np.zeros(niter)
    for ih in range(niter):
        # print("---- hmean=", hmean)
        # geom2=GeometryClass(hmean=hmean)
        # geom2.runGmsh(outdir="Data")
        name = geom.runGmshRefine(number=ih, outdir="Data")

        partion_id = 1
        construct_bdrymeshes = False
        mesh = simfempy.create(meshtype, partion_id, construct_bdrymeshes)
        meshfilename = "Data/"+name+'.msh'
        print("meshfilename",meshfilename)
        mesh.readGmsh(meshfilename)
        mesh.addGeometryObject('MeasureOfCell')
        mesh.addGeometryObject('Normals')
        # mesh.save("Data/"+geom.name+'.h5')

        hmeans[ih] = math.pow(mesh.getNCells(), -1/mesh.getDimension())


        for im,method in enumerate(methods):
            print("####### %s #####" %(method))
            start = time.time()
            solver = simfemcdr.Solver()
            solver.setParameter("method",method)
            solver.setParameter("application",application)
            solver.setParameter("beta",beta)
            solver.setParameter("fem",fem)
            solver.setParameter("diff", diff);
            solver.setParameter("alpha", alpha);
            solver.setParameter("deltasupg", 0.5);
            solver.setMesh(mesh)
            # solver.loadMesh(meshtype, "Data/"+geom.name+'.h5')
            solver.init()
            info = solver.getInfo()
            # print(info)
            d = solver.run()
            errL1[method][ih] = d["L1"]
            errL2[method][ih] = d["L2"]
            errH1[method][ih] = d["H1"]
            # print("d=",d)
            solver.writeXdmf()
        times[method] = time.time()-start
    for method in methods:
        print("errL2 %-20s" %(method), errL2[method])
    print()
    for method in methods:
        print("errL1 %-20s" %(method), errL1[method])
    print()
    for method in methods:
        print("errH1 %-20s" %(method), errH1[method])
    # for method in methods:
    #     errH1[method] = np.sqrt( diff*errH1[method]**2 + errL2[method]**2)
    for method in methods: print("%-20s %10g" %(method,times[method]))

    datatoplot = {}
    datatoplot['L1'] = errL1
    datatoplot['L2'] = errL2
    datatoplot['H1'] = errH1
    methodsnames={}
    methodsnames["traditionalintegration"] = "trad"
    methodsnames["nitscheintegration"] = "nit"
    methodsnames["newnitscheintegration"] = "new"
    methodsnames["traditional"] = "trad"
    methodsnames["nitsche"] = "nit"
    methodsnames["newnitsche"] = "new"

    print("hmeans",hmeans)
    simfemplot = tools.plot.SimFemPlot(methods, params=hmeans, param='hmean')
    simfemplot.ploterrors(datatoplot)


#---------------------------------------------------------#
if __name__ == '__main__':
    main()
