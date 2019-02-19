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

#---------------------------------------------------------#
def main():
    testname = "poisson"
    fem = "P1"
    hmeans=[0.5*0.5**i for i in range(3)]
    methods=["traditional"]

    from mesh.geometries.unitsquare import GeometryClass

    meshtype = "TriangleMesh"
    simfemplot = tools.plot.SimFemPlot()
    simfemplot = tools.plot.SimFemPlot(methods, params=hmeans, param='hmean')
    diff = 1.0
    beta = "zero"
    alpha = 1.0
    application="quadratic"
    errL1 = {}
    errL2 = {}
    errH1 = {}
    times = {}
    for method in methods: times[method] = 0.0
    geom=GeometryClass(hmean=1)
    name = geom.runGmsh(outdir="Data", number=0)
    for ih in range(3):
        print("---- ih=", ih)
        name = geom.runGmshRefine(number=0, outdir="Data")
        print("name", name)
        partion_id = 1
        construct_bdrymeshes = False
        mesh = simfempy.create(meshtype, partion_id, construct_bdrymeshes)
        mesh.readGmsh("Data/"+name+'.msh')
        mesh.addGeometryObject('MeasureOfCell')
        mesh.addGeometryObject('Normals')
        mesh.save("Data/"+name+'.h5')
        print("mesh", mesh)
        for im,method in enumerate(methods):
            print("####### %s #####" %(method))
            errL1[method] =  np.zeros(len(hmeans))
            errL2[method] =  np.zeros(len(hmeans))
            errH1[method] =  np.zeros(len(hmeans))
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
            # solver.loadMesh(meshtype, "Data/"+name+'.h5')
            solver.init()
            info = solver.getInfo()
            print(info)
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

    simfemplot.ploterrors(datatoplot)


#---------------------------------------------------------#
if __name__ == '__main__':
    main()
