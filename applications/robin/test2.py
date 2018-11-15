from __future__ import print_function
import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("@libpythonpath@")

import mesh.geometry
import tools.plot
import numpy as np
import simfem
import simfemrobin
import time

#---------------------------------------------------------#
def main():
    testtype = 'hmean'
    # testtype = 'gamma'
    # if len(sys.argv)>1:
    #     if sys.argv[1] == "1D":
    #         meshtypes = ["LineMesh"]
    #     elif sys.argv[1] == "2D":
    #         meshtypes = ["TriangleMesh"]
    #     elif  sys.argv[1] == "3D":
    #         meshtypes = ["TetrahedralMesh"]
    #     else:
    #         raise KeyError("unknwon argument", sys.argv[1])
    # meshtypes = ["LineMesh", "TriangleMesh", "TetrahedralMesh"]
    meshtype = "TriangleMesh"
    testname = "poisson"
    # testname = "cdrexp"
    # testname = "rdcosh"
    fem = "P1"
    # fem = "CR1"
    application="quadratic"
    paramdict = {}
    paramdict["diff"] = 1.0
    paramdict["gamma"] = 5.0
    paramdict["alpha"] = 1.0
    paramdict["fem"] = fem
    paramdict["application"] = application

    if testtype == 'hmean':
        hmeans=[0.5*0.5**i for i in range(4)]
        methods=["traditional", "nitsche"]
        testerrors(meshtype, paramdict, methods, paramname='hmean', params=hmeans)
    if testtype == 'gamma':
        hmean=0.02
        gammas=[0.1* 4.0**i for i in range(15)]
        methods=["traditional", "nitsche"]
        testerrors(meshtype, paramdict, methods, paramname='gamma', params=gammas, hmean=hmean)


#---------------------------------------------------------#
#---------------------------------------------------------#
def testerrors(meshtype, paramdict, methods, paramname, params, hmean=None):
    from mesh.geometries.unitsquare import GeometryClass
    errLinf = {}
    errL1 = {}
    errL2 = {}
    errH1 = {}
    times = {}
    simfemplot = tools.plot.SimFemPlot()
    if paramname!='hmean':
        geom=GeometryClass(hmean=hmean)
        geom.runGmsh(outdir="Data")
    for method in methods:
        times[method] = 0.0
        errL2[method] =  np.zeros(len(params))
        errL1[method] =  np.zeros(len(params))
        errLinf[method] =  np.zeros(len(params))
        errH1[method] =  np.zeros(len(params))
    for ip,param in enumerate(params):
        print("---- param=", param)
        if paramname=='hmean':
            geom=GeometryClass(hmean=param)
            geom.runGmsh(outdir="Data")
        for im,method in enumerate(methods):
            print("####### %s #####" %(method))
            start = time.time()
            partion_id = 1
            construct_bdrymeshes = False
            mesh = simfem.create(meshtype, partion_id, construct_bdrymeshes)
            mesh.readGmsh("Data/"+geom.name+'.msh')
            mesh.addGeometryObject('MeasureOfCell')
            mesh.addGeometryObject('Normals')
            solver = simfemrobin.Solver()
            solver.setMesh(mesh)
            for key, value in paramdict.items():
                solver.setParameter(key,value)
            solver.setParameter("method",method)
            if paramname!='hmean':
                solver.setParameter(paramname,param)
            solver.init()
            info = solver.getInfo()
            print(info)
            d = solver.run()
            errL1[method][ip] = d["L1"]
            errL2[method][ip] = d["L2"]
            errLinf[method][ip] = d["Linf"]
            errH1[method][ip] = d["H1"]
            # print("d=",d)
            solver.writeXdmf()
        times[method] = time.time()-start
    # for method in methods:
    #     print("errL2 %-20s" %(method), errL2[method])
    # for method in methods:
    #     print("errLinf %-20s" %(method), errLinf[method])
    # for method in methods:
    #     print("errH1 %-20s" %(method), errH1[method])
    # for method in methods:
    #     errH1[method] = np.sqrt( eps*errH1[method]**2 + errL2[method]**2)
    for method in methods: print("%-20s %10g" %(method,times[method]))

    datatoplot = {}
    # datatoplot['L1'] = errL1
    datatoplot['L2'] = errL2
    datatoplot['Linf'] = errLinf
    datatoplot['H1'] = errH1
    datatoplot['L1'] = errL1
    methodsnames={}
    for method in methods: methodsnames[method]=method.lower()
    methodsnames["nitsche"] = "nit"
    methodsnames["traditional"] = "tra"

    simfemplot.ploterrors(datatoplot, methods, params, methodsnames, paramname)


#---------------------------------------------------------#
if __name__ == '__main__':
    main()
