from __future__ import print_function
import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("@libpythonpath@")

import mesh.geometry
import tools.plot
import tools.testerrors
import numpy as np
import simfempy
import simfeminterface
import time


#---------------------------------------------------------#
def main():
    testtype = 'hmean'
    # testtype = 'gamma'
    # testtype = 'kex'
    # testtype = 'kin'
    # testtype = 'xgamma'
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
    paramdict = {}
    paramdict["fem"] = "P1"
    paramdict["application"] = "quadratic_circle"
    # paramdict["application"] = "linear_straight"
    # paramdict["application"] = "quadratic_straight"
    paramdict["kin"] = 1.0
    paramdict["kex"] = 100.0
    paramdict["gamma"] = 5.0
    methods=["nitsche","strong", "newnitsche","P1IFEM","CR1IFEM"]
    # methods=["nitsche","P1IFEM","CR1IFEM"]
    hmean=0.02
    newmeshperiteration=False
    if testtype == 'hmean':
        hmean = None
        params=[0.5*0.5**i for i in range(5)]
    elif testtype == 'gamma':
        params=[0.1* 4.0**i for i in range(15)]
    elif testtype == 'kex':
        params=[0.1* 5.0**i for i in range(12)]
    elif testtype == 'kin':
        paramdict["kex"] = 1.0
        params=[0.1* 5.0**i for i in range(12)]
    elif testtype == 'xgamma':
        newmeshperiteration=True
        paramdict["application"] = "quadratic_straight"
        params=[-0.1 + 0.02*i for i in range(11)]
    errs = ["H1", "E", "L2", "Linf"]
    simfemplot = tools.plot.SimFemPlot(methods, params=params, param=testtype)
    simfemplot.methodsnames["strong"] = "str"
    simfemplot.methodsnames["nitsche"] = "nit"
    simfemplot.methodsnames["P1IFEM"] = "p1iim"
    simfemplot.methodsnames["CR1IFEM"] = "cr1iim"
    simfemplot.methodsnames["nitschestabc"] = "sc"
    simfemplot.methodsnames["nitschestabnc"] = "snc"
    simfemplot.methodsnames["newnitsche"] = "new"
    simfemplot.methodsnames["newnitsche2"] = "new2"
    solver = simfeminterface.Solver()
    for meshtype in meshtypes:
        simfemtesterrors = tools.testerrors.SimFemTestErrors(solver=solver, errs=errs, plotter=simfemplot, meshtype=meshtype, paramdict=paramdict, methods=methods, param=testtype, params=params, hmean=hmean, newmeshperiteration=newmeshperiteration)
        simfemtesterrors.run()


#---------------------------------------------------------#
if __name__ == '__main__':
    main()


#
# #---------------------------------------------------------#
# def main():
#     testtype = 'hmean'
#     testtype = 'gamma'
#     testtype = 'kex'
#     testtype = 'xgamma'
#
#     paramdict = {}
#     paramdict["application"] = "quadratic_circle"
#     # paramdict["application"] = "linear_straight"
#     # paramdict["application"] = "quadratic_straight"
#     paramdict["kin"] = 1.0
#     paramdict["kex"] = 100.0
#     paramdict["gamma"] = 5.0
#
#     if testtype == 'hmean':
#         hmeans=[0.25*0.5**i for i in range(6)]
#         methods=["nitsche","nitschestabc","nitschestabnc","newnitsche"]
#         methods=["nitsche","CR1IFEM","P1IFEM"]
#         methods=["nitsche","newnitsche","newnitsche2"]
#         methods=["nitsche","strong","newnitsche2"]
#         testerrors(paramdict, methods, paramname='hmean', params=hmeans)
#     if testtype == 'gamma':
#         hmean=0.02
#         gammas=[0.1* 4.0**i for i in range(15)]
#         methods=["P1IFEM","CR1IFEM"]
#         methods=["nitsche","newnitsche","newnitsche2"]
#         methods=["nitsche","strong","newnitsche2"]
#         testerrors(paramdict, methods, paramname='gamma', params=gammas, hmean=hmean)
#     if testtype == 'kex':
#         hmean=0.02
#         kexs=[0.1* 5.0**i for i in range(12)]
#         methods=["P1IFEM"]
#         methods=["nitsche","CR1IFEM","P1IFEM"]
#         methods=["nitsche","newnitsche","newnitsche2"]
#         methods=["nitsche","strong","newnitsche2"]
#         testerrors(paramdict, methods, paramname='kin', params=kexs, hmean=hmean)
#     if testtype == 'xgamma':
#         hmean=0.02
#         paramdict["kin"] = 1.0
#         paramdict["kex"] = 10000.0
#         paramdict["application"] = "quadratic_straight"
#         params=[-0.1 + 0.2*i for i in range(11)]
#         methods=["nitsche","newnitsche","newnitsche2"]
#         methods=["nitsche","strong","CR1IFEM","P1IFEM"]
#         methods=["nitsche","strong","newnitsche2"]
#         # methods=["nitsche","strong"]
#         testerrors(paramdict, methods, paramname='xgamma', params=params, hmean=hmean, scale='toto')
#
# #---------------------------------------------------------#
# def testerrors(paramdict, methods, paramname, params, hmean=None, scale='loglog'):
#     meshtype = 'TriangleMesh'
#     from mesh.geometries.unitsquare import GeometryClass
#     errLinf = {}
#     errL2 = {}
#     errH1 = {}
#     errE = {}
#     times = {}
#     simfemplot = tools.plot.SimFemPlot(methods, params=params, param=paramname)
#     if paramname!='hmean':
#         geom=GeometryClass(hmean=hmean)
#         geom.runGmsh(outdir="Data")
#     for method in methods:
#         times[method] = 0.0
#         errL2[method] =  np.zeros(len(params))
#         errLinf[method] =  np.zeros(len(params))
#         errH1[method] =  np.zeros(len(params))
#         errE[method] =  np.zeros(len(params))
#     for ip,param in enumerate(params):
#         print("---- param=", param)
#         if paramname=='hmean':
#             geom=GeometryClass(hmean=param)
#             geom.runGmsh(outdir="Data")
#         for im,method in enumerate(methods):
#             print("####### %s #####" %(method))
#             start = time.time()
#             partion_id = 1
#             construct_bdrymeshes = False
#             mesh = simfempy.create(meshtype, partion_id, construct_bdrymeshes)
#             mesh.readGmsh("Data/"+geom.name+'.msh')
#             solver = simfeminterface.Solver()
#             solver.setMesh(mesh)
#             for key, value in paramdict.items():
#                 solver.setParameter(key,value)
#             solver.setParameter("method",method)
#             if paramname!='hmean':
#                 solver.setParameter(paramname,param)
#             solver.init()
#             info = solver.getInfo()
#             print(info)
#             d = solver.run()
#             errL2[method][ip] = d["L2"]
#             errLinf[method][ip] = d["Linf"]
#             errH1[method][ip] = d["H1"]
#             errE[method][ip] = d["E"]
#             # print("d=",d)
#             solver.writeXdmf()
#         times[method] = time.time()-start
#     for method in methods:
#         print("errL2 %-20s" %(method), errL2[method])
#     print()
#     for method in methods:
#         print("errLinf %-20s" %(method), errLinf[method])
#     print()
#     for method in methods:
#         print("errH1 %-20s" %(method), errH1[method])
#     for method in methods:
#         print("errE %-20s" %(method), errE[method])
#     # for method in methods:
#     #     errH1[method] = np.sqrt( eps*errH1[method]**2 + errL2[method]**2)
#     for method in methods: print("%-20s %10g" %(method,times[method]))
#
#     datatoplot = {}
#     # datatoplot['L1'] = errL1
#     datatoplot['L2'] = errL2
#     datatoplot['Linf'] = errLinf
#     datatoplot['H1'] = errH1
#     datatoplot['E'] = errE
#     methodsnames={}
#     for method in methods: methodsnames[method]=method.lower()
#     methodsnames["P1IFEM"] = "p1iim"
#     methodsnames["CR1IFEM"] = "cr1iim"
#     methodsnames["nitsche"] = "nit"
#     methodsnames["nitschestabc"] = "sc"
#     methodsnames["nitschestabnc"] = "snc"
#     methodsnames["newnitsche"] = "new"
#     methodsnames["newnitsche2"] = "new2"
#
#     simfemplot.ploterrors(datatoplot, scale=scale)
#
#
# #---------------------------------------------------------#
# if __name__ == '__main__':
#     main()
