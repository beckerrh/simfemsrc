from __future__ import print_function
import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("@libpythonpath@")

import tools.plot
import tools.testerrors
import simfempy
import simfemrobin

#---------------------------------------------------------#
def main():
    testtype = 'hmean'
    # testtype = 'robin'
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
    paramdict["application"] = "quadratic"
    paramdict["application"] = "cosinus"
    paramdict["diff"] = 1.0
    paramdict["alpha"] = 1.0
    paramdict["gamma"] = 0.1
    paramdict["robin"] = 100.
    paramdict["beta"] = "zero"
    methods=["traditional", "nitsche"]
    methods=["nitsche", "traditional"]
    # methods=["traditional"]
    # methods=["nitsche"]
    if testtype == 'hmean':
        hmean = None
        params=[0.5*0.5**i for i in range(4)]
    elif testtype == 'robin':
        hmean=0.02
        params=[0.001* 5.0**i for i in range(15)]
    errs = ["L2", "Linf", "H1"]
    simfemplot = tools.plot.SimFemPlot(methods, params=params, param=testtype)
    simfemplot.methodsnames["traditional"] = "trad"
    simfemplot.methodsnames["nitsche"] = "nit"
    solver = simfemrobin.Solver()
    if testtype == 'robin':
        simfemplot.paramname = r"1/$\varepsilon$"
    for meshtype in meshtypes:
        simfemtesterrors = tools.testerrors.SimFemTestErrors(solver=solver, errs=errs, plotter=simfemplot, meshtype=meshtype, paramdict=paramdict, methods=methods, param=testtype, params=params, hmean=hmean)
        simfemtesterrors.run()


#---------------------------------------------------------#
if __name__ == '__main__':
    main()
