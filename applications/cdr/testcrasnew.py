from __future__ import print_function
import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("/Users/becker/Programs/simfem/installdir/lib")

import mesh.geometry
import tools.plot
import tools.testerrors
import numpy as np
import simfempy
import simfemcdr
import time

#---------------------------------------------------------#
def main():
    meshtypes = ["LineMesh"]
    testname = "poisson"
    testname = "cdrexp"
    # testname = "rdcosh"
    testtype = 'hmean'
    testtype = 'diff'
    methods=["newnitscheintegration", "nitscheintegration", "traditionalintegration"]
    paramdict = {}
    paramdict["fem"] = "P1"
    paramdict["deltasupg"] = 0.5
    paramdict["gamma"] = 2.0
    for meshtype in meshtypes:
        test(meshtype, testname, methods, paramdict, testtype)


#---------------------------------------------------------#
def test(meshtype, testname, methods, paramdict, testtype):
    if testname == "cdrexp":
        paramdict["diff"] = 0.001
        paramdict["beta"] = "east"
        paramdict["alpha"] = 0.0
        paramdict["application"]="cdexplayer"
    elif testname == "rdcosh":
        paramdict["symmetric"] = False
        paramdict["diff"] = 0.0001
        paramdict["beta"] = "zero"
        paramdict["alpha"] = 1.0
        paramdict["gamma"] = 2.0
        paramdict["application"]="rdcosh"
    elif testname == "poisson":
        paramdict["diff"] = 1.0
        paramdict["beta"] = "zero"
        paramdict["alpha"] = 0.0
        paramdict["application"]="quadratic"
        # application="cosinus"
        # application="linear"
    solver = simfemcdr.Solver()
    if testtype == 'hmean':
        hmean = None
        params=[0.5*0.5**i for i in range(10)]
        # params=[0.5*0.5**i for i in range(4)]
    elif testtype == 'diff':
        hmean=0.02
        params=[0.1*0.1**i for i in range(5)]
    simfemplot = tools.plot.SimFemPlot(methods, params=params, param=testtype)
    simfemplot.methodsnames["traditionalintegration"] = "trad"
    simfemplot.methodsnames["nitscheintegration"] = "nit"
    simfemplot.methodsnames["newnitscheintegration"] = "new"
    simfemplot.paramname = "h"
    simfemplot.order = False
    simfemplot.filename = testname+'.png'
    errs = ["L1"]
    simfemtesterrors = tools.testerrors.SimFemTestErrors(solver=solver, errs=errs, plotter=simfemplot, meshtype=meshtype, paramdict=paramdict, methods=methods, param=testtype, params=params, hmean=hmean, verbose=0)
    simfemtesterrors.run()


#---------------------------------------------------------#
if __name__ == '__main__':
    main()
