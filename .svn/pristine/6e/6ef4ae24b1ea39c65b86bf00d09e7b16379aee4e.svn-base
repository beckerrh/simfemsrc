from __future__ import print_function
import sys, os, shutil
import tools.plot
import mesh.geometry
import numpy as np
import time
import simfempy
import copy

#------------------------------------------------------------------------
class SimFemTestErrors(object):
    """TestErrors.
    """
#------------------------------------------------------------------------
    def makeMesh(self, hmean):
        if  self.meshtype == "LineMesh":
            from mesh.geometries.unitline import GeometryClass
        elif  self.meshtype == "TriangleMesh":
            from mesh.geometries.unitsquare import GeometryClass
        elif  self.meshtype == "TetrahedralMesh":
            from mesh.geometries.unitcube import GeometryClass
        self.geom=GeometryClass(hmean=hmean)
        self.geom.runGmsh(outdir="Data")
        partion_id = 1
        construct_bdrymeshes = False
        self.mesh = simfempy.create(self.meshtype, partion_id, construct_bdrymeshes)
        self.mesh.readGmsh("Data/"+self.geom.name+'.msh')
        self.mesh.addGeometryObject('MeasureOfCell')
        self.mesh.addGeometryObject('Normals')
        self.mesh.save("Data/"+self.geom.name+'.h5')
#------------------------------------------------------------------------
    def __init__(self, solver, errs, plotter, meshtype, paramdict, methods, param, params, hmean=None, newmeshperiteration=False, verbose=0):
        self.solver = solver
        self.errs = errs
        self.plotter = plotter
        self.meshtype = meshtype
        self.paramdict = paramdict
        self.methods = methods
        self.param = param
        self.params = params
        self.hmean = hmean
        self.newmeshperiteration = newmeshperiteration
        self.verbose = verbose
        if param=='hmean':
            self.newmeshperiteration=True
        self.mesher = self.makeMesh
        self.errors = {}
        for err in self.errs:
            self.errors[err] = {}
        self.times = {}
        for method in self.methods: self.times[method] = 0.0
        for method in self.methods:
            self.times[method] = 0.0
            for err in self.errs:
                self.errors[err][method] =  np.zeros(len(params))
#---------------------------------------------------------#
    def run(self):
        if not self.newmeshperiteration:
            self.mesher(hmean=self.hmean)
            self.solver.loadMesh(self.meshtype, "Data/"+self.geom.name+'.h5')
        for ip,param in enumerate(self.params):
            print("@@@@@@ ", self.param, "=", param, "@@@@@@")
            if self.newmeshperiteration:
                if self.param=='hmean':
                    self.mesher(hmean=param)
                else:
                    self.mesher(hmean=self.hmean)
                self.solver.loadMesh(self.meshtype, "Data/"+self.geom.name+'.h5')
            for im,method in enumerate(self.methods):
                print("\t###### %s ######" %(method))
                start = time.time()
                # solver = self.solver
                for key, value in self.paramdict.items():
                    self.solver.setParameter(key,value)
                self.solver.setParameter("method",method)
                if self.param!='hmean':
                    self.solver.setParameter(self.param,param)
                # solver.setMesh(mesh)
                # self.solver.loadMesh(self.meshtype, "Data/"+self.geom.name+'.h5')
                self.solver.init()
                if self.verbose:
                    meshinfo = self.solver.getMeshInfo()
                    print("meshinfo",meshinfo)
                    if self.verbose>1:
                        info = self.solver.getInfo()
                        print("solverinfo",info)
                d = self.solver.run()
                for err in self.errs:
                    self.errors[err][method][ip] =  d[err]
                # print("d=",d)
                self.solver.writeXdmf()
                self.times[method] += time.time()-start
        for err in self.errs:
            for method in self.methods:
                print("err%s %s %-20s" %( err, method, self.errors[err][method]))
        for method in self.methods: print("%-20s %10g" %(method,self.times[method]))
        datatoplot = {}
        for err in self.errs: datatoplot[err] = self.errors[err]
        self.plotter.ploterrors(datatoplot, keys=self.errs)
