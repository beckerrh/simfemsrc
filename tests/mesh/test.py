from __future__ import print_function
import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("@libpythonpath@")

import mesh.geometry
import numpy as np
import simfem

#---------------------------------------------------------#
def main():
    if len(sys.argv)>1:
        if sys.argv[1] == "1D":
            meshnames = ["LineMesh"]
        elif sys.argv[1] == "2D":
            meshnames = ["TriangleMesh"]
        elif  sys.argv[1] == "3D":
            meshnames = ["TetrahedralMesh"]
        else:
            raise KeyError("unknwon argument", sys.argv[1])
    else:
        meshnames = ["LineMesh", "TriangleMesh", "TetrahedralMesh"]
    for meshname in meshnames:
        test(meshname)

#---------------------------------------------------------#
def test(meshname):
    if  meshname == "LineMesh":
        from mesh.geometries.unitline import GeometryClass
    elif  meshname == "TriangleMesh":
        from mesh.geometries.unitsquare import GeometryClass
    elif  meshname == "TetrahedralMesh":
        from mesh.geometries.unitcube import GeometryClass
    geom=GeometryClass(hmean=0.5)
    geom.runGmsh()
    partion_id = 1
    construct_bdrymeshes = True
    mesh = simfem.create(meshname, partion_id, construct_bdrymeshes)
    mesh.readGmsh(geom.name+'.msh')
    mesh.addGeometryObject('MeasureOfCell')
    mesh.addGeometryObject('Normals')
    mesh.save(geom.name+'.h5')
    import mesh.meshtest as meshtest
    print("%s IO OK ? %s" %(meshname, meshtest.testIo(mesh)))
    vtkname = geom.name+'.vtk'
    vtkbdryname = geom.name+'_boundary.vtk'
    mesh.writeVtk(vtkname)
    mesh.writeBoundaryVtk(vtkbdryname)

#---------------------------------------------------------#
if __name__ == '__main__':
    main()
