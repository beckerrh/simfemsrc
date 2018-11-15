import sys, os, shutil
from argparse import ArgumentParser

sys.path.append("@simfempythonpath@")
sys.path.append("@libpythonpath@")

import mesh.geometry
import numpy as np
import simfem

if __name__ == '__main__':
    geomsdir = os.path.join(os.path.join("@simfempythonpath@", 'mesh'),'geometries')
    print 'geomsdir', geomsdir
    sys.path.insert(0,geomsdir)
    import importlib
    geommodule = importlib.import_module('unitsquare')
    hmean = 1.0
    nparts = 2
    geom = geommodule.GeometryClass(hmean=hmean, npartitions=nparts, datatype='ascii')
    geom.runGmsh()
