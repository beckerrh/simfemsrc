#!/usr/bin/env python3

import sys, os
from argparse import ArgumentParser
from distutils.spawn import find_executable

if any(x not in os.listdir('.') for x in ['doc', 'install.py']):
	print('has to be called from source directory !')
	sys.exit(1)

packages="""
armadillo
boost-py
cmake
gmsh
open-mpi --with-cxx-bindings
doxygen
suite-sparse
vtk
"""
pippackages="""
numpy
scipy
matplotlib
pygmsh
"""

neededs = ['cmake','gmsh']
for needed in neededs:
    if find_executable(needed) is None:
        print("executable '{}' not found".format(needed))
        sys.exit(1)

sys.path.append(os.path.join(os.getcwd(),"python"))
from tools.compile import SimFemCompile

parser = ArgumentParser(description='install.')

startdir = os.getcwd()
startdirup = os.path.dirname(startdir)
installdir = os.path.join(startdirup, 'installdir')
builddir = os.path.join(startdirup, 'simfemsrc.compile')

compiler = SimFemCompile(installdir=installdir, builddir=builddir)
compiler.addArgumentsParser(parser, sys.argv[1:])

# compiler.compile(sourcedir=os.path.join(startdir, 'python'), clean=True)
# compiler.compile(sourcedir=os.path.join(startdir, 'lib'), clean=True)

compiler.compile(sourcedir=os.path.join(startdir, 'python'))
compiler.compile(sourcedir=os.path.join(startdir, 'lib'))
compiler.compile(sourcedir=os.path.join(startdir, 'bin'))

# import tools.simfempath as simfempath
# simfempath.storeSourcePath(installdir=installdir, simfemsourcedir=startdir)
