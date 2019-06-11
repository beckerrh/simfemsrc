import os, sys
simfempythondir = '../simfemsrc/python'
sys.path.append(simfempythondir)
from tools.runsubprocess import run

bindir = '../installdir/bin'

command = ["ReaderMesh ", "BackwardFacingStep", "meshquad", "quad", "ascii"] 
returncode = run(command)
command = ["ReaderMesh ", "meshquad", "meshtri", "quadtotri", "ascii"] 
returncode = run(command)
