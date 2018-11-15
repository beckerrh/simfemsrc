#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser

sys.path.append(os.path.join(os.getcwd(),"python"))
from tools.compile import SimFemCompile

compiler = SimFemCompile(simfemsourcedir=os.getcwd())

compiler.externalInstall(sys.argv[1:])
