from __future__ import print_function
import pygmsh
import os, subprocess

#------------------------------------------------------------------------
class Geometry(pygmsh.built_in.Geometry):
  def __init__(self, **kwargs):
    self.datatype = 'binary'
    self.npartitions = 0
    self.hmean = 1.0
    self.dim = 3
    if 'datatype' in kwargs: self.datatype = kwargs.pop('datatype')
    if 'name' in kwargs: self.name = kwargs.pop('name')
    if 'npartitions' in kwargs: self.npartitions = kwargs.pop('npartitions')
    if 'hmean' in kwargs: self.hmean = kwargs.pop('hmean')
    if 'dim' in kwargs: self.dim = kwargs.pop('dim')
    print("@@@@@@@@  hmean",self.hmean)
    pygmsh.built_in.Geometry.__init__(self, kwargs)
    self.gmsh_executable = pygmsh.helpers._get_gmsh_exe()
    self.defineGeometry()
#------------------------------------------------------------------------
  def runGeo(self, outdir=None, name=None):
    if outdir is None : outdir = os.getcwd()
    if name is None : name = self.name
    # print("@@@@@@@@ name %s" %name)
    if not os.path.isdir(outdir) : os.mkdir(outdir)
    filenamegeo = os.path.join(outdir,name + '.geo')
    # print("@@@@@@@@ filenamegeo %s" %filenamegeo)
    file = open(filenamegeo, "w")
    # file.write(self.get_code().encode())
    file.write(self.get_code())
    file.close()

  def runGmshRefine(self, number, outdir=None, name=None):
    if outdir is None : outdir = os.getcwd()
    if name is None : name = self.name
    namecoarse = name + "_{:03d}".format(number)
    nameref = name + "_{:03d}".format(number+1)
    filenamemsh = os.path.join(outdir,namecoarse+'.msh')
    filenamemshref = os.path.join(outdir,nameref+'.msh')
    print("filenamemsh", filenamemsh)
    print("filenamemshref", filenamemshref)
    cmd = [self.gmsh_executable]
    cmd += ["-refine"]
    cmd += ["{}".format(filenamemsh)]
    cmd += ["-o"]
    cmd += ["{}".format(filenamemshref)]
    cmd += ["-format"]
    cmd += ["{}".format("msh2")]
    p = subprocess.Popen(cmd, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr != "":
      print("cmd=",cmd)
      print("cmd=",' '.join(cmd))
      print("stderr=",stderr)
      raise RuntimeError('Gmsh exited with error (return code %d).' %p.returncode)
    return nameref

  def runGmsh(self, outdir=None, name=None, verbose = True, number=None):
    if outdir is None : outdir = os.getcwd()
    if name is None : name = self.name
    if number is not None: name += "_{:03d}".format(number)
    filenamegeo = os.path.join(outdir,name + '.geo')
    self.runGeo(outdir, name)
    filenamemsh = os.path.join(outdir,name+'.msh')
    # cmd = [self.gmsh_executable, '-'+str(self.dim), filenamegeo, '-o', filenamemsh]
    cmd = [self.gmsh_executable]
    cmd.append("-{:1d}".format(self.dim))
    cmd.append("{}".format(filenamegeo))
    cmd.append("-o")
    cmd.append("{}".format(filenamemsh))
    # cmd.append("-format")
    # cmd.append("{}".format("msh2"))
    if self.npartitions>0:
      cmd.append('-part', str(self.npartitions), '-oneFilePerPart')
    # if self.datatype != 'ascii':
    #   cmd += ['-bin']
    # print("### Format is ASCII for the moment, need to work on format 4")
    p = subprocess.Popen(cmd, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    # if verbose:
    #   while True:
    #     line = p.stdout.readline()
    #     if not line: break
    #     print(line, end='')
    if stderr:
      print("cmd=",cmd)
      print("cmd=",' '.join(cmd))
      # print("stdout=",stdout)
      print("stderr=",stderr)
      raise RuntimeError('Gmsh exited with error (return code %d).' %p.returncode)
    return name
