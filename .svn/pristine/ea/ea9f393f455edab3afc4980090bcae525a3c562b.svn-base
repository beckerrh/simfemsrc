from __future__ import print_function
import sys, os, glob, shutil, subprocess
import argparse
from tools.compile import SimFemCompile

#------------------------------------------------------------------------
class SimFemRun(object):
    def __init__(self, simfemsrcdir):
        # print ('simfemsrcdir', simfemsrcdir)
        self.simfemsrcdir = simfemsrcdir
        parser = argparse.ArgumentParser(
            # description='utility',
            usage='''simfem <command> [<args>]
                    commands are:
                    compile     compile library and/or project
                    fromlib     copy py and cpp files from library/application to project
                    tolib       copy py and cpp files from project to library/application
                    gmsh        run gmsh
                    test        run test
                    progtest    run progtest
                    application run application
                    ''')
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print ('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    def fromlib(self):
        parser = argparse.ArgumentParser(description='Run gmsh')
        parser.add_argument('--nopy', default = False, action="store_true", help='only py')
        parser.add_argument('--nocpp', default = False, action="store_true", help='only cpp')
        name = os.path.basename(os.getcwd())
        args = parser.parse_args(sys.argv[2:])
        self._fromlib(name, vars(args))
    def tolib(self):
        parser = argparse.ArgumentParser(description='Run gmsh')
        parser.add_argument('--nopy', default = False, action="store_true", help='only py')
        parser.add_argument('--nocpp', default = False, action="store_true", help='only cpp')
        parser.add_argument('--dry', default = False, action="store_true", help='only compare')
        name = os.path.basename(os.getcwd())
        args = parser.parse_args(sys.argv[2:])
        self._tolib(name, vars(args))
    def gmsh(self):
        # print 'running gmsh'
        parser = argparse.ArgumentParser(description='Run gmsh')
        parser.add_argument('scriptandargs', help='script and args to launch', nargs='*')
        parser.add_argument('-hmean', default = 0.8, type=float, help='mesh width')
        parser.add_argument('-nparts', default = 0, type=int, help='number of partitions')
        parser.add_argument('-datatype', default = "binary", type=str, help='dataype')
        parser.add_argument('-localdir', default = "", type=str, help="local directory (if '' we use installdir)")
        parser.add_argument('--onlygeo', default = False, action="store_true", help='create only geometry')
        args = parser.parse_args(sys.argv[2:])
        self._gmsh(vars(args))
#------------------------------------------------------------------------
    def compile(self):
        # print 'running compile'
        self.compiler = SimFemCompile(simfemsourcedir=self.simfemsrcdir)
        parser = argparse.ArgumentParser(description='Run compile')
        parser.add_argument('scriptandargs', help='script and args to launch', nargs='*')
        parser.add_argument('--all', default = True, action="store_true", help='compile library and project')
        parser.add_argument('--project', default = False, action="store_true", help='compile project')
        parser.add_argument('--doc', default = False, action="store_true", help='generate doc')
        parser.add_argument('--external', default = False, action="store_true", help='external install')
        parser.add_argument('-installdir', default = "", type=str, help="install directory (if '' we use installdir)")
        self.compiler.addArgumentsParser(parser,sysargs=sys.argv[2:])
        args = parser.parse_args(sys.argv[2:])
        self._compile(vars(args))
#------------------------------------------------------------------------
    def test(self):
        # print 'running test'
        self.compiler = SimFemCompile(simfemsourcedir=self.simfemsrcdir)
        parser = argparse.ArgumentParser(description='Run test')
        parser.add_argument('scriptandargs', help='script and args to launch', nargs='*')
        parser.add_argument('-nmpi', default = 0, type=int, help='number of mpi processes')
        parser.add_argument('--py', default = False, action="store_true", help='test python')
        parser.add_argument('--cpp', default = False, action="store_true", help='test cpp')
        self.compiler.addArgumentsParser(parser,sysargs=sys.argv[2:])
        args = parser.parse_args(sys.argv[2:])
        self._test(vars(args))
#------------------------------------------------------------------------
    def application(self):
        # print 'running test'
        self.compiler = SimFemCompile(simfemsourcedir=self.simfemsrcdir)
        parser = argparse.ArgumentParser(description='Run application')
        parser.add_argument('scriptandargs', help='script and args to launch', nargs='*')
        parser.add_argument('-nmpi', default = 0, type=int, help='number of mpi processes')
        parser.add_argument('--cmake', default = False, action="store_true", help='generate cmake')
        parser.add_argument('--clean', default = False, action="store_true", help='clean local dir')
        parser.add_argument('--py', default = False, action="store_true", help='python')
        self.compiler.addArgumentsParser(parser,sysargs=sys.argv[2:])
        args = parser.parse_args(sys.argv[2:])
        self._application(vars(args))
#------------------------------------------------------------------------
    def progtest(self):
        self.compiler = SimFemCompile(simfemsourcedir=self.simfemsrcdir)
        parser = argparse.ArgumentParser(description='Run progtest')
        parser.add_argument('scriptandargs', help='script and args to launch', nargs='*')
        parser.add_argument('--clean', default = False, action="store_true", help='clean local dir')
        self.compiler.addArgumentsParser(parser,sysargs=sys.argv[2:])
        args = parser.parse_args(sys.argv[2:])
        self._progtest(vars(args))
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    def _gmsh(self, args):
        # print 'Running simfem _gmsh, args=%s' % args
        hmean = args['hmean']
        nparts = args['nparts']
        datatype = args['datatype']
        localdir = args['localdir']
        print ('hmean=', hmean, 'nparts', nparts, 'datatype', datatype)
        import glob
        geomdir = os.path.join(os.path.join(os.path.join(self.simfemsrcdir, 'python'), 'mesh'), 'geometries')
        if len(args['scriptandargs'])==0:
            items = [os.path.basename(x).split('.')[0] for x in glob.glob(geomdir+'/*.py')]
            itemsfiltered = [x for x in items if x[0]!="_"]
            print ('available packages in %s:\n\t%s' %(geomdir, '\n\t'.join(itemsfiltered)))
            sys.exit(1)
        startdirup = os.path.dirname(self.simfemsrcdir)
        if localdir:
            outdir = localdir
        else:
            outdir = os.path.join(os.path.join(startdirup, 'installdir'),'meshes')
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
        geomname = args['scriptandargs'][0]
        try:
            sys.path.insert(0,geomdir)
            import importlib
            module = importlib.import_module(geomname, package=None)
            print ('module', module)
        except:
          print ('### import failed\ngeomname', geomname)
          print ('geomdir', geomdir)
          print ('sys.path = ', sys.path)
          message = "could not find python script for '"+geomname+"'\nfound:\n" + '\n'.join(glob.glob(geomdir+os.path.sep+'*.py'))
          raise ValueError(message)
        # print ('module', module)
        # print ('dir(module)', dir(module))
        geom = module.GeometryClass(hmean=hmean, npartitions=nparts, datatype=datatype)
        if args['onlygeo']:
            geom.runGeo(outdir=outdir)
        else:
            geom.runGmsh(outdir=outdir)
#------------------------------------------------------------------------
    def _compilelibs(self):
        self.compiler.compilelib(sourcedir=os.path.join(self.simfemsrcdir, 'python'))
        self.compiler.compilelib(sourcedir=os.path.join(self.simfemsrcdir, 'lib'))
        # self.compiler.compilelib(sourcedir=os.path.join(self.simfemsrcdir, 'libcpp'))
        # self.compiler.compilelib(sourcedir=os.path.join(self.simfemsrcdir, 'libpy'))
    def _compile(self, args):
        # print 'Running simfem _compile, args=%s' % args
        if args['external']:
            if args['installdir']: self.compiler.installdir = args['installdir']
            self.compiler.externalInstall(args['scriptandargs'])
            return
        if args['doc']:
            self.compiler.generateDoc()
            return
        self._compilelibs()
        if args['project'] or args['all']:
            if 'CMakeLists.txt' not in os.listdir('.'):
                raise KeyError("no local project to compile")
            self.compiler.compiletest(sourcedir=os.getcwd())
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    def _convert(self, filename, replaceDict):
        inputfile = open(filename,'r')
        input = inputfile.read()
        inputfile.close()
        for k, v in replaceDict.items():
            input = input.replace( k, v )
        # print("input replaced", input)
        outputfile = open(filename,'w')
        outputfile.write(input)
        outputfile.close()
    def _mkdirs(self, name, args, cleanlocaldir=True):
        dirsrc = os.path.join(self.simfemsrcdir, name)
        if len(args['scriptandargs'])==0:
          items = next(os.walk(dirsrc))[1]
          print ('available tests in %s:\n%s' %(dirsrc, '\n'.join(items)))
          sys.exit(1)
        dirname = args['scriptandargs'][0]
        startdirup = os.path.dirname(self.simfemsrcdir)
        localdir = os.path.join(startdirup, dirname)
        if cleanlocaldir: shutil.rmtree(localdir, ignore_errors=True)
        from distutils import dir_util
        dir_util.copy_tree(os.path.join(dirsrc, dirname), localdir)
        # shutil.copytree(os.path.join(dirsrc, dirname), localdir)
        os.chdir(localdir)
        libpath = os.path.join(self.compiler.installdir, 'lib')
        simfempythonpath = os.path.join(self.simfemsrcdir, 'python')
        self._compilelibs()
        replaceDict={}
        replaceDict['@simfempythonpath@'] = simfempythonpath
        replaceDict['@libpythonpath@'] = libpath
        replaceDict['@simfeminstallpath@'] = self.compiler.installdir
        return (dirname, localdir, replaceDict)
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    def _fromlib(self, name, args):
        libdir = os.path.join(self.simfemsrcdir, "applications", name)
        projdir = os.getcwd()
        pyfiles = glob.glob(os.path.join(libdir,"*.py"))
        replaceDict={}
        replaceDict['@simfempythonpath@'] = os.path.join(self.simfemsrcdir, 'python')
        for libfile in pyfiles:
            basename = os.path.basename(libfile)
            shutil.copy(libfile, projdir)
            projfile = os.path.join(projdir, basename)
            self._convert(projfile, replaceDict)
        from distutils import dir_util
        dir_util.copy_tree(os.path.join(libdir, "cpp"), os.path.join(projdir, "cpp"))
        dir_util.copy_tree(os.path.join(libdir, "cppy"), os.path.join(projdir, "cppy"))
#------------------------------------------------------------------------
    def _tolib(self, name, args):
        # print("args", args)
        import filecmp, difflib
        if not args['nopy']:
            libdir = os.path.join(self.simfemsrcdir, "applications", name)
            projdir = os.getcwd()
            filestocopy = self._comparefiles(projdir, libdir, suffix='*.py', ignores='sys.path.append')
            # print("filestocopy", filestocopy)
            if not args['dry']:
                replaceDict={}
                replaceDict[os.path.join(self.simfemsrcdir, 'python')] = '@simfempythonpath@'
                for file in filestocopy:
                    projfile = os.path.join(projdir, file)
                    libfile = os.path.join(libdir, file)
                    shutil.copyfile(projfile, libfile)
                    self._convert(libfile, replaceDict)
        if not args['nocpp']:
            appdir = os.path.join(self.simfemsrcdir, "applications", name)
            for dir in ['cpp', 'cppy']:
                libdir = os.path.join(appdir, dir)
                projdir = os.path.join(os.getcwd(), dir)
                filestocopy = self._comparefiles(projdir, libdir, suffix=['*.cpp', '*.hpp'])
                # print("filestocopy", filestocopy)
                if not args['dry']:
                    for file in filestocopy:
                        projfile = os.path.join(projdir, file)
                        libfile = os.path.join(libdir, file)
                        shutil.copyfile(projfile, libfile)
#------------------------------------------------------------------------
    def _comparefiles(self, projdir, libdir, suffix, ignores=[]):
        import filecmp, difflib
        if not isinstance(ignores, (list,)): ignores = list(ignores)
        if isinstance(suffix, (list,)):
            libfiles = [f for suff in suffix for f in glob.glob(os.path.join(libdir,suff))]
        else:
            libfiles = glob.glob(os.path.join(libdir,suffix))
        for libfile in libfiles:
            basename = os.path.basename(libfile)
            projfile = os.path.join(projdir, basename)
            if not os.path.isfile(projfile): print("file '{}'' does not exist in project".format(projfile))
        if isinstance(suffix, (list,)):
            projfiles = [f for suff in suffix for f in glob.glob(os.path.join(projdir,suff))]
        else:
            projfiles = glob.glob(os.path.join(projdir,suffix))
        filestocopy = []
        for projfile in projfiles:
            basename = os.path.basename(projfile)
            libfile = os.path.join(libdir, basename)
            if not os.path.isfile(libfile):
                print("file '{}'' does not exist in lib".format(libfile))
                filestocopy.append(basename)
        filescmp = [os.path.basename(x) for x in projfiles]
        # print("filescmp={}".format(filescmp))
        match, mismatch, errors = filecmp.cmpfiles(projdir, libdir, filescmp)
        # print("match={} mismatch={} errors={}".format(match, mismatch, errors))
        for m in mismatch:
            # print("m={}:".format(m))
            projfile = open(os.path.join(projdir, m))
            libfile = open(os.path.join(libdir, m))
            result = difflib.unified_diff(projfile.readlines(), libfile.readlines(),n=0)
            copifile=False
            for line in result:
                checkline=True
                for ignore in ignores:
                    if ignore in line: checkline=False
                for prefix in ('---', '+++', '@@'):
                    if line.startswith(prefix): checkline=False
                if checkline:
                    print("{}\t {}".format(m,line[1:]), end = '')
                    copifile=True
            if copifile: filestocopy.append(m)
        return filestocopy
#------------------------------------------------------------------------
    def _test_py(self, replaceDict):
        self._convert(filename="test.py", replaceDict=replaceDict)
        command = "python test.py"
        returncode = subprocess.check_call(command, shell=True)
    def _test_cpp(self, args, replaceDict, testname, localdir):
        self._convert(filename="test.cpp", replaceDict=replaceDict)
        self.compiler.writeCmakeTest(testname=testname)
        self.compiler.compiletest(sourcedir=localdir)
        nmpi = args['nmpi']
        if nmpi:
            cmd = ['mpirun', '--oversubscribe', '-n', str(nmpi), './'+testname]
        else:
            cmd = ['./'+testname]
        print ('command', ' '.join(cmd))
        returncode = subprocess.check_call(' '.join(cmd), shell=True)
    def _test(self, args):
        (testname, localdir, replaceDict) = self._mkdirs(name='tests', args=args)
        os.chdir(localdir)
        if args['py']:
            self._test_py(replaceDict=replaceDict)
        if args['cpp']:
            self._test_cpp(args=args, replaceDict=replaceDict, testname=testname, localdir=localdir)
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    def _application_buildcpp(self, args, applicationname, cmake, localdir):
        withcppy = os.path.isdir('cppy')
        if cmake:
            self.compiler.writeCmakeApplication(name=applicationname, withcppy=withcppy)
        self.compiler.compiletest(sourcedir=localdir)
    def _application_runcpp(self, args, applicationname):
        nmpi = args['nmpi']
        if nmpi:
            cmd = ['mpirun', '--oversubscribe', '-n', str(nmpi), './'+applicationname]
        else:
            cmd = ['./'+applicationname]
        print ('command', ' '.join(cmd))
        returncode = subprocess.check_call(' '.join(cmd), shell=True)
    def _application_py(self):
        command = "python test.py"
        returncode = subprocess.check_call(command, shell=True)
    def _application(self, args):
        (applicationname, localdir, replaceDict) = self._mkdirs(name='applications', args=args,cleanlocaldir=args['clean'])
        self._convert(filename="cpp/main.cpp", replaceDict=replaceDict)
        os.chdir(localdir)
        cmake = args['clean'] or args['cmake'] or not os.path.isfile('CMakeLists.txt')
        self._application_buildcpp(args, applicationname, cmake, localdir)
        for file in glob.glob("*.py"):
            self._convert(filename=file, replaceDict=replaceDict)
        if args['py']:
            self._application_py()
        else:
            self._application_runcpp(args, applicationname)
#------------------------------------------------------------------------
#------------------------------------------------------------------------
    def _progtest(self, args):
        (testname, localdir, replaceDict) = self._mkdirs(name='progtest', args=args, cleanlocaldir=args['clean'])
        dirsrc = os.path.join(self.simfemsrcdir, 'progtest')
        srcdir = os.path.join(dirsrc, testname)
        os.chdir(localdir)
        if os.path.isfile("nocmake"):
            pass
        else:
            self.compiler.writeCmakeProgtest(name=testname, srcdir=srcdir)
        self.compiler.compiletest(sourcedir=localdir)
