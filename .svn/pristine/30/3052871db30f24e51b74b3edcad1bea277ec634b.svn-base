from __future__ import print_function
import sys, os, shutil
import tools.runsubprocess as runsubprocess
from argparse import ArgumentParser

#------------------------------------------------------------------------
class SimFemCompile(object):
    """Perform different compiling tasks.
    Keyword arguments:
    installdir         -- directory where to install our results
    builddir           -- directory where to build
    simfemsourcedir    -- the directory of the simfem sources
    simfeminstalldir   -- the directory where simfem is installed
    debug              -- try to help development
    runoptions         -- control ouput
    """

#------------------------------------------------------------------------
    def __init__(self, installdir=None, builddir=None, simfemsourcedir=None, simfeminstalldir=None, debug=False):
        self.CMAKE_BUILD_TYPES = ['Release','RelWithDebInfo','Debug','DebugFull','Profile']
        self.debug=debug
        self.installdir = installdir
        self.builddir = builddir
        self.simfemsourcedir = simfemsourcedir
        self.simfeminstalldir = simfeminstalldir
        self.runoptions = None
        if simfemsourcedir:
            startdirup = os.path.dirname(simfemsourcedir)
            if not builddir: self.builddir = os.path.join(startdirup, 'simfemsrc.compile')
            if not installdir: self.installdir = os.path.join(startdirup, 'installdir')
            if not simfeminstalldir: self.simfeminstalldir = os.path.join(startdirup, 'installdir')
        self.LIBRARYCPPY = "simfempy"

#------------------------------------------------------------------------
    def addArgumentsParser(self, parser, sysargs=''):
        parser.add_argument('-t', default = self.CMAKE_BUILD_TYPES[0], help='build type', choices=self.CMAKE_BUILD_TYPES)
        parser.add_argument('--cleanlib', default = False, action="store_true", help='clean library (build)')
        parser.add_argument('--debug', default = False, action="store_true", help='debug compile')
        parser.add_argument('--verbose', default = False, action="store_true", help='compile blabbing')
        parser.add_argument('--nopy', default = False, action="store_true", help='without python wrapping')
        args = vars(parser.parse_args(sysargs))
        self.build_type = args['t']
        self.cleanlib = args['cleanlib']
        self.debug = args['debug']
        self.withoutpy = args['nopy']
        # print("args['nopy']", args['nopy'])
        if not args['verbose'] : self.runoptions="silent"
        if self.debug:
            print ('SimFemCompile addArgumentsParser() args', args)

#------------------------------------------------------------------------
    def externalInstall(self, packages):
        externalsdir = os.path.join(self.simfemsourcedir, 'External')
        if len(packages)==0:
            import glob
            items = [os.path.basename(x).split("_install.py")[0] for x in glob.glob(externalsdir+os.path.sep+"*_install.py")]
            print ('available packages in %s:\n\t%s' %(externalsdir, '\n\t'.join(items)))
            sys.exit(1)
        import importlib
        sys.path.insert(0,externalsdir)
        for package in packages:
            packageinstall = package+'_install'
            print("package", packageinstall)
            print("externalsdir", externalsdir)
            try:
                module = importlib.import_module(packageinstall,externalsdir)
            except:
                raise ValueError("could not find python script for %s" %packageinstall)
            module.install(externalsdir=externalsdir, installdir=self.installdir, builddir=self.builddir)

#------------------------------------------------------------------------
    def compilelib(self, sourcedir):
        assert self.installdir and self.builddir
        builddir = os.path.join(self.builddir, os.path.basename(sourcedir))
        builddir = os.path.join(builddir, self.build_type)
        if self.debug:
            print ('compile: sourcedir', sourcedir)
            print ('compile: self.installdir', self.installdir)
            print ('compile: self.builddir', self.builddir)
            print ('compile: builddir', builddir)
        if self.cleanlib:
            shutil.rmtree(builddir, ignore_errors=True)
        try:
            os.makedirs(builddir)
        except:
            pass
        startdir = os.getcwd()
        os.chdir(builddir)
        cmakeoptions = " -DCMAKE_BUILD_TYPE="+self.build_type + " -DCMAKE_INSTALL_PREFIX="+self.installdir
        # print("self.withoutpy", self.withoutpy)
        if self.withoutpy : cmakeoptions = cmakeoptions +  " -DWITHOUTPYTHON=TRUE"
        else : cmakeoptions = cmakeoptions +  " -DWITHOUTPYTHON=FALSE"
        command = "cmake " + sourcedir + cmakeoptions
        returncode = runsubprocess.run(command, options=self.runoptions)
        command = "make -j4"
        returncode = runsubprocess.run(command, options=self.runoptions)
        command = "make install"
        returncode = runsubprocess.run(command, options=self.runoptions)
        os.chdir(startdir)
#------------------------------------------------------------------------
    def compiletest(self, sourcedir, localbuilddir="build"):
        assert self.simfemsourcedir
        assert self.simfeminstalldir
        if localbuilddir is None:
            assert self.builddir
            localbuilddir = self.builddir
        if self.debug:
            print ('compile: sourcedir', sourcedir)
            print ('compile: localinstalldir', localinstalldir)
            print ('compile: localbuilddir', localbuilddir)
        builddir = os.path.join(localbuilddir, self.build_type)
        startdir = os.getcwd()
        if not os.path.isdir(builddir):
            os.makedirs(builddir)
            os.chdir(builddir)
            command = "cmake " + sourcedir + " -DSIMFEM_INSTALL_DIR="+self.simfeminstalldir+ " -DCMAKE_BUILD_TYPE="+self.build_type
            returncode = runsubprocess.run(command, options=self.runoptions)
        else:
            os.chdir(builddir)
        command = "make -j4"
        returncode = runsubprocess.run(command, options=self.runoptions)
        command = "make install"
        returncode = runsubprocess.run(command, options=self.runoptions)
        os.chdir(startdir)
#------------------------------------------------------------------------
    def replaceinfile(self, infilename, outfilename, wordstochange):
        file = open(infilename, 'r')
        oldlines = file.readlines()
        newlines = []
        for line in oldlines:
            # for ow, nw in wordstochange.iteritems():
            for ow, nw in wordstochange.items():
                # print 'ow, nw', ow, nw
                line = line.replace(ow,nw)
                # print 'line', line
            newlines.append(line)
        file.close()
        file = open(outfilename, 'w')
        for line in newlines:
            file.write(line)
        file.close()

#------------------------------------------------------------------------
    def generateDoc(self):
        sourcedir = os.path.join(self.simfemsourcedir, 'lib')
        installdir = os.path.join(self.installdir, 'doc')
        sourcedirdoxygen = os.path.join(self.simfemsourcedir, 'doc')
        builddirdoxygen = os.path.join(self.builddir, 'doc')
        try:
            os.makedirs(builddirdoxygen)
        except:
            pass
        wordstochange={}
        wordstochange["INPUT                  ="] = "INPUT                  =%s" %(sourcedir)
        wordstochange["OUTPUT_DIRECTORY       ="] = "OUTPUT_DIRECTORY       =%s" %(installdir)
        infilename = os.path.join(sourcedirdoxygen, 'Doxyfile.in')
        outfilename = os.path.join(builddirdoxygen, 'Doxyfile')
        self.replaceinfile(infilename, outfilename, wordstochange)
        startdir = os.getcwd()
        os.chdir(builddirdoxygen)
        command = "doxygen Doxyfile"
        returncode = runsubprocess.run(command, options=self.runoptions)
        os.chdir(startdir)
#------------------------------------------------------------------------
    def _cmakeHeader(self, testname):
        cmaketext = "\
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)\n\
SET(CMAKE_MACOSX_RPATH 1)\n\
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)\n\
SET(SIMFEM_SRC_DIR %s)\n\
SET(SIMFEM_INSTALL_DIR %s)\n\
SET(PROJECT_NAME %s)\n\
SET(SIMFEMCPPLIBRARY ${SIMFEM_INSTALL_DIR}/lib/libSimFem${CMAKE_BUILD_TYPE}.dylib)\n\
PROJECT(${PROJECT_NAME})\n\
SET(CMAKE_MODULE_PATH ${SIMFEM_SRC_DIR}/CMakeModules)\n\
INCLUDE(SimFemOptions)\n\
INCLUDE(SimFemLoadVariables)\n\
INCLUDE_DIRECTORIES(${SIMFEM_SRC_DIR}/libcpp ${SIMFEM_INCLUDES})\n" % (self.simfemsourcedir, self.simfeminstalldir, testname)
        return cmaketext
#------------------------------------------------------------------------
    def writeCmakeTest(self, testname):
        outputfile = open("CMakeLists.txt",'w')
        outputfile.write(self._cmakeHeader(testname))
        cmaketext = "\
ADD_EXECUTABLE(${PROJECT_NAME} test)\n\
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${SIMFEM_EXTERNAL_LINKLIBS} ${SIMFEMCPPLIBRARY})\n\
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}\n\
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)"
        outputfile.write(cmaketext)
        outputfile.close()
        return
#------------------------------------------------------------------------
    def writeCmakeApplication(self, name, withcppy, cppsrc="cpp", cppysrc="cppy"):
        outputfile = open("CMakeLists.txt",'w')
        outputfile.write(self._cmakeHeader(name))
        cmaketextcpp = "\
SET(CPPSRC %s)\n\
AUX_SOURCE_DIRECTORY(${CPPSRC} SRCCPP)\n\
list(REMOVE_ITEM SRCCPP ${CPPSRC}/main.cpp)\n\
SET(LIBRARYCPP lib${PROJECT_NAME})\n\
ADD_LIBRARY(${LIBRARYCPP} ${SRCCPP})\n\
TARGET_LINK_LIBRARIES(${LIBRARYCPP} ${SIMFEM_EXTERNAL_LINKLIBS} ${SIMFEMCPPLIBRARY})\n\
INSTALL(TARGETS ${LIBRARYCPP} DESTINATION ${SIMFEM_INSTALL_DIR}/lib)\n\
ADD_EXECUTABLE(${PROJECT_NAME} ${CPPSRC}/main.cpp)\n\
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${SIMFEM_EXTERNAL_LINKLIBS} ${SIMFEMCPPLIBRARY} ${LIBRARYCPP})\n\
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}\n\
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)\n\
        " %(cppsrc)
        cmaketextpy = "\
SET(CPPYSRC %s)\n\
AUX_SOURCE_DIRECTORY(${CPPYSRC} SRCPY)\n\
SET(LIBRARYPY simfem${PROJECT_NAME})\n\
INCLUDE(SimFemWrapper)\n\
INCLUDE_DIRECTORIES(${CPPSRC})\n\
INCLUDE_DIRECTORIES(${SIMFEM_SRC_DIR}/lib/libcppy)\n\
ADD_LIBRARY(${LIBRARYPY} ${SRCPY})\n\
SET(SIMFEMPYLIBRARY ${SIMFEM_INSTALL_DIR}/lib/%s.so)\n\
TARGET_LINK_LIBRARIES(${LIBRARYPY} ${SIMFEM_EXTERNAL_LINKLIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${SIMFEMCPPLIBRARY} ${SIMFEMPYLIBRARY} ${LIBRARYCPP})\n\
set_target_properties(${LIBRARYPY} PROPERTIES PREFIX \"\" )\n\
set_target_properties(${LIBRARYPY} PROPERTIES SUFFIX .so)\n\
INSTALL(TARGETS ${LIBRARYPY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})\n\
        " %(cppysrc, self.LIBRARYCPPY)
        outputfile.write(cmaketextcpp)
        if withcppy and not self.withoutpy: outputfile.write(cmaketextpy)
        outputfile.close()
        return
#------------------------------------------------------------------------
    def writeCmakeProgtest(self, name, srcdir, cppysrc="."):
        outputfile = open("CMakeLists.txt",'w')
        outputfile.write(self._cmakeHeader(name))
        cmaketext = "\
SET(CPPSRC %s)\n\
SET(CPPYSRC %s)\n\
AUX_SOURCE_DIRECTORY(${CPPYSRC} SRCPY)\n\
SET(LIBRARYPY simfem${PROJECT_NAME})\n\
INCLUDE(SimFemWrapper)\n\
INCLUDE_DIRECTORIES(${SIMFEM_SRC_DIR}/lib/libcppy)\n\
# ADD_LIBRARY(${LIBRARYPY} ${SRCPY})\n\
# TARGET_LINK_LIBRARIES(${LIBRARYPY} ${SIMFEM_EXTERNAL_LINKLIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${SIMFEMCPPLIBRARY})\n\
# set_target_properties(${LIBRARYPY} PROPERTIES PREFIX \"\" )\n\
# set_target_properties(${LIBRARYPY} PROPERTIES SUFFIX .so)\n\
# INSTALL(TARGETS ${LIBRARYPY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})\n\
        " %(cppysrc,cppysrc)
        outputfile.write(cmaketext)
        srcdircmake = os.path.join(srcdir, "CMakeLists.txt")
        if os.path.isfile(srcdircmake):
            inputfile = open(srcdircmake)
            text = inputfile.read()
            outputfile.write(text)
        outputfile.close()
        return
