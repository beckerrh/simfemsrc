import sys, os, shutil, subprocess

class PackageInstall(object):
    def __init__(self, sourcepath, installdir, builddir):
        self.sourcepath = sourcepath
        self.builddir = os.path.join(builddir, 'External')
        try: os.makedirs(self.builddir)
        except : pass
        self.installdir = os.path.join(installdir, 'External')
        try: os.makedirs(self.installdir)
        except : pass
        os.chdir(self.builddir)

        command = "tar tzf " + sourcepath
        proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, bufsize=1)
        with proc.stdout:
            for line in iter(proc.stdout.readline, b''):
                self.packagename = line.split('\n')[0].split('/')[0].strip('/')
        print 'self.packagename', self.packagename
        returncode = proc.wait()
        if returncode:
          raise ValueError("could not determine packagename from " + sourcepath)

        command = "tar xf " + sourcepath
        subprocess.call(command, shell=True)
        self.builddir = os.path.join(self.builddir, self.packagename)

        self.builddircomp = os.path.join(self.builddir, self.packagename+'.comp')
        if(os.path.isdir(self.builddircomp)) : shutil.rmtree(self.builddircomp)
        os.makedirs(self.builddircomp)

    def install_make(self, add):
        print("self.builddir=", self.builddir)
        os.chdir(self.builddir)
        command = "make -j4 install " + add
        subprocess.call(command, shell=True)

    def install_cmake(self, cmakeoptions=""):
        os.chdir(self.builddircomp)
        # command = "cmake -DCMAKE_INSTALL_PREFIX=%s %s ../%s" %(self.installdir, cmakeoptions, self.packagename)
        command = "cmake -DCMAKE_INSTALL_PREFIX=%s %s .." %(self.installdir, cmakeoptions)
        subprocess.call(command, shell=True)
        command = "make -j4"
        subprocess.call(command, shell=True)
        command = "make install"
        subprocess.call(command, shell=True)

    def install_configure(self, configureoptions):
        print 'configureoptions', configureoptions
        print 'self.builddir', self.builddir
        print 'os.getcwd()', os.getcwd()
        os.chdir(self.builddir)
        command = "./configure --prefix=%s %s" %(self.installdir, configureoptions)
        print 'command', command
        subprocess.call(command, shell=True)
        command = "make -j4 install"
        print 'command', command
        subprocess.call(command, shell=True)
        # shutil.copy('lib'+self.name+'.a', self.install_dir+'/lib' )
