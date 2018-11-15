import os
from packageinstall import PackageInstall

def install(externalsdir, installdir, builddir):
    pi = PackageInstall(os.path.join(externalsdir, "SuiteSparse-5.1.0.tgz"), installdir, builddir)
    #INSTALL=
    add = "INSTALL=%s" %(pi.installdir)
    pi.install_make(add)
