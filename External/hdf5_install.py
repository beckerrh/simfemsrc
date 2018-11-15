import os, shutil, subprocess
from packageinstall import PackageInstall

def install(externalsdir, installdir, builddir):
    pi = PackageInstall(os.path.join(externalsdir, "hdf5-1.10.1.tgz"), installdir, builddir)
    pi.install_configure(configureoptions="--enable-cxx")
