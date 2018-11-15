import os
from packageinstall import PackageInstall

def install(externalsdir, installdir, builddir):
  # installdir = "/usr/local"
  pi = PackageInstall(os.path.join(externalsdir, "Xdmf_13_02_2018.tgz"), installdir, builddir)
  os.environ["HDF5_ROOT"] = "/usr/local"
  #XDMF_WRAP_PYTHON
  cmakeoptions="-DBUILD_SHARED_LIBS=1 -DXDMF_WRAP_PYTHON=1 -Wno-dev -DXDMF_BUILD_TESTING=1"
  pi.install_cmake(cmakeoptions=cmakeoptions)
