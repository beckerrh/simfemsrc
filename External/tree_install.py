import os, shutil
from packageinstall import PackageInstall

def install(externalsdir, installdir, builddir):
    pi = PackageInstall(os.path.join(externalsdir, "tree-3.1.tar"), installdir, builddir)
    # print("pi.packagename", pi.packagename)
    # print("pi.builddir", pi.builddir)
    treefile = pi.builddir + '/src/tree.hh'
    includedir = installdir+'/include'
    if not os.path.exists(includedir): os.mkdir(includedir)
    includedir += '/tree'
    if not os.path.exists(includedir): os.mkdir(includedir)
    shutil.copy(treefile, includedir )
    