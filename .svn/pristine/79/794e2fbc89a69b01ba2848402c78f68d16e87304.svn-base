
#------------------------------------------------------------------------
def testIo(mesh):
    import tempfile
    import subprocess, os, shutil
    dirname = tempfile.mkdtemp()
    # print 'dirname', dirname
    th1 = dirname + os.path.sep + 'h1'
    th2 = dirname + os.path.sep + 'h2'
    ta1 = dirname + os.path.sep + 'a1'
    ta2 = dirname + os.path.sep + 'a2'
    mesh.save(th1)
    mesh.load(th1)
    mesh.save(th2)
    subprocess.check_call("h5dump " + th1 + " >> " + ta1, shell=True)
    subprocess.check_call("h5dump " + th2 + " >> " + ta2, shell=True)
    cmd = "diff -C 0 %s %s"%(ta1, ta2)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    out, err = process.communicate()
    # print(out)
    skips = [ta1, ta2, th1, th2, '*** 1 ****', '--- 1 ----','***************']
    ok = True
    for line in out.strip().split('\n'):
        conditions = [skip in line for skip in skips]
        if any(conditions): continue
        ok = False
        print "line ", line
    shutil.rmtree(dirname)
    return ok


def showVtk(mesh, name="meshtest"):
    import vtk
    import mesh.simfemvtk
    vtkname = name+'.vtk'
    vtkbdryname = name+'_boundary.vtk'
    mesh.writeVtk(vtkname)
    mesh.writeBoundaryVtk(vtkbdryname)
    simfemvtk = mesh.simfemvtk.SimFemVtk()
    simfemvtk.addMeshActor(vtkname)
    simfemvtk.addMeshActor(vtkbdryname)
    simfemvtk.addAxes()
    simfemvtk.run()
