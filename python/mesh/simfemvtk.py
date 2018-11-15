import vtk

#------------------------------------------------------------------------
class SimFemVtk(object):
    def __init__(self):
        windowsize = 1200
        self.ren= vtk.vtkRenderer()
        self.ren.SetBackground( 0.9, 0.92, 0.94 )
        self.renwin = vtk.vtkRenderWindow()
        self.renwin.AddRenderer(self.ren)
        self.renwin.SetSize( windowsize, windowsize )
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.renwin)
        self.actors = []
#------------------------------------------------------------------------
    def run(self):
        self.iren.Initialize()
        self.iren.Start()
#------------------------------------------------------------------------
    def addActor(self, actor):
        self.actors.append(actor)
        self.ren.AddActor(actor)
#------------------------------------------------------------------------
    def addMeshActor(self, filename, wireframe=True):
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllScalarsOn()
        reader.Update()
        vtkdata = reader.GetOutput()
        scalar_range = vtkdata.GetScalarRange()
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(vtkdata)
        mapper.SetScalarRange(scalar_range)
        meshActor = vtk.vtkActor()
        meshActor.SetMapper(mapper)
        if wireframe: meshActor.GetProperty().SetRepresentationToWireframe()
        meshActor.GetProperty().SetLineWidth(3);
        self.addActor(meshActor)
#------------------------------------------------------------------------
    def addAxes(self):
        axes = vtk.vtkAxesActor()
#  The axes are positioned with a user transform
        transform = vtk.vtkTransform()
        transform.Translate(1.0, 0.0, 0.0)
        axes = vtk.vtkAxesActor()
        axes.SetUserTransform(transform)
        axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,0,0);
        axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,1,0);
        axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,1);
        # axes.SetXAxisLabelText("test");
        self.addActor(axes)
