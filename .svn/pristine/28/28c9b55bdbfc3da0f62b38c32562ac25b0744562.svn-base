
from mesh import geometry

class GeometryClass(geometry.Geometry):
  def __init__(self, **kwargs):
    kwargs['name'] = 'unitline'
    kwargs['dim'] = 1
    geometry.Geometry.__init__(self,**kwargs)

  def defineGeometry(self):
    h = self.hmean
    p0 =  self.add_point([-1.0, 0.0, 0.0], h)
    p1 =  self.add_point([1.0, 0.0, 0.0], h)
    l0 =  self.add_line(p0, p1)
    self.add_physical_point(p0, label=11)
    self.add_physical_point(p1, label=22)
    self.add_physical_line(l0, label=99)
