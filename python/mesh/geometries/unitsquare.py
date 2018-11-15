
from mesh import geometry

class GeometryClass(geometry.Geometry):
  def __init__(self, **kwargs):
    kwargs['name'] = 'unitsquare'
    kwargs['dim'] = 2
    geometry.Geometry.__init__(self,**kwargs)

  def defineGeometry(self):
    h = self.hmean
    p0 =  self.add_point([-1.0, -1.0, 0.0], h)
    p1 =  self.add_point([1.0, -1.0, 0.0], h)
    p2 =  self.add_point([1.0, 1.0, 0.0], h)
    p3 =  self.add_point([-1.0, 1.0, 0.0], h)
    l0 =  self.add_line(p0, p1)
    l1 =  self.add_line(p1, p2)
    l2 =  self.add_line(p2, p3)
    l3 =  self.add_line(p3, p0)
    ll =  self.add_line_loop([l0, l1, l2, l3])
    surf =  self.add_plane_surface(ll)
    self.add_physical_point(p0, label=11)
    self.add_physical_point(p1, label=22)
    self.add_physical_point(p2, label=33)
    self.add_physical_point(p3, label=44)
    self.add_physical_line(l0, label=111)
    self.add_physical_line(l1, label=222)
    self.add_physical_line(l2, label=333)
    self.add_physical_line(l3, label=444)
    self.add_physical_surface(surf, label=1111)
