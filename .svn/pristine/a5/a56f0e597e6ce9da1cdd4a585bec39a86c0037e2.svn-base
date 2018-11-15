
from mesh import geometry

class GeometryClass(geometry.Geometry):
    def __init__(self, **kwargs):
        kwargs['name'] = 'unitcube'
        geometry.Geometry.__init__(self,**kwargs)

    def defineGeometry(self):
        h = self.hmean
        p0 =  self.add_point([-1.0, -1.0, -1.0], h)
        p1 =  self.add_point([ 1.0, -1.0, -1.0], h)
        p2 =  self.add_point([ 1.0,  1.0, -1.0], h)
        p3 =  self.add_point([-1.0,  1.0, -1.0], h)
        p4 =  self.add_point([-1.0, -1.0,  1.0], h)
        p5 =  self.add_point([ 1.0, -1.0,  1.0], h)
        p6 =  self.add_point([ 1.0,  1.0,  1.0], h)
        p7 =  self.add_point([-1.0,  1.0,  1.0], h)
        self.add_physical_point(p0, label=11)
        self.add_physical_point(p1, label=12)
        self.add_physical_point(p2, label=13)
        self.add_physical_point(p3, label=14)
        self.add_physical_point(p4, label=15)
        self.add_physical_point(p5, label=16)
        self.add_physical_point(p6, label=17)
        self.add_physical_point(p7, label=18)
        l0 =  self.add_line(p0, p1)
        l1 =  self.add_line(p1, p2)
        l2 =  self.add_line(p2, p3)
        l3 =  self.add_line(p3, p0)
        l4 =  self.add_line(p4, p5)
        l5 =  self.add_line(p5, p6)
        l6 =  self.add_line(p6, p7)
        l7 =  self.add_line(p7, p4)
        l8 =  self.add_line(p0, p4)
        l9 =  self.add_line(p1, p5)
        l10 =  self.add_line(p2, p6)
        l11 =  self.add_line(p3, p7)
        self.add_physical_line(l0, label=111)
        self.add_physical_line(l1, label=112)
        self.add_physical_line(l2, label=113)
        self.add_physical_line(l3, label=114)
        self.add_physical_line(l4, label=115)
        self.add_physical_line(l5, label=116)
        self.add_physical_line(l6, label=117)
        self.add_physical_line(l7, label=118)
        self.add_physical_line(l8, label=119)
        self.add_physical_line(l9, label=120)
        self.add_physical_line(l10, label=121)
        self.add_physical_line(l11, label=122)
        ll1 =  self.add_line_loop([l0, l1, l2, l3])
        surf1 =  self.add_plane_surface(ll1)
        self.add_physical_surface(surf1, label=1111)
        ll2 =  self.add_line_loop([l4, l5, l6, l7])
        surf2 =  self.add_plane_surface(ll2)
        self.add_physical_surface(surf2, label=1112)
        ll3 =  self.add_line_loop([l0, l9, -l4, -l8])
        surf3 =  self.add_plane_surface(ll3)
        self.add_physical_surface(surf3, label=1113)
        ll4 =  self.add_line_loop([l1, l10, -l5, -l9])
        surf4 =  self.add_plane_surface(ll4)
        self.add_physical_surface(surf4, label=1114)
        ll5 =  self.add_line_loop([l2, l11, -l6, -l10])
        surf5 =  self.add_plane_surface(ll5)
        self.add_physical_surface(surf5, label=1115)
        ll6 =  self.add_line_loop([l3, l8, -l7, -l11])
        surf6 =  self.add_plane_surface(ll6)
        self.add_physical_surface(surf6, label=1116)
        sl = self.add_surface_loop([surf1, surf2, surf3, surf4, surf5, surf6])
        vol = self.add_volume(sl)
        self.add_physical_volume(vol, label=9999)




    # quad = self.add_polygon([
    #   [-1.0,-1.0, 0.0],
    #   [ 1.0,-1.0, 0.0],
    #   [ 1.0, 1.0, 0.0],
    #   [-1.0, 1.0, 0.0],
    #   ],
    #   h)
    # # print 'quad', vars(quad)
    # # print 'quad.surface', vars(quad.surface)
    # # print 'quad.line_loop', vars(quad.line_loop)
    # self.add_physical_line(quad.line_loop.lines[0], label=11)
    # self.add_physical_line(quad.line_loop.lines[1], label=22)
    # self.add_physical_line(quad.line_loop.lines[2], label=33)
    # self.add_physical_line(quad.line_loop.lines[3], label=44)
    # self.add_physical_surface(quad.surface, label=111)
    # axis = [0, 0, 2]
    # top, vol, ext = self.extrude(quad.surface, axis)
    # print 'vol', vars(vol)
    # print 'top', vars(top)
    # print 'top.id', top.id
    # print 'ext[0]', vars(ext[0])
    # self.add_physical_surface(top, label=666)
    # self.add_physical_surface(ext[0], label=222)
    # self.add_physical_surface(ext[1], label=333)
    # self.add_physical_surface(ext[2], label=444)
    # self.add_physical_surface(ext[3], label=555)
    # self.add_physical_volume(vol, label=9999)
