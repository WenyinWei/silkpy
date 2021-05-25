
from .geometry_map import GeometryMap as _GeometryMap
class CoordTransform(_GeometryMap):
    pass


class Cylindrical2Cartesian(CoordTransform):
    def __init__(self, syms=None):
        from sympy import symbols, Array, cos, sin
        if syms is None:
            syms = symbols('R, Z, phi', real=True)
        R, Z, phi = syms[0], syms[1], syms[2]
        CoordTransform.__init__(self, [R*cos(phi), R*sin(phi), Z], syms)

class Cartesian2Cylindrical(CoordTransform):
    def __init__(self, syms=None):
        from sympy import symbols, Array, cos, sin, sqrt, atan2
        if syms is None:
            syms = symbols('x, y, z', real=True)
        x, y, z = syms[0], syms[1], syms[2]
        CoordTransform.__init__(self, [sqrt(x**2+y**2), z, atan2(x, y)], syms)


