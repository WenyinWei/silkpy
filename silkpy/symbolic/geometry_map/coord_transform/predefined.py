
from .coord_transform import CoordTransform as _CoordTransform
class Cylindrical2Cartesian(_CoordTransform):
    def __init__(self, syms=None):
        from sympy import symbols, Array, cos, sin
        if syms is None:
            _R = symbols('R', positive=True)
            _Z, _phi = symbols('Z, phi', real=True)
            syms = (_R, _Z, _phi)
        R, Z, phi = syms[0], syms[1], syms[2]
        _CoordTransform.__init__(self, [R*cos(phi), R*sin(phi), Z], syms)

class Cartesian2Cylindrical(_CoordTransform):
    def __init__(self, syms=None):
        from sympy import symbols, Array, cos, sin, sqrt, atan2
        if syms is None:
            syms = symbols('x, y, z', real=True)
        x, y, z = syms[0], syms[1], syms[2]
        _CoordTransform.__init__(self, [sqrt(x**2+y**2), z, atan2(y, x)], syms)


