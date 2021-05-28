
from .geometry_map import GeometryMap as _GeometryMap
class CoordTransform(_GeometryMap):

    def chain(self, other):
        from ..curve import ParametricCurve
        from ..surface import ParametricSurface

        # assert( self.expr().rank() == 1)
        assert( len(other.sym()) == int(self.expr().shape.args[0]) )

        if isinstance(other, ParametricCurve):
            return ParametricCurve(other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))}), 
                self._syms)
        elif isinstance(other, ParametricSurface):
            return ParametricSurface(other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))}), 
                self._syms)
        elif isinstance(other, CoordTransform):
            return CoordTransform(other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))}), 
                self._syms)
        else:
            raise TypeError("The chain succession must be a GeometryMap.")
    def __or__(self, other):
        return self.chain(other)

class Cylindrical2Cartesian(CoordTransform):
    def __init__(self, syms=None):
        from sympy import symbols, Array, cos, sin
        if syms is None:
            _R = symbols('R', positive=True)
            _Z, _phi = symbols('Z, phi', real=True)
            syms = (_R, _Z, _phi)
        R, Z, phi = syms[0], syms[1], syms[2]
        CoordTransform.__init__(self, [R*cos(phi), R*sin(phi), Z], syms)

class Cartesian2Cylindrical(CoordTransform):
    def __init__(self, syms=None):
        from sympy import symbols, Array, cos, sin, sqrt, atan2
        if syms is None:
            syms = symbols('x, y, z', real=True)
        x, y, z = syms[0], syms[1], syms[2]
        CoordTransform.__init__(self, [sqrt(x**2+y**2), z, atan2(y, x)], syms)


