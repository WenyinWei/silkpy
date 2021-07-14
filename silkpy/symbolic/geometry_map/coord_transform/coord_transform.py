
from ..geometry_map import GeometryMap as _GeometryMap
class CoordTransform(_GeometryMap):
    def __or__(self, other):
        from silkpy.symbolic.curve.curve import ParametricCurve
        from silkpy.symbolic.surface.surface import ParametricSurface

        # assert( self.expr().rank() == 1)
        assert( len(other.sym()) == int(self.expr().shape.args[0]) )

        if isinstance(other, ParametricCurve):
            return ParametricCurve(self._syms, other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))})
                )
        elif isinstance(other, ParametricSurface):
            return ParametricSurface(self._syms, other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))})
                )
        elif isinstance(other, CoordTransform):
            return CoordTransform(self._syms, other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))})
                )
        else:
            raise TypeError("The chain succession must be a GeometryMap.")