
from ..geometry_map import GeometryMap as _GeometryMap
class CoordTransform(_GeometryMap):

    def chain(self, other):
        from silkpy.curve import ParametricCurve
        from silkpy.surface import ParametricSurface

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
