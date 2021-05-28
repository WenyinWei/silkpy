from sympy import Array as _Array
from sympy import diff as _diff
from sympy import sqrt as _sqrt


from ..geometry_map import GeometryMap
class ParametricCurve(GeometryMap):
    def __init__(self, exprs, sym):
        if not isinstance(sym, list): sym = [sym]
        GeometryMap.__init__(self, exprs, sym)
    def subs(self, subs_arg):
        return ParametricCurve(self._exprs.subs(subs_arg), self._syms)

    def chain(self, other):
        from ..geometry_map.coord_transform import CoordTransform

        assert( len(other.sym()) == int(self.expr().shape.args[0]) )
        
        if isinstance(other, CoordTransform):
            return ParametricCurve(other._exprs.subs({
                other.sym(i): self.expr(i) for i in range(len(self.sym()))}), 
                self._syms)
        else:
            raise TypeError("The chain succession must be a GeometryMap.")
    def __or__(self, other):
        return self.chain(other)


#     def r_t(self):
#         return self._r.diff(self.sym(0))
    
    # def __str__(self):
    #     return f"A curve = {self._r} in {self._coord} coordinate system, with {self.sym(0)} in {self.sym(0)_limit}."
    
    def curvature(self): # k(t)
        from silkpy.sympy_utility import norm, cross 
        drdt = self.expr().diff(self.sym(0))
        d2rdt2 = drdt.diff(self.sym(0))
        return norm(cross(drdt, d2rdt2)) / norm(drdt)**3
        
    def torsion(self): # \tau(t)
        from silkpy.sympy_utility import norm, cross, triple_prod
        drdt = self.expr().diff(self.sym(0))
        d2rdt2 = drdt.diff(self.sym(0))
        d3rdt3 = d2rdt2.diff(self.sym(0))
        
        numerator = triple_prod(drdt, d2rdt2, d3rdt3)
        denominator = norm( cross( drdt, d2rdt2 ) )**3
        return numerator / denominator

    def T_N_B(self): # \vec{T}(t), \vec{N}(t), \vec{B}(t)
        from silkpy.sympy_utility import norm, cross
        drdt = self.expr().diff(self.sym(0))
        d2rdt2 = self.expr().diff(self.sym(0), 2)

        T = (drdt / norm(drdt)).simplify()
        N = (d2rdt2 / norm(d2rdt2)).simplify()
        B = cross(T, N).simplify()
        return T, N, B
                


