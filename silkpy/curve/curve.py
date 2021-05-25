from sympy import Array as _Array
from sympy import diff as _diff
from sympy import sqrt as _sqrt


from ..geometry_map import GeometryMap
class ParametricCurve(GeometryMap):
    def __init__(self, exprs, sym):
        if not isinstance(sym, list): sym = [sym]
        GeometryMap.__init__(self, exprs, sym)
        

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

    def unit_tangent_vec(self): # \vec{T}(t)
        from silkpy.sympy_utility import norm
        drdt = self.expr().diff(self.sym(0))
        return (drdt / norm(drdt)).simplify()
    def unit_normal_vec(self): # \vec{N}(t)
        from silkpy.sympy_utility import norm
        d2rdt2 = self.expr().diff(self.sym(0), 2)
        return (d2rdt2 / norm(d2rdt2)).simplify()
    def unit_subnormal_vec(self): # \vec{B}(t)
        from silkpy.sympy_utility import cross
        return cross( self.unit_tangent_vec(), self.unit_normal_vec() ).simplify()



