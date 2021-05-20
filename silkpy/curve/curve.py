from sympy import Array as _Array
from sympy import diff as _diff
from sympy import sqrt as _sqrt


from ..geometry_map import GeometryMap
class ParametricCurve(GeometryMap):
    def __init__(self, exprs, syms):
        
        GeometryMap.__init__(self, exprs, syms)

#     def r_t(self):
#         return self._r.diff(self._t)
    
    # def __str__(self):
    #     return f"A curve = {self._r} in {self._coord} coordinate system, with {self._t} in {self._t_limit}."
    
    def curvature(self):
        drdt = self._r.diff(self._t)
        d2rdt2 = drdt.diff(self._t)
        return norm(cross(drdt, d2rdt2)) / norm(drdt)**3
        
    def torsion(self):
        drdt = self._r.diff(self._t)
        d2rdt2 = drdt.diff(self._t)
        d3rdt3 = d2rdt2.diff(self._t)
        
        numerator = triple_prod(drdt, d2rdt2, d3rdt3)
        denominator = norm( cross( drdt, d2rdt2 ) )**3
        return numerator / denominator



def param_transform(old_curve, newt=None, t_expr=None):
    from sympy import S, solveset, Eq
    if newt is None: 
        from sympy import Symbol
        newt = Symbol('s', real=True)
        
        from sympy import integrate
        drdt = norm(self.r_t())
        newt_expr = integrate(drdt, self._t).simplify()
        solset = solveset(Eq(newt, newt_expr), t, domain=S.Reals)
        if len(solset) != 1:
            raise RuntimeError(f"Sympy is not smart enough to inverse s(t) into t(s).\
            It found these solutions: {solset}.\
            Users need to choose from them or deduce manually, and then set it by obj.param_norm(s_symbol, t_expressed_by_s")
        t_expr = next(iter(solset))
    return ParametricCurve(
        self._r.applyfunc(lambda x: x.subs(self._t, t_expr)), 
        (newt, newt_expr.subs(self._t, self._t_limit[0]), newt_expr.subs(self._t, self._t_limit[1])), self._sys)