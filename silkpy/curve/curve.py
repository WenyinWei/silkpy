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
    
    def curvature(self):
        from silkpy.sympy_utility import norm, cross 
        drdt = self.expr().diff(self.sym(0))
        d2rdt2 = drdt.diff(self.sym(0))
        return norm(cross(drdt, d2rdt2)) / norm(drdt)**3
        
    def torsion(self):
        from silkpy.sympy_utility import norm, cross, triple_prod
        drdt = self.expr().diff(self.sym(0))
        d2rdt2 = drdt.diff(self.sym(0))
        d3rdt3 = d2rdt2.diff(self.sym(0))
        
        numerator = triple_prod(drdt, d2rdt2, d3rdt3)
        denominator = norm( cross( drdt, d2rdt2 ) )**3
        return numerator / denominator



def param_transform(old_curve, newt=None, t_expr=None):
    from sympy import S, solveset, Eq
    if newt is None: 
        from sympy import Symbol
        newt = Symbol('s', real=True)
        
        from sympy import integrate
        drdt = norm(old_curve.r_t())
        newt_expr = integrate(drdt, old_curve._t).simplify()
        solset = solveset(Eq(newt, newt_expr), t, domain=S.Reals)
        if len(solset) != 1:
            raise RuntimeError(f"Sympy is not smart enough to inverse s(t) into t(s).\
            It found these solutions: {solset}.\
            Users need to choose from them or deduce manually, and then set it by obj.param_norm(s_symbol, t_expressed_by_s")
        t_expr = next(iter(solset))
    return ParametricCurve(
        old_curve._r.applyfunc(lambda x: x.subs(old_curve._t, t_expr)), 
        (newt, newt_expr.subs(old_curve._t, old_curve._t_limit[0]), newt_expr.subs(old_curve._t, old_curve._t_limit[1])), old_curve._sys)