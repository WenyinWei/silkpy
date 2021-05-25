from .curve import ParametricCurve as _ParametricCurve
from sympy import Symbol as _Symbol

def curve_normalization(
    other:_ParametricCurve,
    new_var=_Symbol('s', real=True)):
    from sympy import S, solveset, Eq
    from sympy import integrate
    from silkpy.sympy_utility import norm
    
    drdt = norm(other.expr().diff(other.sym(0)))
    new_var_in_old = integrate(drdt, other.sym(0)).simplify()
    solset = solveset(Eq(new_var, new_var_in_old), other.sym(0), domain=S.Reals).simplify()
    try:
        if len(solset) != 1:
            raise RuntimeError(f"We have not yet succedded in inverse s(t) into t(s).\
            It found these solutions: {solset}.\
            Users need to choose from them or deduce manually, and then set it by obj.param_norm(s_symbol, t_expressed_by_s")
    except:
        raise RuntimeError(f"We have not yet succedded in inverse s(t) into t(s). Try the curve_param_transform function instead and set the transform relation manually.")
    else:
        old_var_in_new = next(iter(solset))
    return _ParametricCurve(
        other.expr().subs(other.sym(0), old_var_in_new), 
        (new_var, 
         new_var_in_old.subs(other.sym(0), other.sym_limit(0)[0]),
         new_var_in_old.subs(other.sym(0), other.sym_limit(0)[1])))

def curve_param_transform(old_curve, newt, t_expr=None):
    from sympy import S, solveset, Eq
    return _ParametricCurve(
        old_curve._r.applyfunc(lambda x: x.subs(old_curve._t, t_expr)), 
        (newt, newt_expr.subs(old_curve._t, old_curve._t_limit[0]), newt_expr.subs(old_curve._t, old_curve._t_limit[1])), old_curve._sys)