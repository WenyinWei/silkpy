
class GeometryMap:
    """Base class for geometry map.
    """
    def __init__(self, exprs, syms) -> None:
        from sympy import Array, oo
        if not isinstance(exprs, list):
            raise ValueError("The ctor arg of GeometryMap: exprs should be a list of equations, or a list of sympy expressions.")
        self._exprs = exprs

        if not isinstance(syms, list):
            raise ValueError("The ctor arg of GeometryMap: syms should be a list of sympy variables, or (sympy.Symbol, inf, sup)")
        self._syms = []
        for sym in syms:
            if isinstance(sym, tuple):
                assert(len(sym)==3)
                self._syms.append(sym)
            else:
                self._syms.append( (sym, -oo, +oo) )

    def expr(self, i=None):
        if i is None:
            return self._exprs
        else:
            return self._exprs[i]

    def sym(self, i=None):
        if i is None:
            return [sym[0] for sym in self._syms]
        else:
            return self._syms[i][0]

    def sym_limit(self, i=None):
        if i is None:
            return [sym[1:] for sym in self._syms]
        else:
            return self._syms[i][1:]