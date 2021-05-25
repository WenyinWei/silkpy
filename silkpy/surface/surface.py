from sympy import Array as _Array


from ..geometry_map import GeometryMap
class ParametricSurface(GeometryMap):
    
    def __init__(self, exprs, syms):
        GeometryMap.__init__(self, exprs, syms)

    # def r(self): return self.expr()
    # def u(self): return self._u
    # def v(self): return self._v
    # def u_limit(self): return self._u_limit
    # def v_limit(self): return self._v_limit
    
    # def __str__(self):
    #     return f"A surface = {self.expr()}, with {self.u} domain {self._u_limit}, {self.v} domain {self._v_limit}."
        
    def E_F_G(self):
        from silkpy.sympy_utility import dot
        r_u = self.expr().diff(self.sym(0))
        r_v = self.expr().diff(self.sym(1))
        return dot(r_u, r_u), dot(r_u, r_v), dot(r_v, r_v)
    def metric(self):
        from einsteinpy.symbolic import MetricTensor
        E, F, G = self.E_F_G()
        return MetricTensor([[E, F], [F, G]], self.sym(), config='ll')
            
    def normal(self):
        from silkpy.sympy_utility import cross, norm
        from einsteinpy.symbolic.vector import GenericVector
        r_u = self.expr().diff(self.sym(0))
        r_v = self.expr().diff(self.sym(1))
        r_u_x_r_v = cross(r_u, r_v)
        # cross product of r_u and r_v
        return _Array(r_u_x_r_v / norm(r_u_x_r_v))
    def L_M_N(self):
        from silkpy.sympy_utility import dot
        n = self.normal()
        r_uu = self.expr().diff(self.sym(0), 2)
        r_uv = self.expr().diff(self.sym(0), self.sym(1))
        r_vv = self.expr().diff(self.sym(1), 2)
        return dot(r_uu, n), dot(r_uv, n), dot(r_vv, n)
    def Omega(self, i,j):
        from einsteinpy.symbolic.tensor import BaseRelativityTensor
        L, M, N = self.L_M_N()
        return BaseRelativityTensor(
                    [[L, M], [M, N]], 
                    self.sym(), 
                    config='ll', 
                    parent_metric=self.metric() # TODO: check the metric.
                )
        
    def Christoffel(self):
        from einsteinpy.symbolic import ChristoffelSymbols
        return ChristoffelSymbols.from_metric(self.metric())
    def RiemannCurvature(self):
        from einsteinpy.symbolic import RiemannCurvatureTensor
        return RiemannCurvatureTensor.from_christoffels(self.Christoffel())

    def omega(self):
        E, F, G = self.E_F_G()
        ans = 0
        for l in range(1, 2+1):
            if (j,l) == (1,1):
                g__jl = G
            elif ((j,l) == (1,2)) or ((j,l) == (2,1)):
                g__jl =-F
            elif (j,l) == (2,2):
                g__jl = E
            g__jl = E*G-F**2
            
            ans = ans + g__jl * self.Omega(l, i)
        return ans
    
    def WeingartenMatrix(self):
        from sympy import Matrix
        W = Matrix(
            (self.omega(1, 2), self.omega(1, 2)), 
            (self.omega(2, 1), self.omega(2, 2)))
        return W
    
    def prin_curvature_dir(self):
        from silkpy.sympy_utility import norm
        W = self.WeingartenMatrix()
        r_u, r_v = self.r_u(), self.r_v()
        eigen = W.eigenvects()
        if eigen[0][1] == 2: # if the two eigenvalues are identical
            self._k1 = self._k2 = eigen[0][0]
            er1 = eigen[0][2][0] * r_u +  eigen[0][2][1] * r_v
            self._e1 = er1 / norm(er1)
            self._e2 = self._e1
        else:
            self._k1 = eigen[0][0]
            self._k2 = eigen[1][0]
            er1 = eigen[0][2][0] * r_u +  eigen[0][2][1] * r_v
            self._e1 = er1 / norm(er1)
            er2 = eigen[1][2][0] * r_u +  eigen[1][2][1] * r_v
            self._e2 = er2 / norm(er2)
        return (self._k1, self._e1), (self._k2, self._e2)
    def KH(self):
        return self._k1 * self._k2, (self._k1 + self._k2) / 2
    
    def Weingarten_transform(self, v): 
        """
        Args:
        v: planar vector in tangent plane.
        """
        from sympy.solvers.solveset import linsolve
        from sympy import Matrix, symbols
        c_ru, c_rv = symbols('c_ru, c_rv', real=True)
        r_u, r_v = self.r_u(),  self.r_v()
        solset = linsolve(Matrix(
            ((r_u[0], r_v[0], v[0]), 
             (r_u[1], r_v[1], v[1]), 
             (r_u[2], r_v[2], v[2]))), (c_ru, c_rv))
        if len(solset) != 1:
            raise RuntimeError(f"Sympy is not smart enough to decompose the v vector with r_u, r_v as the basis.\
            It found these solutions: {solset}.\
            Users need to choose from them or deduce manually, and then set it by arguments.")
        c_ru, c_rv = next(iter(solset))
        W_r_u = self.omega(1, 1) * r_u + self.omega(1, 2) * r_v
        W_r_v = self.omega(2, 1) * r_u + self.omega(2, 2) * r_v
        return c_ru * W_r_u + c_rv * W_r_v
        
def param_transform(other_surface, newu=None, newv=None, u_expr=None, v_expr=None):
    from sympy import S, solveset, Eq
    from silkpy.sympy_utility import norm

    if newu is None or newv is None: 
        from sympy import Symbol
        newu = Symbol('u', real=True)
        newv = Symbol('v', real=True)
        
        from sympy import integrate
        drdt = norm(other_surface.expr().diff(other_surface._t))
        s_expr = integrate(drdt, other_surface._t).simplify()
        solset = solveset(Eq(s, s_expr), t, domain=S.Reals)
        if len(solset) != 1:
            raise RuntimeError(f"Sympy is not smart enough to inverse s(t) into t(s).\
            It found these solutions: {solset}.\
            Users need to choose from them or deduce manually, and then set it by obj.param_norm(s_symbol, t_expressed_by_s")
        t_expr = next(iter(solset))
    return Surface(
        other_surface.expr().applyfunc(lambda x: x.subs(other_surface._u, u_expr)).applyfunc(lambda x: x.subs(other_surface._v, v_expr)), 
        (newu, newu_expr.subs(other_surface._u, other_surface._u_limit[0]), newu_expr.subs(other_surface._u, other_surface._u_limit[1])),
        (newv, newv_expr.subs(other_surface._v, other_surface._v_limit[0]), newv_expr.subs(other_surface._v, other_surface._v_limit[1])), other_surface._sys)