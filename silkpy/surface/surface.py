from sympy import Array as _Array


from ..geometry_map import GeometryMap
class ParametricSurface(GeometryMap):
    def __init__(self, exprs, syms):
        from sympy import oo

        GeometryMap.__init__(self, exprs, syms)

    def r(self): return self._r
    def u(self): return self._u
    def v(self): return self._v
    def u_limit(self): return self._u_limit
    def v_limit(self): return self._v_limit
    
    def __str__(self):
        return f"A surface = {self._r}, with {self.u} domain {self._u_limit}, {self.v} domain {self._v_limit}."
    
#     def r_u(self):
#         return self._r.diff(self._u)
#     def r_v(self):
#         return self._r.diff(self._v)
#     def r_uu(self):
#         return self._r.diff(self._u).diff(self._u)
#     def r_uv(self):
#         return self._r.diff(self._u).diff(self._v)
#     def r_vv(self):
#         return self._r.diff(self._v).diff(self._v)
#     def n(self):
#         r_u_x_r_v = cross(self.r_u(), self.r_v())
#         return r_u_x_r_v / norm(r_u_x_r_v)
    
    def EFG(self):
        r_u = self.r().diff(self.u())
        r_v = self.r_v()
        return dot(r_u, r_u), dot(r_u, r_v), dot(r_v, r_v)
    def g(i,j):
        E, F, G = self.EFG()
        if (i,j) == (1,1):
            return E
        elif ((i,j) == (1,2)) or ((i,j) == (2,1)):
            return F
        elif (i,j) == (2,2):
            return G
        else:
            raise ValueError("please input i, j in \{1, 2\}.")
            
    def LMN(self):
        n = self.n()
        return dot(self.r_uu(), n), dot(self.r_uv(), n), dot(self.r_vv(), n)
    def Omega(self, i,j):
        L, M, N = self.LMN()
        if (i,j) == (1,1):
            return L
        elif ((i,j) == (1,2)) or ((i,j) == (2,1)):
            return M
        elif (i,j) == (2,2):
            return N
        else:
            raise ValueError("please input i, j in \{1, 2\}.")
        
    def Christoffel(self, k, i, j):
        ans = 0
        E, F, G = self.EFG()
        for m in range(1, 2+1):
            if (k,m) == (1,1):
                g__km = G
            elif ((k,m) == (1,2)) or ((k,m) == (2,1)):
                g__km =-F
            elif (k,m) == (2,2):
                g__km = E
            g__km = E*G-F**2
            ans = ans + 1/2 * g__km * (g(m, j).diff(i) + g(i, m).diff(j) - g(i, j).diff(m))
        return ans
    
    def omega(self, i, j):
        E, F, G = self.EFG()
        ans = 0
        for l in range(1, 2+1):
            if (j,l) == (1,1):
                g__jl = G
            elif ((j,l) == (1,2)) or ((j,l) == (2,1)):
                g__jl =-F
            elif (j,l) == (2,2):
                g__jl = E
            g__jl = E*G-F**2
            
            ans = ans + g__jl * Omega(l, i)
        return ans
    
    def WeingartenMatrix(self):
        from sympy import Matrix
        W = Matrix(
            (omega(1, 2), omega(1, 2)), 
            (omega(2, 1), omega(2, 2)))
        return W
    
    def prin_curvature_dir(self):
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
        
    def param_transform(self, newu=None, newv=None, u_expr=None, v_expr=None):
        from sympy import S, solveset, Eq
        if newu is None or newv is None: 
            from sympy import Symbol
            newu = Symbol('u', real=True)
            newv = Symbol('v', real=True)
            
            from sympy import integrate
            drdt = norm(self._r.diff(self._t))
            s_expr = integrate(drdt, self._t).simplify()
            solset = solveset(Eq(s, s_expr), t, domain=S.Reals)
            if len(solset) != 1:
                raise RuntimeError(f"Sympy is not smart enough to inverse s(t) into t(s).\
                It found these solutions: {solset}.\
                Users need to choose from them or deduce manually, and then set it by obj.param_norm(s_symbol, t_expressed_by_s")
            t_expr = next(iter(solset))
        return surface(
            self._r.applyfunc(lambda x: x.subs(self._u, u_expr)).applyfunc(lambda x: x.subs(self._v, v_expr)), 
            (newu, newu_expr.subs(self._u, self._u_limit[0]), newu_expr.subs(self._u, self._u_limit[1])),
            (newv, newv_expr.subs(self._v, self._v_limit[0]), newv_expr.subs(self._v, self._v_limit[1])), self._sys)