from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from matplotlib import interactive
import numpy as np
class MEG(object):
    """
    Solve time harmonic Maxwell's equation. 
    curl (mu^{-1} curl E) - C E = F     in \Omega,
                         n x E  = 0     on dbc,
            mu^{-1}curl E - B E = g     on ibc.

    dbc, ibc refer Dirichlet boundary (PEC) and Impedance
    boundary (1st order approximation to ABC). 

    The functions E, C, B, F, g are complex-valued.
    The relative permeability mu is real-valued.

    B = i omega 
    C = eps_r omega^2 + i (sigma/epsilon_0) omega

    For 2D, the curl curl should be understood by extending
    the z-component as zero.
    """
    def __init__(self, args):
        self.frequency    = args.frequency            # frequency number
        self.permittivity = args.permittivity         # relative permittivity
        self.permeability = args.permeability         # relative permeability
        self.conductivity = args.conductivity         # conductivity
        
        self.source_r     = args.source_r             # real part of source F
        self.source_c     = args.source_c             # imag part of source F
        self.imp_r        = args.imp_r                # real part of impedance g
        self.imp_c        = args.imp_c                # imag part of impedance g


        self.Lambda       = args.Lambda               # impedance parameter

        self.mesh         = args.mesh                 # geometry 
        self.dim          = args.dim                  # dimension of mesh

        self.dbc          = args.dbc                  # Dirichlet boundary func
        
        self.set_finite_element(args.scalar_elem, args.vector_elem, args.deg)
        self.set_coefficient()
        self.set_boundary(args.ibc)

    def set_finite_element(self, scalar_elem, vector_elem, deg):
        self.s_elem = FiniteElement(
            scalar_elem, cell=self.mesh.ufl_cell(), degree=deg)

        self.scalar_space = FunctionSpace(self.mesh, self.s_elem)

        self.v_elem = FiniteElement(
            vector_elem, cell=self.mesh.ufl_cell(), degree=deg)

        self.vector_space = VectorFunctionSpace(self.mesh, scalar_elem, degree=deg)

        mix_elem = MixedElement([self.v_elem, self.v_elem])

        self.tensor_space = FunctionSpace(self.mesh, mix_elem)

        self.n = FacetNormal(self.mesh)

    def set_coefficient(self):
        # The coefficients are 'ufl.algebra.Product'
        self.C_r = self.frequency * self.frequency * self.permittivity
        self.C_c = self.frequency * self.conductivity

        self.B_r = 0.0
        self.B_c = self.Lambda * self.frequency 

    def set_boundary(self, _ibc):
        boundary_markers = MeshFunction('size_t', self.mesh, dim = 1)

        class ibc(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary  and _ibc(x)

        _Impedance = ibc()
        _Impedance.mark(boundary_markers, 0)

        self.ds = Measure('ds', \
                          domain=self.mesh, \
                          subdomain_data=boundary_markers)

    def extend(self, v):
        # In order to use 'cross', in 3D, this is not necessary.
        if self.dim == 2:
            return as_vector( (v[0], v[1], 0.) )
        elif self.dim == 3:
            return v
        else:
            raise ValueError('dimension can only be 2 or 3.')

    def trace(self, v):
        return cross(\
            self.extend(self.n), \
            cross(self.extend(v), self.extend(self.n))
            )

    def construct_variational_form(self):
        # real, complex components
        E_r, E_c = TrialFunctions(self.tensor_space)
        P,   Q   = TestFunctions (self.tensor_space)

        A_r = ((1.0/self.permeability)     \
            * (inner(curl(E_r), curl(P)))  \
            - self.C_r * (inner(E_r, P))   \
            + self.C_c * (inner(E_c, P))) * dx  \
            + (\
                self.B_c * inner(self.trace(E_c), self.trace(P)) - \
                self.B_r * inner(self.trace(E_r), self.trace(P))   \
                ) * self.ds(0)

        A_c = ((1.0/self.permeability)      \
            * (inner(curl(E_c), curl(Q)))   \
            - self.C_r * (inner(E_c, Q))    \
            - self.C_c * (inner(E_r, Q)) ) * dx  \
            - (\
                self.B_r * inner(self.trace(E_c), self.trace(Q)) + \
                self.B_c * inner(self.trace(E_r), self.trace(Q))   \
                ) * self.ds(0)

        A = A_r + A_c

        b_r = inner(self.source_r, P) * dx + \
              inner(self.trace(P), self.extend(self.imp_r) \
                    ) * self.ds(0)
        
        b_c = inner(self.source_c, Q) * dx + \
              inner(self.trace(Q), self.extend(self.imp_c) \
                    ) * self.ds(0)

        b = b_r + b_c
        
        return A, b

    def solve(self, A, b):
        u = Function(self.tensor_space)

        if self.dbc == None:
            bcs = []
        else:
            bcs = [ DirichletBC(self.tensor_space.sub(0), \
                            Constant((0., 0.)), \
                            self.dbc), \
                    DirichletBC(self.tensor_space.sub(1),
                            Constant((0., 0.)),
                            self.dbc) \
                ]

        solve(A == b, u, bcs)

        return u.split()

    def visualize(self, _u, _v):
        # visualization only for 2D.
        if self.dim == 3:
            raise ValueError('3 dimension not implemented.')
        else:
            u = interpolate(_u, self.vector_space)
            v = interpolate(_v, self.vector_space)

            u_ax0 = plt.subplot(221)
            u_c0  = plot(u.sub(0), mode='color', cmap='jet')
            u_ax0.set_title('x component of real part')
            u_cbar0 = plt.colorbar(u_c0)
            u_cbar0.formatter.set_powerlimits((0, 0))

            u_ax1 = plt.subplot(222)
            u_c1  = plot(u.sub(1), mode='color', cmap='jet')
            u_ax1.set_title('y component of real part')
            u_cbar1 = plt.colorbar(u_c1)
            u_cbar1.formatter.set_powerlimits((0, 0))

            v_ax0 = plt.subplot(223)
            v_c0  = plot(v.sub(0), mode='color', cmap='jet')
            v_ax0.set_title('x component of imag part')
            v_cbar0 = plt.colorbar(v_c0)
            v_cbar0.formatter.set_powerlimits((0, 0))

            v_ax1 = plt.subplot(224)
            v_c1 = plot(v.sub(1), mode='color', cmap='jet')
            v_ax1.set_title('x component of real part')
            v_cbar1 = plt.colorbar(v_c1)
            v_cbar1.formatter.set_powerlimits((0, 0))




