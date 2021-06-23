import netgen.gui
from netgen.geom2d import SplineGeometry
from ngsolve import *

import matplotlib.pyplot as plt
import numpy as np


class UMMEG(object):
    """
    Solve time harmonic Maxwell's equation.
    curl (mu^{-1} curl E) - C E = F     in R^2

    with Silver-Muller boundary condition.

    The functions E, C, B, F, g are complex-valued.
    The relative permeability mu is real-valued.

    C = eps omega^2 + i sigma omega

    For 2D, the curl curl should be understood by extending
    the z-component as zero.
    """

    def __init__(self, args):
        self.frequency    = args.frequency     # omega / speed of light
        self.permittivity = args.permittivity  # relative permittivity
        self.permeability = args.permeability  # relative permeability
        self.conductivity = args.conductivity  # conductivity

        self.source       = args.source        # real part of source F
        self.mesh         = args.mesh          # geometry

        self._set_finite_element(args.deg)
        self._set_coefficient()

    def _set_finite_element(self, deg):
        self.finite_element_space = HCurl(self.mesh, complex=True, order=deg)

    def _set_coefficient(self):
        self.C_r = self.frequency * self.frequency * self.permittivity
        self.C_c = self.frequency * self.conductivity


    def construct_variational_form(self):
        U = self.finite_element_space.TrialFunction()
        V = self.finite_element_space.TestFunction()

        A = BilinearForm(self.finite_element_space, symmetric=True)
        A += SymbolicBFI((1.0/self.permeability) * curl(U) * curl(V))
        A += SymbolicBFI(-self.C_r * U * V)
        A += SymbolicBFI(-1J * self.C_c * U * V)
        A += SymbolicBFI(-1J * self.frequency * U.Trace() * V.Trace(), BND)

        b = LinearForm (self.finite_element_space)
        b += SymbolicLFI(1J * self.frequency * self.source * V)

        A.Assemble()
        b.Assemble()

        return A, b

    def solve(self, A, b, filename_str):
        SOL = GridFunction(self.finite_element_space)
        SOL.vec.data = A.mat.Inverse() * b.vec

        # VTKOutput object
        vtk = VTKOutput(ma=self.mesh,
                        coefs=[SOL.real, SOL.imag],
                        names=["solution"],
                        filename=filename_str,
                        subdivision=3)
        # Exporting the results:
        vtk.Do()

        return SOL