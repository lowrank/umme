import netgen.gui
from netgen.geom2d import SplineGeometry
from ngsolve import *

import matplotlib.pyplot as plt
import numpy as np


class UMMEG(object):
    """
    Solve time harmonic Maxwell's equation.
    curl (mu^{-1} curl E) - C E = F     in R^2

    with PEC boundary.

    The functions E, C, B, F, g are complex-valued.
    The relative permeability mu is real-valued.

    C = eps omega^2 + i sigma omega

    For 2D, the curl curl should be understood by extending
    the z-component as zero.
    """

    def __init__(self, args):
        self.frequency    = args.frequency     # omega / speed of light
        self.permittivity = args.permittivity  # relative permittivity
        self.eps0         = args.eps0          # background relative permittivity
        self.permeability = args.permeability  # relative permeability
        self.conductivity = args.conductivity  # conductivity

        self.source       = args.source        # real part of source F
        self.mesh         = args.mesh          # geometry
        self.imp          = args.imp           # impedance

        self._set_finite_element(args.deg)
        self._set_coefficient()

    def _set_finite_element(self, deg):
        self.finite_element_space = HCurl(self.mesh, complex=True, order=deg)

    def _set_coefficient(self):
        self.C_r = self.frequency * self.frequency * self.permittivity
        self.C_c = self.frequency * self.conductivity

    def set_coefficient(self, C_r, C_c):
        self.C_r = C_r
        self.C_c = C_c

    def construct_variational_form(self):
        U = self.finite_element_space.TrialFunction()
        V = self.finite_element_space.TestFunction()

        A = BilinearForm(self.finite_element_space, symmetric=True)
        A += (1.0/self.permeability) * curl(U) * curl(V) * dx 
        A += -self.C_r * U * V * dx 
        A += -1J * self.C_c * U * V * dx 
        A += -1J *sqrt(self.eps0) *  self.frequency * U.Trace() * V.Trace() * ds

        b = LinearForm (self.finite_element_space)
        b += SymbolicLFI(1J * self.frequency * self.source * V)
        b += self.imp * V.Trace() * ds

        A.Assemble()
        b.Assemble()

        return A, b

    def construct_linear_form(self):
        U = self.finite_element_space.TrialFunction()
        V = self.finite_element_space.TestFunction() 

        b = BilinearForm(self.finite_element_space)
        b += U * V * dx

        b.Assemble()

        return b
        
    def solve(self, A, b, filename_str):
        SOL = GridFunction(self.finite_element_space)
        SOL.vec.data = A.mat.Inverse() * b.vec

        # VTKOutput object

        if len(filename_str) > 0:
            vtk = VTKOutput(ma=self.mesh,
                            coefs=[SOL.real],
                            names=["solution_real"],
                            filename=filename_str + '_real')
            # Exporting the results:
            vtk.Do()

            # VTKOutput object
            vtk = VTKOutput(ma=self.mesh,
                            coefs=[SOL.imag],
                            names=["solution_imag"],
                            filename=filename_str + '_imag')
            # Exporting the results:
            vtk.Do()

        return SOL
