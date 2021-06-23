import netgen.gui
from netgen.geom2d import SplineGeometry
from ngsolve import *
import matplotlib.pyplot as plt

import numpy as np

# Construct the domain and PML.
geo = SplineGeometry()
geo.AddCircle ( (0, 0), r=1.3, leftdomain=1, rightdomain=2, bc="domainBC")
geo.AddCircle ( (0, 0), r=1.6, leftdomain=2, bc="pmlBC")
geo.SetMaterial(1, "domain")
geo.SetMaterial(2, "pml")

ngmesh = geo.GenerateMesh(maxh=0.1)
mesh = Mesh(ngmesh)
mesh.Curve(3)

mesh.SetPML(pml.Radial(origin=(0,0), rad=1.3, alpha=0.1j), definedon=2)

mu0  = 4*np.pi*1e-7
eps0 = 8.854e-12

mu = CoefficientFunction(mu0)

freq =1E9 #1000 MHz
w    =2*np.pi*freq
k0   =w*sqrt(mu0*eps0)
print("k0=",k0)


fes = HCurl(mesh, complex=True, order=4)
u = fes.TrialFunction()
v = fes.TestFunction()

uin=CoefficientFunction( (0, IfPos(1-2*x, 1, 0) * IfPos(2*x+1, 1, 0) * IfPos(1-2*y, 1, 0) * IfPos(1+2*y,1,0) ) )
uscat = GridFunction (fes)
uscat.Set (uin, definedon=mesh.Materials("domain"))

a = BilinearForm(fes, symmetric=True)
a += SymbolicBFI(curl(u)*curl(v) - k0*k0*u*v)
a += SymbolicBFI(-1J * k0 * u.Trace() * v.Trace(), BND)

f = LinearForm (fes)
f += SymbolicLFI(uscat * v)

a.Assemble()
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse() * f.vec

Draw (gfu, mesh, "gfu")



