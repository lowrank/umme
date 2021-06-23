from UMMEG import *

'''
In the following example, the domain is set in a unit disk and extended with PML of width 0.4. 

The maximum h of each element is set as 0.1, see the code.
'''
geometry = SplineGeometry()

geometry.AddCircle ((0, 0), r=1.0, leftdomain=1, rightdomain=2, bc="domainBC")
geometry.AddCircle ((0, 0), r=1.4, leftdomain=2, bc="pmlBC")

geometry.SetMaterial(1, "domain")
geometry.SetMaterial(2, "pml")

ngmesh = geometry.GenerateMesh(maxh=0.03)
mesh = Mesh(ngmesh)
mesh.Curve(3)

mesh.SetPML(pml.Radial(origin=(0,0), rad=1.0, alpha=1j), "pml")


mu  = CoefficientFunction(1.0)
eps = CoefficientFunction(1.0 +   \
                          IfPos( (x-0.1)**2 + (y-0.3)**2 - 0.16 , \
                                 1.0 + cos( sqrt( (x-0.1)**2 + (y-0.3)**2 ) * 5.0 * np.pi / 2 ) , 0))
sigma = CoefficientFunction(1.0)

freq = 3 * 2 * np.pi

J_r    = CoefficientFunction(IfPos( (x-0.15)**2 + (y+0.7)**2 - 0.04 , 2, 0)) + \
       CoefficientFunction(IfPos( (x+0.35)**2 + (y-0.4)**2 - 0.04 , 4, 0))

J_c   = CoefficientFunction(IfPos( (x-0.35)**2 + (y-0.5)**2 - 0.04 , 2, 0)) + \
       CoefficientFunction(IfPos( (x+0.65)**2 + (y-0.2)**2 - 0.04 , 4, 0))

class Arg(object):
    pass

args = Arg()
args.frequency    = freq
args.permittivity = eps
args.permeability = mu
args.conductivity = sigma
args.source       = CoefficientFunction((J_r, J_c))
args.mesh         = mesh
args.deg          = 1


m = UMMEG(args)
A, b = m.construct_variational_form()
Sol = m.solve(A, b, "result")