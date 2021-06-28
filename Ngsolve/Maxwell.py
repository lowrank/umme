from UMMEG import *

'''
In the following example, the domain is set in a unit disk and extended with PML of width 0.4. 

The maximum h of each element is set as 0.1, see the code.
'''
geometry = SplineGeometry()

geometry.AddCircle ((0, 0), r=2.0, leftdomain=1)
# geometry.AddCircle ((0, 0), r=5, leftdomain=2, bc="pmlBC")

geometry.SetMaterial(1, "domain")

ngmesh = geometry.GenerateMesh(maxh=0.03)
mesh = Mesh(ngmesh)
mesh.Curve(3)

# mesh.SetPML(pml.Radial(origin=(0,0), rad=1.0, alpha=0.1j), definedon=2)


mu  = CoefficientFunction(1.0)
eps = CoefficientFunction(1.0 +   \
                          0.4 *IfPos( (x-0.1)**2 + (y-0.3)**2 - 0.16 , \
                                 1.0 + cos( sqrt( (x-0.1)**2 + (y-0.3)**2 ) * 5.0 * np.pi / 2 ) , 0))
sigma = CoefficientFunction(
    0.8 * IfPos((x-0.05)**2 + (y-0.2)**2 - 0.16,
                1.0 + cos(sqrt((x-0.05)**2 + (y-0.2)**2) * 5.0 * np.pi / 2), 0))

freq = 4 * 2 * np.pi

J_r    = CoefficientFunction(
                    0.5 * IfPos((x-0.15)**2 + (y+0.2)**2 - 0.25,
                                1.0 + cos(sqrt((x-0.15)**2 + (y+0.2)**2 ) * 4.0 * np.pi / 2), 0))

J_c = CoefficientFunction(
    0.2 * IfPos((x-0.35)**2 + (y-0.5)**2 - 0.04,
                1.0 + cos(sqrt((x-0.35)**2 + (y-0.5)**2) * 10.0 * np.pi / 2), 0))


class Arg(object):
    pass

args = Arg()
args.frequency    = freq
args.permittivity = eps
args.permeability = mu
args.conductivity = sigma
args.source       = CoefficientFunction((J_r, J_c))
args.mesh         = mesh
args.deg          = 3


m = UMMEG(args)
A, b = m.construct_variational_form()
Sol = m.solve(A, b, "result")
