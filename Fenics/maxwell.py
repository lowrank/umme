from MEG import *

# domain 
mesh = RectangleMesh(Point(0., 0.), Point(1., 1.), 200, 200)

freq    = 4 * (2 * np.pi)
Lbd     = 1.0
dim     = 2

epsilon = Constant(1)
mu      = Constant(1.0)
sigma   = Constant(0.0000)


class MyFunctionExpression(UserExpression):
    def eval(self, values, x):
        if (x[0]-0.15)**2 + (x[1]-0.3)**2 - 0.04 > 0:
            values[0] = 0.5 * ( 1.0 + cos(sqrt((x[0]-0.15)**2 + (x[1]-0.3)**2) * 10.0 * np.pi / 2))
        else:
            values[0] = 0

        if (x[0]-0.35)**2 + (x[1]-0.5)**2 - 0.04 > 0:
            values[1] = 0.2 * (1.0 + cos(sqrt((x[0]-0.35)**2 + (x[1]-0.5)**2) * 10.0 * np.pi / 2))
        else:
            values[1] = 0

    def value_shape(self):
        return (2,)

fr      = Constant((0.0, 0.0))
fc      = Constant((0.0, 0.0))
gr      = Constant((0.0, 0.0)) # Expression(('x[1]*x[1]', 'x[0]*x[0]'), degree=2)
gc      = Constant((0.0, 0.0))


def dbc(x):
    return x[1] < DOLFIN_EPS 

def ibc(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS

class Arg(object):
    pass 

args = Arg()
args.frequency    = freq
args.permittivity = epsilon 
args.permeability = mu
args.conductivity = sigma 
args.source_r     = fr 
args.source_c     = fc 
args.imp_r        = gr
args.imp_c        = gc 
args.Lambda       = Lbd
args.mesh         = mesh 
args.dim          = dim 

# args.dbc          = dbc
# args.ibc          = ibc

args.dbc = None
args.ibc = lambda x: True

args.scalar_elem  = 'CG'
args.vector_elem  = 'N1curl'
args.deg          = 3

m = MEG(args)

m.source_c = freq * MyFunctionExpression(element=m.vector_space.ufl_element())

A, b = m.construct_variational_form()

Er, Ec = m.solve(A, b)


plt.figure(1)
m.visualize(Er, Ec)
plt.show()
interactive(True)
