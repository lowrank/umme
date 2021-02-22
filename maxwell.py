from MEG import *

# domain 
mesh = RectangleMesh(Point(0, 0), Point(1, 1), 160, 160)

freq    = 30 * (2 * pi)
Lbd     = 1.0
dim     = 2

epsilon = Constant(1e-9)
mu      = Constant(1.0)
sigma   = Constant(0.0001)

fr      = Constant((0.0, 0.0))
fc      = Constant((0.0, 0.0))
gr      = Expression(( 'x[1]*x[1]', 'x[0]*x[0]' ), degree = 2)
gc      = Constant((1.0, 1.0))

def dbc(x):
    return x[1] < DOLFIN_EPS 

def ibc(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

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
args.deg          = 2 

m = MEG(args)

A, b = m.construct_variational_form()

Er, Ec = m.solve(A, b)


plt.figure(1)
m.visualize(Er, Ec)
plt.show()
interactive(True)
