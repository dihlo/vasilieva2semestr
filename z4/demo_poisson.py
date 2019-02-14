from dolfin import *

mesh = Mesh( 'mesh.xml' )
subdomains = MeshFunction ('size_t', mesh, 'mesh_physical_region.xml' )
boundaries = MeshFunction ('size_t', mesh, 'mesh_facet_region.xml' )

class HeterExpression(UserExpression):
	def setvalues(self, u1, u2):
		self.u1 = u1
		self.u2 = u2
	def eval_cell (self, value, x, ufc_cell):
		if subdomains [ufc_cell.index] == 1:
			value [0] = self.u1
		else:	
			value [0] = self.u2
k = HeterExpression(degree = 0)
k.setvalues(200.0, 10.0)
mu = HeterExpression(degree = 0)
mu.setvalues(0.8e9 , 0.5e9)
lmbda = HeterExpression(degree = 0)
lmbda.setvalues(1.25e9 , 0.65e9)
beta = 1.0e-5 * (3 * lmbda + 2* mu)
c = 1.0e6
alpAir = 100.0

TAir = 10.0
THot = 150.0

tmax = 100000
dt = tmax/50

def epsilon (u):
	return 0.5 * (grad(u) + grad(u).T)
def sigma (u):
	return lmbda * div(u) * Identity (2) + 2 * mu * epsilon (u)

V = VectorElement("CG", mesh.ufl_cell(), 1)
Q = FiniteElement("CG", mesh.ufl_cell(), 1)
TH = V * Q

W = FunctionSpace(mesh, TH)

un_val = Constant((0.0, 0.0))
Tn_val = Constant(TAir)

un = interpolate(un_val , W.sub(0).collapse ())
Tn = interpolate(Tn_val , W.sub(1).collapse ())

w = Function(W)
g0 = Constant(0.0)
bc2 = DirichletBC(W.sub(0).sub (0), g0, boundaries, 2)
bc3 = DirichletBC(W.sub(0).sub (1), g0, boundaries, 1)
bcT1 = DirichletBC(W.sub(1), THot, boundaries, 3)
bcT3 = DirichletBC (W.sub(1), TAir, boundaries, 1)
bcs = [bc2 , bc3 , bcT1 , bcT3]

u, T = TrialFunctions(W)
v, q = TestFunctions(W)
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

a = inner(sigma (u), epsilon (v)) * dx \
+ beta * inner(grad(T), v) * dx \
+ beta * Tn_val * div(u) * q * dx + c * T * q * dx \
+ dt * inner(k * grad(T), grad(q)) * dx \
+ dt * alpAir * T * q * ds (4)
L = c * Tn * q * dx + beta * Tn_val * div(un) * q * dx + dt * alpAir * TAir * q * ds (4)
w.rename('u', '0')
fileu = File("./results/u.pvd")
fileT = File("./results/T.pvd")
t = 0
while t < tmax:
	print(t)
	t += dt
	solve(a == L, w, bcs)
	(uu , TT) = w. split ( deepcopy = True)
	un.assign (uu)
	Tn.assign (TT)
	print( ' w(%g, %g) ' % (w. vector (). min (), w. vector (). max ()))
	fileu << uu
	fileT << TT