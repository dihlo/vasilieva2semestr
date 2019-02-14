from dolfin import *

mesh = Mesh( 'mesh.xml' )
subdomains = MeshFunction ('size_t', mesh, 'mesh_physical_region.xml' )
boundaries = MeshFunction ('size_t', mesh, 'mesh_facet_region.xml' )

mu = 0.8 * 1.0e11
lmbda = 1.25 * 1.0e11
f = Constant((0.0 , -1.0e8))
g = Constant(0.0)

def epsilon (u):
	return 0.5 * (grad(u) + grad(u).T)
def sigma (u):
	return lmbda * div(u) * Identity(2) + 2*mu*epsilon(u)

V = VectorFunctionSpace(mesh , "CG", 1)
bc2 = DirichletBC(V.sub (0) , g, boundaries , 2)
bc3 = DirichletBC(V.sub (1) , g, boundaries , 3)
bcs = [bc2 , bc3]

u = TrialFunction(V)
v = TestFunction(V)
dx = Measure('dx', domain=mesh, subdomain_data = subdomains)
ds = Measure('ds', domain=mesh, subdomain_data = boundaries)
k_1 = Constant(20)
k_2 = Constant(1)
a = k_1*inner(sigma(u), epsilon(v))*dx(1) + k_2*inner(sigma(u), epsilon(v))*dx(2)
L = inner(f,v)*ds(1)


u = Function(V)
solve (a == L, u, bcs)
print('u(%g, %g)' % (u.vector().min(), u.vector().max ()))
u. rename ('u', '0')
file = File("./results/u.pvd")
file << u
W = TensorFunctionSpace (mesh , "DG", 0)
stress = project(sigma(u), V=W)
stress.rename('u', '0')
File("./results/s.pvd") << stress