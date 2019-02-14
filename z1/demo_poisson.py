from dolfin import *

# Load mesh
mesh = Mesh("mesh.xml")
subdomains = MeshFunction("size_t", mesh, "mesh_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "mesh_facet_region.xml")
print(mesh)

# define function space
V = FunctionSpace(mesh, "CG", 1)

#define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f1 = Constant(50.0)
f2 = Constant(50.0)
f3 = Constant(50.0)
f4 = Constant(50.0)
k1 = Constant(10.0)
k2 = Constant(100.0)
k3 = Constant(100.0)
k4 = Constant(100.0)
alpha = Constant(1000.0)
g = Constant(10.0)
beta = Constant(100)
h = Constant(40)
dx = Measure('dx', domain=mesh, subdomain_data = subdomains)
ds = Measure('ds', domain=mesh, subdomain_data = boundaries)


#a = inner(k1*grad(u), grad(v))*dx(1) + inner(k2*grad(u), grad(v))*dx(2) + inner(k3*grad(u), grad(v))*dx(3) + inner(k4*grad(u), grad(v))*dx(4)
#L = f1*v*dx(1) + f2*v*dx(2)+ f3*v*dx(3) + f4*v*dx(4)


#a = inner(k1*grad(u), grad(v))*dx(1) + inner(k2*grad(u), grad(v))*dx(2)+alpha*u*v*ds(1)+beta*u*v*ds(3)
#L = f1*v*dx(1) +f2*v*dx(2)+alpha*g*v*ds(1)+beta*h*v*ds(3)

a = inner(k1*grad(u), grad(v))*dx(1) + inner(k2*grad(u), grad(v))*dx(2)+alpha*u*v*ds(1)
L = f1*v*dx(1) +f2*v*dx(2)+alpha*g*v*ds(1)


#Define boundary connidition
#niz
g0 = Constant(0.0)
#pravo
g1 = Constant(-5.0)
#vverh
g2 = Constant(0.0)
#levo
g3 = Constant(100.0)
#bc1 = DirichletBC(V, g0, boundaries, 1)
bc2 = DirichletBC(V, g1, boundaries, 2)
#bc3 = DirichletBC(V, g2, boundaries, 3)
bc4 = DirichletBC(V, g3, boundaries, 4)
#bcs = [bc1, bc3]
bcs = [bc2, bc4]

# Compute solution
u = Function(V)
solve(a == L, u, bcs)
print 'u(%g, %g)' % (u.vector().min(), u.vector().max())

#Save solution in VTK format
file = File("2_poisson.pvd")
file << u

# Plot solution
#import matplotlib.pyplot as plt
#plot(u)
#plt.show()
