from dolfin import *
mesh = Mesh( 'mesh.xml' )
domains = MeshFunction ('size_t', mesh, 'mesh_physical_region.xml' )
boundaries = MeshFunction ('size_t', mesh, 'mesh_facet_region.xml' )
V = FunctionSpace (mesh ,'CG', 1)
dx = Measure ('dx')(subdomain_data = domains)
ds = Measure ('ds')(subdomain_data = boundaries)
f = Constant (0.0)
g = Constant (80.0)
u0_val = Constant (5.0)
u0 = interpolate (u0_val , V)

#Define boundary connidition
#niz
g0 = Constant(0.0)
#pravo
g1 = Constant(0.0)
#vverh
g2 = Constant(0.0)
#levo
g3 = Constant(100.0)
#bc1 = DirichletBC(V, g0, boundaries, 1)
bc2 = DirichletBC(V, g1, boundaries, 2)
#bc3 = DirichletBC(V, g2, boundaries, 3)
bc4 = DirichletBC(V, g3, boundaries, 4)
#bcs = [bc1, bc3]
bcs = [bc4]

k_1 = Constant (20)
k_2 = Constant (1)
C_1 = Constant (150)
C_2 = Constant (100)
T = 60.0
N = 100
tau = T/N
u = TrialFunction (V)
v = TestFunction (V)
a = (C_1 / tau) * u * v * dx (1) + (C_2 / tau) * u * v * dx (2) + k_1 * inner (grad(u), grad(v)) * dx (1) + k_2 * inner (grad(u), grad(v)) * dx (2)
L = (C_1 / tau) * u0 * v * dx (1) + (C_2 / tau) * u0 * v * dx (2) + f * v * dx (1) + f * v * dx (2)
u = Function (V)
file = File('./results/time_dep.pvd')
t = 0
while t < T:
	t += tau
	solve (a == L, u, bcs)
	file << u
	u0. assign (u)