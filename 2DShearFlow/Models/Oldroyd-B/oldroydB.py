from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Use UFLACS to speed-up assembly and limit quadrature degree
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4

fileName = "Output"

class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 2.0
        y[1] = x[1]

def writeToFile(w, t):
    g = HDF5File(mesh.mpi_comm(),"out/result.hdf5",'w')
    g.write(w,'u')

    (a,b,c) = w.split()

    a.rename("v", "velocity")
    b.rename("p", "pressure")
    c.rename("Bp", "Btensor")


    # Create files for storing solution
    vfile = XDMFFile("{}/velo{}.xdmf".format(fileName, round(t,2)))
    pfile = XDMFFile("{}/pres{}.xdmf".format(fileName, round(t,2)))
    Bfile = XDMFFile("{}/Btens{}.xdmf".format(fileName, round(t,2)))
   
    vfile.write(a)
    pfile.write(b)
    Bfile.write(c)


mesh = RectangleMesh(Point(0,0), Point(2, 1), 20, 10, "crossed")


inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 2)'
topWall  = 'near(x[1], 1)'
bottomWall = 'near(x[1], 0)'



Ep = FiniteElement("CG",mesh.ufl_cell(),1)
Ev = VectorElement("CG",mesh.ufl_cell(),2)
Eb = TensorElement("CG",mesh.ufl_cell(),1)
#En = VectorElement("CG",mesh.ufl_cell(),2)

# Build function spaces (Taylor-Hood)
W = FunctionSpace(mesh, MixedElement([Ev, Ep, Eb]), constrained_domain=PeriodicBoundary())



# Define boundary conditions
bcu_movement  = DirichletBC(W.sub(0), Constant((2.0, 0)), topWall)
bcu_noslip  = DirichletBC(W.sub(0), Constant((0, 0)), bottomWall)
bcp_in = DirichletBC(W.sub(1), Constant(0.0), outflow)
#bcn_nconst = DirichletBC(W.sub(3), Constant((1,0)), topWall)
bcu = [bcu_movement, bcu_noslip, bcp_in]

v, q, Bpt = TestFunctions(W)
w = Function(W)
u, p, Bp = split(w)

#old timestep
w0 = Function(W)
u0, p0, Bp0 = split(w0)
# Time-stepping parameters
steps = 500
dt = 0.05
t = 0.0

# Facet normal, identity tensor and boundary measure
n = FacetNormal(mesh)
I = Identity(mesh.geometry().dim())
nu = Constant(1)
nu1 = Constant(1)
nu2 = Constant(1)
mu = Constant(1)
c1 = Constant(1)
c2 = Constant(1)
#pozor na gradient v KKK
L = grad(u)

# Define variational forms
T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp - I*tr(Bp))
BEq = nu2*(Bp - I) - L*Bp - Bp*L.T
#NEq = (c1*I + c2*Bp)*np - L*np
F = Constant(1.0/dt)*inner(u-u0, v)*dx + inner(T, grad(v))*dx - q*div(u)*dx + inner(grad(u)*u, v)*dx \
    + nu1*Constant(1.0/dt)*inner(Bp-Bp0, Bpt)*dx + inner(BEq, Bpt)*dx + inner(dot(grad(Bp),u),Bpt)*dx

#assign initial conditions
Bp0 = Expression(("0.0","0.0","0.0","1.0","0.0","0.0","1.0"), degree = 1)
w.assign(interpolate(Bp0, W))
w0.assign(interpolate(Bp0, W))

writeToFile(w,0)

sheerStressFile = open("sheerStress", "w+")

for i in range(1, steps):
    solve(F == 0, w, bcu,  solver_parameters={"newton_solver":{"relative_tolerance":1e-10},"newton_solver":{"maximum_iterations":15}})
    
    t += dt
    print (i)
    if i % 1  < 0.5:
        writeToFile(w, i)

    #Sheer stress
    a_1 = Point(1, 0.5)
    T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp - I*tr(Bp))
    pT = project(T[0,1])
    sheerStressFile.write('{}\t{}\n'.format(round(t,4), pT(a_1)))

    w0.assign(w)

sheerStressFile.close()
#plot(Bp[0,0])

#plt.show()







