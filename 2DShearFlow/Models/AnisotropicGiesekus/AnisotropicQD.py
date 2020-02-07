from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Use UFLACS to speed-up assembly and limit quadrature degree
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4

class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 2.0
        y[1] = x[1]

fileName = "Output"

def writeToFile(w, t):
    (a,b,c,d) = w.split()

    a.rename("v", "velocity")
    b.rename("p", "pressure")
    c.rename("Bp", "Btensor")
    d.rename("n", "n vector")
  

    # Create files for storing solution
    vfile = XDMFFile("{}/velo{}.xdmf".format(fileName, round(t,2)))
    pfile = XDMFFile("{}/pres{}.xdmf".format(fileName, round(t,2)))
    Bfile = XDMFFile("{}/Btens{}.xdmf".format(fileName, round(t,2)))
    nfile = XDMFFile("{}/nvect{}.xdmf".format(fileName, round(t,2)))

    vfile.write(a)
    pfile.write(b)
    Bfile.write(c)
    nfile.write(d)

mesh = RectangleMesh(Point(0,0), Point(2, 1), 20, 10, "crossed")


inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 2)'
topWall  = 'near(x[1], 1)'
bottomWall = 'near(x[1], 0)'



Ep = FiniteElement("CG",mesh.ufl_cell(),1)
Ev = VectorElement("CG",mesh.ufl_cell(),2)
Eb = TensorElement("CG",mesh.ufl_cell(),2)
En = VectorElement("CG",mesh.ufl_cell(),2)

# Build function spaces (Taylor-Hood)
W = FunctionSpace(mesh, MixedElement([Ev, Ep, Eb, En]), constrained_domain=PeriodicBoundary())

# Define boundary conditions
bcu_movement  = DirichletBC(W.sub(0), Constant((2.0, 0)), topWall)
bcu_noslip  = DirichletBC(W.sub(0), Constant((0, 0)), bottomWall)
bcp_in = DirichletBC(W.sub(1), Constant(0.0), outflow)
#bcn_nconst = DirichletBC(W.sub(4), Constant((0.0,1.0)), topWall)
#bcn_nconst2 = DirichletBC(W.sub(4), Constant((1.0,0.0)), bottomWall)
bcu = [bcu_movement, bcu_noslip, bcp_in]#, bcn_nconst, bcn_nconst2]

v, q, Bpt, npt = TestFunctions(W)
w = Function(W)
u, p, Bp, np = split(w)

#old timestep
w0 = Function(W)
u0, p0, Bp0, np0 = split(w0)
# Time-stepping parameters
steps = 500
dt = 0.05
t = 0.0




# Facet normal, identity tensor and boundary measure
n = FacetNormal(mesh)
I = Identity(mesh.geometry().dim())
nu = Constant(1)
nu1 = Constant(1)
nu2 = Constant(0.1)
mu = Constant(1)
const = Constant(1)
#pozor na gradient v KKK

def f(n):
    return (Constant(1)-Constant(1)/inner(n,n))

#def f(n):
#    return (inner(n,n) - Constant(1))

# Define steady part of the equation
def a(u, v, p, q, Bp, np):
    T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp -I) + const*f(np)*outer(np, np)
    return inner(T, grad(v))*dx  + inner(grad(u)*u, v)*dx - q*div(u)*dx

def b(Bp, Bpt, u, np):
    L = grad(u)
    BEq = mu*(Bp - I)*Bp - nu1*(L*Bp + Bp*L.T) + const*f(np)*outer(np,np)
    return inner(BEq, Bpt)*dx + nu1*inner(dot(grad(Bp), u), Bpt)*dx

def c(np, npt, Bp, u):
    L = grad(u)
    neq = mu*nu2/nu1*(Bp - I)*np  + const*f(np)/2/nu1*((nu2+nu1)*Bp + (nu2-nu1)*inner(np,np)*I)*np - 2*nu2*L*np
    return inner(neq, npt)*dx + 2*nu2*inner(dot(grad(np), u), npt)*dx



theta = Constant(1.0) # Euler scheme 
#theta = Constant(0.5)  # Crank-Nicolson scheme


F = Constant(1.0/dt)*inner(u-u0, v)*dx + theta*a(u, v, p, q, Bp, np) + (1.0-theta)*a(u0, v, p0, q, Bp0, np0) \
    + nu1*Constant(1.0/dt)*inner(Bp-Bp0, Bpt)*dx + theta*b(Bp, Bpt, u, np) + (1.0-theta)*b(Bp0, Bpt, u0, np0) \
    + 2*nu2*Constant(1.0/dt)*inner(np - np0, npt)*dx + theta*c(np, npt, Bp, u) \
    + (1.0 - theta)*c(np0, npt, Bp0, u0)

#assign initial conditions
initC = Expression(("0.0","0.0",\
    "0.0", \
    "1.0","0.0","0.0","1.0", \
    "0.0","1.0"), degree = 2)
w0.assign(interpolate(initC, W))
w.assign(interpolate(initC, W))


with open("sheerStress", "w+") as sheerStressFile:
    sheerStressFile.write('{}\t{}\n'.format(round(t,2), 0))
    writeToFile(w0, 0)

    for i in range (1, steps+1):
        
        t = t + dt

        solve(F == 0, w, bcu,  solver_parameters={"newton_solver":{"linear_solver":"mumps"},"newton_solver":{"relative_tolerance":1e-15},"newton_solver":{"maximum_iterations":15}})

        print (i)
        if i % 10  < 0.5:
            writeToFile(w, i)

        #Sheer stress
        a_1 = Point(1, 0.5)
        T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp -I) + const*f(np)*outer(np, np)
        pT = project(T[0,1])
        sheerStressFile.write('{}\t{}\n'.format(round(t,2), pT(a_1)))

        w0.assign(w)
   

