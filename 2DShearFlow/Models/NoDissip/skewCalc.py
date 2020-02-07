from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Use UFLACS to speed-up assembly and limit quadrature degree
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4

fileName = "output"

class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 2.0
        y[1] = x[1]

def writeToFile(w, np, t):
    (a,b,c,d, e) = w.split()

    a.rename("v", "velocity")
    b.rename("p", "pressure")
    c.rename("Bp", "Btensor")
    d.rename("a", "Alpha")
    e.rename("l", "Length")
    np.rename("n", "nvector")

    # Create files for storing solution
    vfile = XDMFFile("{}/velo{}.xdmf".format(fileName, round(t,3)))
    pfile = XDMFFile("{}/pres{}.xdmf".format(fileName, round(t,3)))
    Bfile = XDMFFile("{}/Btens{}.xdmf".format(fileName, round(t,3)))
    afile = XDMFFile("{}/alpha{}.xdmf".format(fileName, round(t,3)))
    lfile = XDMFFile("{}/length{}.xdmf".format(fileName, round(t,3)))
    nfile = XDMFFile("{}/nvect{}.xdmf".format(fileName, round(t,3)))

    vfile.write(a)
    pfile.write(b)
    Bfile.write(c)
    afile.write(d)
    lfile.write(e)
    nfile.write(np)
   


mesh = RectangleMesh(Point(0,0), Point(2, 1), 20, 10, "crossed")


inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 2)'
topWall  = 'near(x[1], 1)'
bottomWall = 'near(x[1], 0)'



Ep = FiniteElement("CG",mesh.ufl_cell(),1)
Ev = VectorElement("CG",mesh.ufl_cell(),2)
Eb = TensorElement("CG",mesh.ufl_cell(),2)
Ea = FiniteElement("CG",mesh.ufl_cell(),2)
El = FiniteElement("CG",mesh.ufl_cell(),2)

# Build function spaces (Taylor-Hood)
W = FunctionSpace(mesh, MixedElement([Ev, Ep, Eb, Ea, El]), constrained_domain=PeriodicBoundary())

# Define boundary conditions
bcu_movement  = DirichletBC(W.sub(0), Constant((2.0, 0)), topWall)
bcu_noslip  = DirichletBC(W.sub(0), Constant((0, 0)), bottomWall)
bcp_in = DirichletBC(W.sub(1), Constant(0.0), outflow)
#bcn_nconst = DirichletBC(W.sub(4), Constant((0.0,1.0)), topWall)
#bcn_nconst2 = DirichletBC(W.sub(4), Constant((1.0,0.0)), bottomWall)
bcu = [bcu_movement, bcu_noslip, bcp_in]#, bcn_nconst, bcn_nconst2]

v, q, Bpt, at, lt = TestFunctions(W)
w = Function(W)
u, p, Bp, al, lam = split(w)

#old timestep
w0 = Function(W)
u0, p0, Bp0, al0, lam0 = split(w0)
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
const = Constant(1)
#pozor na gradient v KKK

def n(alpha):
    return as_vector((cos(alpha), sin(alpha)))

def f(n):
    return (Constant(1)-Constant(1)/inner(n,n))

def cross2D(aV, bV):
    return aV[0]*bV[1] - aV[1]*bV[0]

#def f(n):
#    return (inner(n,n) - Constant(1))

# Define steady part of the equation
def a(u, v, p, q, Bp, np):
    T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp -I) + const*f(np)*outer(np, np)
    return inner(T, grad(v))*dx  + inner(grad(u)*u, v)*dx - q*div(u)*dx

def b(Bp, Bpt, u, np):
    L = grad(u)
    BEq = mu*(Bp - I) - nu1*(L*Bp + Bp*L.T) + const*f(np)*outer(np,np)
    return inner(BEq, Bpt)*dx + inner(dot(grad(Bp), u), Bpt)*dx

def c(Bp, ns, at):
    return inner(cross2D(dot(Bp,ns), ns), at)*dx

def df(Bp, ns, lam, lt):
    return inner( inner(dot(Bp,ns), ns) - lam*lam, lt)*dx

#theta = Constant(0.5) # Crank-Nicolson scheme
theta = Constant(1.0) # Implicit Euler scheme

ns = n(al)
ns0 = n(al0)
np = lam*ns
np0 = lam*ns0
print(np)

F = Constant(1.0/dt)*inner(u-u0, v)*dx + theta*a(u, v, p, q, Bp, np) + (1.0-theta)*a(u0, v, p0, q, Bp0, np0) \
    + nu1*Constant(1.0/dt)*inner(Bp-Bp0, Bpt)*dx + theta*b(Bp, Bpt, u, np) + (1.0-theta)*b(Bp0, Bpt, u0, np0) \
    + c(Bp, ns, at) + df(Bp, ns, lam, lt)

#assign initial conditions
initC = Expression(("0.5","0.0",\
    "0.0", \
    "0.842287","0.23427","0.23427","0.681137", \
    #"1.0","0.0","0.0","1.0", \
    #str(1/sqrt(2.0)), str(1/sqrt(2.0)), \
    #str(1/sqrt(2.0)), str(1/sqrt(2.0)), \
    "0.619765", "1.0047"), degree = 2) \
    #"0.0", "1.0", \
    #"0.0","0.5","-0.5","0.0"), degree = 2)
w0.assign(interpolate(initC, W))
w.assign(interpolate(initC, W))

#problem = NonlinearVariationalProblem(lhs(F), rhs(F), u, bcu)
#solver  = NonlinearVariationalSolver(problem)
#solver.solve()

sheerStressFile = open("sheerStress", "w+")
angleFile = open("angleFile", "w+")
periodFile = open("periodFile", "w+")

sheerStressFile.write("0\t0\n")

writeToFile(w0, project(lam*n(al)), 0)

for i in range (1, steps):
    
    t = t + dt

    solve(F == 0, w, bcu,  solver_parameters={"newton_solver":{"linear_solver":"mumps"},"newton_solver":{"relative_tolerance":1e-15},"newton_solver":{"maximum_iterations":15}})

    print (i)
    if i % 1  < 0.5:
        writeToFile(w, project(lam*n(al)), i)
    
    

    #Sheer stress
    a_1 = Point(1, 0.5)
    np = lam*n(al)
    T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp -I) + const*f(np)*outer(np, np)
    pT = project(T[0,1])
    sheerStressFile.write('{}\t{}\n'.format(round(t,4), pT(a_1)))

    #angle
    nAlpha = 0
    nNumeric = project(np)
    nNumeric = nNumeric(a_1)
    if nNumeric[0] < 10*DOLFIN_EPS:
        if nNumeric[1] > 0: nAlpha = pi/4
        if nNumeric[1] < 0: nAlpha = 3*pi/4
    nAlpha = atan(nNumeric[1]/nNumeric[0])
    angleFile.write('{}\t{}\n'.format(round(t,4), round(nAlpha,7)))
    

    #Period of rotation
    nAlpha0 = 0
    nNumeric0 = project(np0)
    nNumeric0 = nNumeric0(a_1)
    if nNumeric0[0] < 10*DOLFIN_EPS:
        if nNumeric0[1] > 0: nAlpha0 = pi/4
        if nNumeric0[1] < 0: nAlpha0 = 3*pi/4
    nAlpha0 = atan(nNumeric0[1]/nNumeric0[0])
    periodFile.write('{}\t{}\n'.format(round(t,4), round(abs(nAlpha - nAlpha0)/dt,7)))

    w0.assign(w)



