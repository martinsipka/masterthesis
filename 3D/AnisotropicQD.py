from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Use UFLACS to speed-up assembly and limit quadrature degree
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4

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

class TopDomain(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1],0.01)

# Define finite elements
mesh = Mesh("/NHCoarse/coneplate.xml")
bndry = MeshFunction("size_t", mesh, "/NHCoarse/coneplate_facet_region.xml")
#plot(mesh, "3DMesh")



Ep = FiniteElement("CG",mesh.ufl_cell(),1)
Ev = VectorElement("CG",mesh.ufl_cell(),2)
Eb = TensorElement("CG",mesh.ufl_cell(),2)
En = VectorElement("CG",mesh.ufl_cell(),2)

# Build function spaces (Taylor-Hood)
W = FunctionSpace(mesh, MixedElement([Ev, Ep, Eb, En]))

#Boundary parts
top_boundary = TopDomain()

#boundary_parts = FacetFunction("size_t", mesh)
top_boundary.mark(bndry,1)


# Define boundary conditions

bcu_noslip  = DirichletBC(W.sub(0), Constant((0, 0, 0)), bndry, 3)

omega = 0.1
v_in = Expression(("angularVelocity * x[1]","-1.0 * angularVelocity * x[0]", "0.0"), angularVelocity = omega , degree = 1)  
bcu_movement = DirichletBC(W.sub(0), v_in, bndry, 1)
bcp_in = DirichletBC(W.sub(1), Constant(0.0), bndry, 2)
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
steps = 250
dt = 0.05
t = 0.0


# Facet normal, identity tensor and boundary measure
ds = Measure("ds")(subdomain_data=bndry)
n = FacetNormal(mesh)
I = Identity(mesh.geometry().dim())
nu = Constant(1)
nu1 = Constant(1)
nu2 = Constant(1)
mu = Constant(1)
const = Constant(1)

#pozor na gradient v KKK




def f(n):
    return (Constant(1)-Constant(1)/inner(n,n))

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
    neq = mu*(Bp - I)*np  + const*f(np)/nu1*Bp*np - 2*nu2*L*np
    return inner(neq, npt)*dx + 2*nu2*inner(dot(grad(np), u), npt)*dx



#theta = Constant(0.5)  # Crank-Nicolson scheme 
theta = Constant(1)  # Implicit Euler scheme 


F = Constant(1.0/dt)*inner(u-u0, v)*dx + theta*a(u, v, p, q, Bp, np) + (1.0-theta)*a(u0, v, p0, q, Bp0, np0) \
    + nu1*Constant(1.0/dt)*inner(Bp-Bp0, Bpt)*dx + theta*b(Bp, Bpt, u, np) + (1.0-theta)*b(Bp0, Bpt, u0, np0) \
    + 2*nu2*Constant(1.0/dt)*inner(np - np0, npt)*dx + theta*c(np, npt, Bp, u) \
    + (1.0 - theta)*c(np0, npt, Bp0, u0)

#assign initial conditions
initC = Expression(("0.0","0.0","0.0",\
    "0.0", \
    "1.0","0.0", "0.0" , "0.0", "1.0", "0.0", "0.0", "0.0","1.0", \
    "0.0", "0.0","1.0"), degree = 2)
w0.assign(interpolate(initC, W))
w.assign(interpolate(initC, W))


with open("sheerStress", "w+") as sheerStressFile: 
    writeToFile(w0, 0)
    sheerStressFile.write('{}\t{}\n'.format(round(t,2), 0.0))

    for i in range (0, steps):
        
        t = t + dt

        solve(F == 0, w, bcu,  solver_parameters={"newton_solver":{"linear_solver":"mumps"},"newton_solver":{"relative_tolerance":1e-15},"newton_solver":{"maximum_iterations":15}})

        print(i)
        if i % 1  < 0.5:
            writeToFile(w, i)

        

        #Sheer stress
        T = -p*I  + 2 * nu * sym(grad(u)) + mu*(Bp -I) + const*f(np)*outer(np, np)
        surfaceNormal = Constant(("0.0", "0.0", "1.0"))
        velNormal = Expression(("x[1]", "-1.0*x[0]", "0.0"), degree=2)
        velNormal = velNormal/sqrt(inner(velNormal, velNormal))
        #print(velNormal)
        sheerStress = inner(dot(T, surfaceNormal), velNormal)
        sheerForce = assemble(sheerStress*ds(1,subdomain_data=bndry))
        sheerRate = omega/0.04
        visc = sheerForce/sheerRate
        
        sheerStressFile.write('{}\t{}\n'.format(round(t,2), visc))
        sheerStressFile.flush()

        w0.assign(w)


