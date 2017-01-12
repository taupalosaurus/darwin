from firedrake import *

class LinearShallowWaterSolver(object):
    def __init__(self, mesh, depth, dt,
            theta=1.0, # constant in theta time-integration scheme: 0=explicit, 1=implicit
            g=9.81, # gravity constant
            ):
        """Set up a linear shallow water problem on the given mesh."""

        mesh = mesh
        # vector function space for velocity: P1_DG, i.e. piecewise linear and discontinuous
        V = VectorFunctionSpace(mesh, "DG", 1)
        # function space for free surface elevation: P2, i.e. piecewise quadratic continuous
        H = FunctionSpace(mesh, "CG", 2)
        # combined function space for (velocity, pressure) solution
        Z = V * H

        u, eta = TrialFunctions(Z)
        v, q = TestFunctions(Z)
        # z0 = (u0, eta0) contains the solution from the previous time-step
        self.z0 = Function(Z)
        self.u0, self.eta0 = self.z0.split()
        self.u0.rename("Velocity")
        self.eta0.rename("Freesurface")

        # we need to use a Constant() instead of directly using the float,
        # because firedrake gets confused if theta==0, or 1-theta==0
        theta_const = Constant(theta)

        # normal to the boundary:
        n = FacetNormal(mesh)

        # theta-weighted variables:
        eta_theta = theta * eta + (1.0-theta) * self.eta0
        u_theta = theta * u + (1.0-theta) * self.u0

        # the complete weak form of the equations:
        F = (
             # time derivative of depth-averaged momentum eqn:
             dot(v, u - self.u0) * dx
             # pressure gradient term
             + dt * theta_const * dot(v, g*grad(eta_theta)) * dx
             # eta=0 boundary condition on boundary marked with id 1
             # (which is the left boundary when using RectangleMesh)
             + dt * dot(v, n) * g * (eta_theta - 0.0) * ds(1)
             # time derivative of continuity equation:
             + q * (eta - self.eta0) * dx
             # q * div(depth * u) integrated by parts:
             - dt * theta_const * dot(grad(q), depth * u_theta) * dx
             # integration by parts gives the following boundary term on the
             # open boudaries:
             + dt * q * dot(n, u_theta) * ds(1)
             # on the other boundaries (closed) we impose u.n=0 by substitution
             # thus these boundary terms drop out
            )

        # split in lhs and rhs
        a, L = lhs(F), rhs(F)

        # the problem we want to solve, note that the solution is directly
        # written into z0 so that it can be used as the "previous timestep" solution
        # in the next solve
        prob = LinearVariationalProblem(a, L, self.z0)
        self.solver = LinearVariationalSolver(prob)

        M = TensorFunctionSpace(mesh, "CG", 1);
        sigma = TestFunction(M)
        self.H = Function(M)
        Lh = inner(sigma, self.H)*dx + inner(div(sigma), grad(self.eta0))*dx
        Lh -= (sigma[0, 1]*n[1]*self.eta0.dx(0) + sigma[1, 0]*n[0]*self.eta0.dx(1))*ds
        H_prob = NonlinearVariationalProblem(Lh, self.H)
        self.H_solv = NonlinearVariationalSolver(H_prob, solver_parameters={'snes_rtol': 1e-2,
                                                                   'ksp_rtol': 1e-5,
                                                                   'ksp_gmres_restart': 20,
                                                                   'pc_type': 'sor',
                                                                   'snes_monitor': True,
                                                                   'snes_view': False,
                                                                   'ksp_monitor_true_residual': False,
                                                                   'snes_converged_reason': True,
                                                                   'ksp_converged_reason': True})


    def get_solution(self):
        """Get the velocity and elevation that contain the current solution.
        These can also be set to different values, for instance to provide 
        an initial solution."""
        return self.u0, self.eta0

    def iterate(self):
        """Perform one timestep."""
        self.solver.solve()

    def get_hessian(self):
        return self.H

    def compute_hessian(self):

        self.H_solv.solve()

