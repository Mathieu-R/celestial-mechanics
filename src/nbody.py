import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

from mpl_toolkits.mplot3d import Axes3D

from .edo import (hamiltonian, compute_angular_momentum, n_body_dqdt, n_body_dpdt)
from .solvers import (heun, euler_symp, stormer_verlet)
from consts import (au_to_meter, day_to_second)

from utils import (set_size, set_size_square_plot)

class NBodySimulation():
  def __init__(self, bodies, t0, tN, dt, options):
    self.options = options
    self.bodies = bodies
    self.t0 = t0
    self.tN = tN
    self.dt = dt

    self.total_mass = sum(body.mass for body in self.bodies)
    self.mean_pos = sum(body.initial_positions * body.mass for body in self.bodies) / self.total_mass
    self.mean_vel = sum(body.initial_impulsions for body in self.bodies) / self.total_mass

    # shift the coordinate frame so that the barycenter is at rest.
    for body in self.bodies:
      body.initial_positions -= self.mean_pos
      body.initial_impulsions -= body.mass * self.mean_vel

    # number of time step
    self.nt = int((self.tN - self.t0) / self.dt)
    self.time_mesh = np.linspace(start=self.t0, stop=self.tN, num=self.nt)
    self.legends = ["Heun (RK2)", "Euler Symplectique", "Stormer-Verlet"]

    self.solvers = [
      {"call": heun, "name": "Heun (RK2)", "color": "teal", "bodies": self.bodies},
      {"call": euler_symp, "name": "Euler Symplectique", "color": "darkorange", "bodies": self.bodies},
      {"call": stormer_verlet, "name": "Stormer Verlet", "color": "crimson", "bodies": self.bodies}
    ]

    self.solvers2 = {
      "rk2": {
        "call": heun,
        "name": "Heun (RK2)",
        "color": "teal",
        "bodies": self.bodies
      },
      "euler-sympectic": {
        "call": euler_symp,
        "name": "Euler Symplectique",
        "color": "darkorange",
        "bodies": self.bodies
      },
      "stormer-verlet": {
        "call": stormer_verlet,
        "name": "Stormer Verlet",
        "color": "crimson",
        "bodies": self.bodies
      }
    }

    if self.options["save"]:
      self.figure_options = set_size(width="full-size", subplots=(1,3))
    else:
      self.figure_options = 8,8

    self.results = []

  def solve(self, solver, dt, nt, bodies):
    # positions and impulsions state vectors
    # each body is in \R^3 (3D space x, y, z)
    q = np.zeros((self.nt, len(self.bodies) * 3))
    p = np.zeros((self.nt, len(self.bodies) * 3))

    # energy mesh -- computed from the hamiltonian
    energy = np.zeros(self.nt)
    # angular momentum mesh
    angular_momentum = np.zeros((self.nt, 3))

    # set initial conditions
    q[0] = np.concatenate(np.array([body.initial_positions for body in bodies]))
    p[0] = np.concatenate(np.array([body.initial_impulsions for body in bodies]))

    # set initial energy
    #print(energy)
    #print(hamiltonian(qk=q[0], pk=p[0], bodies=bodies))
    energy[0] = hamiltonian(qk=q[0], pk=p[0], bodies=bodies)
    angular_momentum[0] = compute_angular_momentum(qk=q[0], pk=p[0], bodies=bodies)

    return solver(dqdt=n_body_dqdt, dpdt=n_body_dpdt, q=q, p=p, dt=dt, nt=nt, bodies=bodies, energy=energy, angular_momentum=angular_momentum)

  def simulate(self):
    self.results = []
    for solver in self.solvers:
      q, p, energy, angular_momentum = self.solve(solver=solver["call"], dt=self.dt, nt=self.nt, bodies=self.bodies)
      self.results.append({"solver": solver["name"], "color": solver["color"], "q": q, "p": p, "energy": energy, "angular_momentum": angular_momentum})

  def plot2D(self):
    self.fig = plt.figure(figsize=(set_size_square_plot(width="full-size", subplots=(2,2))))
    # loop for each result (corresponding to a specific solving method)
    for (index, result) in enumerate(self.results):
      # create a 2D plot
      ax = self.fig.add_subplot(2, 2, index + 1)
      # set plot parameters
      ax.set_title(result["solver"], fontsize=8)
      #ax.text(0.5, 1.05, f"({index + 1})", ha="center", transform=ax.transAxes, size=8)

      # plotting each body
      for (ind, body) in enumerate(self.bodies):
        x_index = (ind * 3)
        y_index = (ind * 3) + 1

        x, y = result["q"][:,x_index], result["q"][:,y_index]
        ax.plot(x, y, c=body.color, label=body.name, marker=body.marker)

      # set labels
      ax.set_xlabel("x (AU)")
      ax.set_ylabel("y (AU)")

    plt.tight_layout()

    if self.options["save"]:
      self.fig.savefig("figures/orbital-plot2d.png")

    plt.show()

  def plot3D(self):
    self.fig = plt.figure(figsize=(set_size_square_plot(width="column-size", subplots=(1,1))))
    # create subplot grid
    # loop for each result (corresponding to a specific solving method)
    for (index, result) in enumerate(self.results):
      # create a 3D plot
      ax = self.fig.add_subplot(1, 1, 1, projection="3d")
      # set plot parameters
      #ax.set_title(result["solver"])

      max_range = 0

      # plotting each body
      for (ind, body) in enumerate(self.bodies):
        x_index = (ind * 3)
        y_index = (ind * 3) + 1
        z_index = (ind * 3) + 2

        x, y, z = result["q"][:,x_index], result["q"][:,y_index], result["q"][:,z_index]

        max_dim = max(max(x), max(y), max(z))
        if max_dim > max_range:
          max_range = max_dim

        ax.plot(xs=x, ys=y, zs=z, c=body.color, label=body.name)

      # limiting plot
      ax.set_xlim([-max_range,max_range])
      ax.set_ylim([-max_range,max_range])
      ax.set_zlim([-max_range,max_range])

      # change tick font
      for tick in ax.get_xticklabels(which="both"):
        tick.set_fontname("Roboto Condensed")

      # set labels
      ax.set_xlabel("x (AU)")
      ax.set_ylabel("y (AU)")
      ax.set_zlabel("z (AU)")

      ax.set_xticklabels([])
      ax.set_yticklabels([])
      ax.set_aspect('auto')

      if self.options["save"]:
        self.fig.savefig(f"figures/orbital-plot3d--{result['solver']}.pdf")

    plt.show()

  def plot_energy(self):
    self.fig = plt.figure(figsize=(set_size(width="full-size")))
    ax = self.fig.add_subplot(1, 1, 1)

    ax.set_title(f"Energie", fontsize=8)

    for (index, result) in enumerate(self.results):
      solver_name = result["solver"]
      # kg * au^2 / day^2 => kg * m^2 / s^2
      energy = result["energy"] * ((au_to_meter ** 2) / (day_to_second ** 2))
      ax.plot(self.time_mesh / 365.25, energy, c=result["color"])
      ax.set_xlabel("Temps [années]")
      ax.set_ylabel("Energie totale [$J$]")


    if self.options["save"]:
      self.fig.savefig('figures/orbital-energy.pdf')

    plt.show()

  def plot_angular_momentum(self):
    self.fig = plt.figure(figsize=(set_size(width="full-size", subplots=(2,2))))

    for (index, result) in enumerate(self.results):
      ax = self.fig.add_subplot(2, 2, index + 1)

      solver_name = result["solver"]
      # kg * au^2 / day => kg * m^2 / s
      angular_momentum = result["angular_momentum"] * ((au_to_meter ** 2) / (day_to_second))

      for (ind, body) in enumerate(self.bodies):
        ax.plot(self.time_mesh / 365.25, angular_momentum[:,ind], c=body.color, label=body.name)

      ax.tick_params(axis='x', labelsize=10)
      ax.tick_params(axis='y', labelsize=10)

      ax.set_xlabel("Temps [années]", fontsize=10)
      ax.set_ylabel("L [$kg.m^2.s^{-1}$]", fontsize=10)
      #ax.text(0.5, 1.05, f"({index + 1})", ha="center", transform=ax.transAxes, size=8)
      ax.set_title(f"Moment angulaire: {solver_name}", fontsize=8)

    plt.tight_layout()

    if self.options["save"]:
      self.fig.savefig('figures/orbital-angular-momentum.pdf')

    plt.show()

  def update(self, index):
    for i in range(len(self.bodies)):
      x_index = (i * 3)
      y_index = (i * 3) + 1
      z_index = (i * 3) + 2

      x = self.q[:index,x_index]
      y = self.q[:index,y_index]
      z = self.q[:index,z_index]

      self.lines[i].set_data_3d(x, y, z)

      # point is the "head of the line"
      self.points[i].set_data_3d(x[-1:], y[-1:], z[-1:])

    return self.lines + self.points

  def init(self):
    for line, point in zip(self.lines, self.points):
      line.set_data_3d([], [], [])
      point.set_data_3d([], [], [])

    return self.lines + self.points

  def animate(self, solver_name, save):
    self.fig = plt.figure(figsize=(8, 8))
    self.axes = self.fig.add_subplot(projection="3d")

    # parameters
    self.axes.set_facecolor((1, 1, 1))
    self.axes.grid(False)
    self.axes.set_xticklabels([])
    self.axes.set_yticklabels([])
    self.axes.set_xlabel("x [a.u]")
    self.axes.set_ylabel("y [a.u]")
    self.axes.set_zlabel("z [a.u]")

    #self.time_text = self.axes.text()

    self.lines = sum([
      self.axes.plot(
        [], [], [],
        '-',
        color=body.color,
        linewidth=2,
        markevery=10000,
        markersize=body.markersize,
        markerfacecolor=body.color,
        label=body.name
      ) for body in self.bodies
    ], [])

    self.points = sum([
      self.axes.plot(
        [], [], [],
        'o',
        color=body.color,
        linewidth=2,
        markevery=10000,
        markersize=body.markersize,
        markerfacecolor=body.color,
        label=body.name
      ) for body in self.bodies
    ], [])

    self.axes.legend()

    # solve for that specific solver
    solver = self.solvers2[solver_name]
    self.q, self.p, energy = self.solve(solver=solver["call"], dt=self.dt, nt=self.nt, bodies=self.bodies)

    max_range = self.limit_plot(self.q)

    # limiting plot
    self.axes.set_xlim([-max_range,max_range])
    self.axes.set_ylim([-max_range,max_range])
    self.axes.set_zlim([-max_range,max_range])

    ani = animation.FuncAnimation(
      fig=self.fig,
      func=self.update,
      init_func=self.init,
      frames=4000,
      interval=5,
      blit=True
    )

    if save:
      ani.save(filename=f"../n-body-{solver_name}.gif", writer="imagemagick", fps=60)

    plt.show()

  def limit_plot(self, q):
    max_range = 0

    for (ind, body) in enumerate(self.bodies):
      x_index = (ind * 3)
      y_index = (ind * 3) + 1
      z_index = (ind * 3) + 2

      x, y, z = q[:,x_index], q[:,y_index], q[:,z_index]

      max_dim = max(max(x), max(y), max(z))
      if max_dim > max_range:
        max_range = max_dim

    return max_range


