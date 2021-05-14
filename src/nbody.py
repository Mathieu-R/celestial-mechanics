import os
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
      "heun": {
        "call": heun,
        "name": "Heun (RK2)",
        "color": "teal",
        "bodies": self.bodies
      },
      "euler-symplectic": {
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

    script_dir = os.path.dirname(__file__)
    self.figures_dir = os.path.join(script_dir, f"../report/figures/{int(self.tN / 365.25)}_years")

    print(script_dir, self.figures_dir)

    # create figures folder if does not exist
    if not os.path.isdir(self.figures_dir):
      os.makedirs(self.figures_dir)

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
    angular_momentum = np.zeros(self.nt)
    # total area swept mesh
    area_swept = np.zeros(self.nt)

    # set initial conditions
    q[0] = np.concatenate(np.array([body.initial_positions for body in bodies]))
    p[0] = np.concatenate(np.array([body.initial_impulsions for body in bodies]))

    # set initial energy
    #print(energy)
    #print(hamiltonian(qk=q[0], pk=p[0], bodies=bodies))
    energy[0] = hamiltonian(qk=q[0], pk=p[0], bodies=bodies)
    angular_momentum[0] = compute_angular_momentum(qk=q[0], pk=p[0], bodies=bodies)
    area_swept[0] = 0

    return solver(dqdt=n_body_dqdt, dpdt=n_body_dpdt, q=q, p=p, dt=dt, nt=nt, bodies=bodies, energy=energy, angular_momentum=angular_momentum, area_swept=area_swept)

  def simulate(self):
    self.results = []
    for solver in self.solvers:
      q, p, energy, angular_momentum, area_swept = self.solve(solver=solver["call"], dt=self.dt, nt=self.nt, bodies=self.bodies)
      self.results.append({"solver": solver["name"], "color": solver["color"], "q": q, "p": p, "energy": energy, "angular_momentum": angular_momentum, "area_swept": area_swept})

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
      self.fig.savefig(f"{self.figures_dir}/orbital-plot2d.png", dpi=600)

    plt.show()

  def plot3D(self):
    self.fig = plt.figure(figsize=(set_size_square_plot(width="column-size", subplots=(2,2))))
    # loop for each result (corresponding to a specific solving method)
    for (index, result) in enumerate(self.results):
      # create a 3D plot
      ax = self.fig.add_subplot(2, 2, index + 1, projection="3d")
      # set plot parameters
      ax.set_title(result["solver"], fontsize=8)

      # plotting each body
      for (ind, body) in enumerate(self.bodies):
        x_index = (ind * 3)
        y_index = (ind * 3) + 1
        z_index = (ind * 3) + 2

        x, y, z = result["q"][:,x_index], result["q"][:,y_index], result["q"][:,z_index]
        ax.plot(xs=x, ys=y, zs=z, c=body.color, marker=body.marker, label=body.name)

      # set labels
      ax.set_xlabel("x (AU)", fontsize=6)
      ax.set_ylabel("y (AU)", fontsize=6)
      ax.set_zlabel("z (AU)", fontsize=6)

    plt.tight_layout()

    if self.options["save"]:
      self.fig.savefig(f"{self.figures_dir}/orbital-plot3d.png", dpi=600)

    plt.show()

  def plot_energy(self):
    self.fig = plt.figure(figsize=(set_size(width="full-size")))
    ax = self.fig.add_subplot(1, 1, 1)

    ax.set_title(f"Energie du système à {len(self.bodies)}-corps", fontsize=10)

    for (index, result) in enumerate(self.results):
      #solver_name = result["solver"]
      # kg * au^2 / day^2 => kg * m^2 / s^2
      energy = result["energy"] * ((au_to_meter ** 2) / (day_to_second ** 2))
      ax.plot(self.time_mesh / 365.25, energy, c=result["color"])
      ax.set_xlabel("Temps [années]")
      ax.set_ylabel("Energie totale [$J$]")


    if self.options["save"]:
      self.fig.savefig(f"{self.figures_dir}/orbital-energy.pdf")

    plt.show()

  def plot_angular_momentum(self):
    self.fig = plt.figure(figsize=(set_size(width="full-size")))
    ax = self.fig.add_subplot(1, 1, 1)

    ax.set_title(f"Moment angulaire du système à {len(self.bodies)}-corps", fontsize=10)

    for (index, result) in enumerate(self.results):
      #ax = self.fig.add_subplot(2, 2, index + 1)
      #solver_name = result["solver"]

      # kg * au^2 / day => kg * m^2 / s
      angular_momentum = result["angular_momentum"] * ((au_to_meter ** 2) / (day_to_second))
      ax.plot(self.time_mesh / 365.25, angular_momentum, c=result["color"], label=result["solver"])

      #for (ind, body) in enumerate(self.bodies):
      #  ax.plot(self.time_mesh / 365.25, angular_momentum[:,ind], c=body.color, label=body.name)

      ax.tick_params(axis='both', which='major', labelsize=8)
      ax.tick_params(axis='both', which='minor', labelsize=6)

      ax.set_xlabel("Temps [années]", fontsize=8)
      ax.set_ylabel("L [$kg.m^2.s^{-1}$]", fontsize=8)
      #ax.text(0.5, 1.05, f"({index + 1})", ha="center", transform=ax.transAxes, size=8)
      #ax.set_title(f"{solver_name}", fontsize=8)

    plt.tight_layout()

    if self.options["save"]:
      self.fig.savefig(f"{self.figures_dir}/orbital-angular-momentum.pdf")

    plt.show()

  def plot_area_swept(self):
    self.fig = plt.figure(figsize=(set_size(width="full-size")))
    ax = self.fig.add_subplot(1, 1, 1)

    ax.set_title(f"Surface totale balayée pour le système à {len(self.bodies)}-corps", fontsize=10)

    for (index, result) in enumerate(self.results):
      area_swept = result["area_swept"]
      ax.plot(self.time_mesh / 365.25, area_swept, c=result["color"], label=result["solver"])

      #ax.tick_params(axis='both', which='major', labelsize=8)
      #ax.tick_params(axis='both', which='minor', labelsize=6)

      ax.set_xlabel("Temps [années]", fontsize=8)
      ax.set_ylabel("Aire totale balayée [$(AU)^2$]", fontsize=8)

    plt.tight_layout()

    if self.options["save"]:
      self.fig.savefig(f"{self.figures_dir}/orbital-total_area_swept.pdf")

    plt.show()

  def update(self, index):
    for i in range(len(self.bodies)):
      x_index = (i * 3)
      y_index = (i * 3) + 1
      z_index = (i * 3) + 2

      # show a portion of the path
      trail = 40
      beginning = max(1, index - trail)

      x = self.q[index:beginning:-1,x_index]
      y = self.q[index:beginning:-1,y_index]
      z = self.q[index:beginning:-1,z_index]

      # update elapsed time
      #print(f"Temps écoulé: {round(self.time_mesh[index] / 365.25, 1)} ans")
      self.elapsed_text.set_text(f"Temps écoulé: {round(self.time_mesh[index] / 365.25, 1)} ans")

      self.lines[i].set_data_3d(x, y, z)

      # point is the "head of the line"
      #self.points[i].set_data_3d(x[-1:], y[-1:], z[-1:])

    return tuple(self.lines) + (self.elapsed_text,) #+ self.points

  def init(self):
    for line in self.lines: #line, point in zip(self.lines, self.points):
      line.set_data_3d([], [], [])
      #point.set_data_3d([], [], [])

    self.elapsed_text.set_text("")

    return tuple(self.lines) + (self.elapsed_text,) #+ self.points

  def animate(self, solver_name):
    self.fig = plt.figure(figsize=(8, 8))
    self.axes = self.fig.add_subplot(projection="3d")

    # parameters
    self.axes.set_facecolor((1, 1, 1))
    self.axes.grid(False)
    self.axes.set_xticklabels([])
    self.axes.set_yticklabels([])
    self.axes.set_zticklabels([])
    self.axes.set_xlabel("x (AU)")
    self.axes.set_ylabel("y (AU)")
    self.axes.set_zlabel("z (AU)")

    self.lines = sum([
      self.axes.plot(
        [], [], [],
        body.marker_anim,
        #'-',
        color=body.color,
        linewidth=2,
        markevery=10000,
        markersize=body.markersize,
        markerfacecolor=body.color,
        label=body.name
      ) for body in self.bodies
    ], [])

    # self.points = sum([
    #   self.axes.plot(
    #     [], [], [],
    #     'o',
    #     color=body.color,
    #     linewidth=2,
    #     markevery=10000,
    #     markersize=body.markersize,
    #     markerfacecolor=body.color,
    #     label=body.name
    #   ) for body in self.bodies
    # ], [])

    self.axes.legend()

    # solve for that specific solver
    solver = self.solvers2[solver_name]
    self.q, self.p, energy, angular_momentum, area_swept = self.solve(solver=solver["call"], dt=self.dt, nt=self.nt, bodies=self.bodies)

    max_range = self.limit_plot(self.q)

    # limiting plot
    self.axes.set_xlim([-max_range,max_range])
    self.axes.set_ylim([-max_range,max_range])
    self.axes.set_zlim([-max_range,max_range])

    # text, annotations,...
    self.elapsed_text = self.axes.text2D(
      x=0.5,
      y=0.85,
      s='',
      bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
      transform = self.axes.transAxes,
      ha='center'
    )

    ani = animation.FuncAnimation(
      fig=self.fig,
      func=self.update,
      init_func=self.init,
      frames=4000,
      interval=5,
      blit=True
    )

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


