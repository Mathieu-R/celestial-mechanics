import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from mpl_toolkits.mplot3d import Axes3D

from .edo import (n_body_dqdt, n_body_dpdt)
from .solvers import (heun, rk4, euler_symp, stormer_verlet)

class NBodySimulation():
  def __init__(self, bodies, t0, tN, dt):
    self.bodies = bodies
    self.t0 = t0
    self.tN = tN
    self.dt = dt
    self.legends = ["Heun (RK2)", "RK4", "Euler Symplectique", "Stormer-Verlet"]

    self.solvers = [
      {"call": heun, "name": "Heun (RK2)"},
      {"call": rk4, "name": "RK4"},
      {"call": euler_symp, "name": "Euler Symplectique"},
      {"call": stormer_verlet, "name": "Stormer Verlet"}
    ]

    self.results = []

  def solve(self, solver):
    # number of time step
    nt = int((self.tN - self.t0) / self.dt)

    # positions and impulsions state vectors
    # each body is in \R^3 (3D space x, y, z)
    q = np.zeros((nt, len(self.bodies) * 3))
    p = np.zeros((nt, len(self.bodies) * 3))

    # set initial conditions
    q[0] = np.concatenate(np.array([body.initial_positions for body in self.bodies]))
    p[0] = np.concatenate(np.array([body.initial_impulsions for body in self.bodies]))

    return solver(dqdt=n_body_dqdt, dpdt=n_body_dpdt, q=q, p=p, dt=self.dt, nt=nt, bodies=self.bodies)

  def simulate(self):
    self.results = []
    for solver in self.solvers:
      q, p = self.solve(solver["call"])
      self.results.append({"solver": solver["name"], "q": q, "p": p})

  def plot(self):
    self.fig = plt.figure(figsize=(8,8))
    colors = ['r', 'b', 'g', 'y', 'm', 'c']

    # loop for each result (corresponding to a specific solving method)
    for (index, result) in enumerate(self.results):
      print(result["solver"], result["q"][-1])
      # create a 3D plot
      ax = self.fig.add_subplot(2, 2, index + 1, projection="3d")
      # set plot parameters
      ax.set_title(result["solver"])

      max_range = 0

      # plotting each body
      for (ind, body) in enumerate(self.bodies):
        x_index = (ind * 3)
        y_index = (ind * 3) + 1
        z_index = (ind * 3) + 2

        x, y, z = result["q"][:,x_index], result["q"][:,y_index], result["q"][:,z_index]

        max_dim = max(max(x), max(y), max(z))
        #print(max(x), max(y), max(z))
        if max_dim > max_range:
          max_range = max_dim

        ax.plot(xs=x, ys=y, zs=z, c=random.choice(colors), label=body.name)

      # limiting plot
      ax.set_xlim([-max_range,max_range])
      ax.set_ylim([-max_range,max_range])
      ax.set_zlim([-max_range,max_range])

      #ax.legend()

    plt.show()

  def animate(self):
    pass
