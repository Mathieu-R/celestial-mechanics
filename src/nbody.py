from consts import DATA_SUB_INTERVAL_LENGTH
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
    # number of time step
    self.nt = int((self.tN - self.t0) / self.dt)
    self.legends = ["Heun (RK2)", "RK4", "Euler Symplectique", "Stormer-Verlet"]

    self.solvers = [
      {"call": heun, "name": "Heun (RK2)", "bodies": self.bodies},
      {"call": rk4, "name": "RK4", "bodies": self.bodies},
      {"call": euler_symp, "name": "Euler Symplectique", "bodies": self.bodies},
      {"call": stormer_verlet, "name": "Stormer Verlet", "bodies": self.bodies}
    ]

    self.solvers2 = {
      "rk2": {
        "call": heun,
        "name": "Heun (RK2)",
        "bodies": self.bodies
      },
      "rk4": {
        "call": rk4,
        "name": "RK4",
        "bodies": self.bodies
      },
      "euler-sympectic": {
        "call": euler_symp,
        "name": "Euler Symplectique",
        "bodies": self.bodies
      },
      "stormer-verlet": {
        "call": stormer_verlet,
        "name": "Stormer Verlet",
        "bodies": self.bodies
      }
    }

    self.results = []

  def solve(self, solver, dt, nt, bodies):
    # positions and impulsions state vectors
    # each body is in \R^3 (3D space x, y, z)
    q = np.zeros((self.nt, len(self.bodies) * 3))
    p = np.zeros((self.nt, len(self.bodies) * 3))

    # set initial conditions
    q[0] = np.concatenate(np.array([body.initial_positions for body in bodies]))
    p[0] = np.concatenate(np.array([body.initial_impulsions for body in bodies]))

    return solver(dqdt=n_body_dqdt, dpdt=n_body_dpdt, q=q, p=p, dt=dt, nt=nt, bodies=bodies)

  def simulate(self):
    self.results = []
    for solver in self.solvers:
      q, p = self.solve(solver=solver["call"], dt=self.dt, nt=self.nt, bodies=self.bodies)
      self.results.append({"solver": solver["name"], "q": q, "p": p})

  def plot2D(self):
    self.fig, ax = plt.subplots(1, 1, figsize=(8, 8))


    # set plot parameters
    result = self.results[0]
    ax.set_title(result["solver"])
    ax.set_xlabel("x [a.u]")
    ax.set_ylabel("y [a.u]")

    max_range = 0

    # plotting each body
    for (ind, body) in enumerate(self.bodies):
      x_index = (ind * 3)
      y_index = (ind * 3) + 1

      x, y= result["q"][:,x_index], result["q"][:,y_index]

      max_dim = max(max(x), max(y))
      #print(max(x), max(y), max(z))
      if max_dim > max_range:
        max_range = max_dim

      # plot Sun at the center (but should move a little bit over a long time ?)
      if body.name == "Sun":
        ax.plot(0, 0, "o", color="orange", markersize=3)

      ax.plot(x, y, c=body.color, label=body.name)

    # limiting plot
    ax.set_xlim([-max_range,max_range])
    ax.set_ylim([-max_range,max_range])

    #ax.legend()

    plt.show()

  def plot3D(self):
    self.fig = plt.figure(figsize=(8,8))

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

        ax.plot(xs=x, ys=y, zs=z, c=body.color, label=body.name)

      # limiting plot
      ax.set_xlim([-max_range,max_range])
      ax.set_ylim([-max_range,max_range])
      ax.set_zlim([-max_range,max_range])

      #ax.legend()

    plt.show()

  def update(self, index, solver_name):
    print(index, solver_name)
    if index > self.tN:
      return

    sub_interval_t0 = (index - 1) * DATA_SUB_INTERVAL_LENGTH
    sub_interval_tk = (index) * DATA_SUB_INTERVAL_LENGTH

    last_index = int(sub_interval_tk / self.dt) - 1
    print("last index", last_index)

    nt = int((sub_interval_tk - sub_interval_t0) / self.dt)
    print("nt", nt)

    # get the solver
    solver = self.solvers2[solver_name]

    # simulate on the sub-interval
    q, p = self.solve(solver=solver["call"], dt=self.dt, nt=nt, bodies=solver["bodies"])

    for (body_index, body) in enumerate(solver["bodies"]):
      # N bodies : [x1, y1, z1, ..., xN, yN, zN]
      x_index = body_index * 3
      y_index = (body_index * 3) + 1
      z_index = (body_index * 3) + 2

      # update initial conditions for the next sub-interval
      body.initial_positions = q[last_index, x_index:z_index+1]
      body.initial_impulsions = p[last_index, x_index:z_index+1]

      #print(q, p)
      print(body.initial_positions)
      print(body.initial_impulsions)

      # update plot
      body.line.set_xdata(q[0:last_index,x_index])
      body.line.set_ydata(q[0:last_index,y_index])
      # https://stackoverflow.com/questions/46685326/how-to-set-zdata-for-a-line3d-object-in-python-matplotlib
      body.line.set_3d_properties(q[0:last_index,z_index])

  def animate(self, solver_name):
    self.fig = plt.figure(figsize=(8, 8))
    self.axes = self.fig.add_subplot(projection="3d")

    # parameters
    self.axes.set_facecolor((0.5, 0.5, 0.5)) # 50% gray
    self.axes.grid(False)
    self.axes.set_xticklabels([])
    self.axes.set_yticklabels([])
    self.axes.set_xlabel("x [a.u]")
    self.axes.set_xlabel("y [a.u]")
    self.axes.set_xlabel("z [a.u]")

    #self.time_text = self.axes.text()

    for body in self.bodies:
      body.line, = self.axes.plot(
        xs=np.zeros((self.nt)),
        ys=np.zeros((self.nt)),
        zs=np.zeros((self.nt)),
        color=body.color,
        linewidth=2,
        markevery=10000,
        markersize=body.markersize,
        markerfacecolor=body.color,
        label=body.name
      )

    ani = animation.FuncAnimation(
      fig=self.fig,
      func=self.update,
      fargs=(solver_name,),
      frames=range(1, int(self.tN / DATA_SUB_INTERVAL_LENGTH)),
      interval=self.dt,
      repeat=False
    )

    plt.show()


