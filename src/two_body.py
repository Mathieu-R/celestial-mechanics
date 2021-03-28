import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from mpl_toolkits.mplot3d import Axes3D

from .edo import (two_body_dhdr, two_body_dhdp)
from .solvers import (heun, rk4, euler_symp, stormer_verlet)

from consts import (sun_position0, sun_impulsion0, jupiter_position0, jupiter_impulsion0)
from consts import (t0, tN, dt)

# I think that this command, also remove precision
# in real array in memory not only in the print
#np.set_printoptions(precision=3, suppress=True)

class TwoBody():
  def __init__(self):
    self.initial_positions = [sun_position0, jupiter_position0]
    self.initial_impulsions = [sun_impulsion0, jupiter_impulsion0]
    self.t0 = t0
    self.tN = tN
    self.dt = dt
    self.legends = ["Heun (RK2)", "RK4", "Euler Symplectique", "Stormer-Verlet (Symplectique)"]

  def solve(self):
    nt = int((self.tN - self.t0) / self.dt)

    time_mesh = np.linspace(start=self.t0, stop=self.tN, num=nt)

    # positions and impulsions state vectors
    q = np.ones((nt, 2, 3))
    p = np.ones((nt, 2, 3))

    # set initial conditions
    q[0] = self.initial_positions
    p[0] = self.initial_impulsions

    q_heun, p_heun = heun(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q.copy(), p=p.copy(), dt=self.dt, nt=nt)
    q_rk4, p_rk4 = rk4(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q.copy(), p=p.copy(), dt=self.dt, nt=nt)
    q_euler, p_euler = euler_symp(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q.copy(), p=p.copy(), dt=self.dt, nt=nt)
    q_stormer, p_stormer = stormer_verlet(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q.copy(), p=p.copy(), dt=self.dt, nt=nt)

    return q_heun, q_rk4, q_euler, q_stormer

  def plot(self):
    q_heun, q_rk4, q_euler, q_stormer = self.solve()

    self.fig = plt.figure(figsize=(8,8))
    self.heun_axes = self.fig.add_subplot(2, 2, 1, projection="3d")
    self.rk4_axes = self.fig.add_subplot(2, 2, 2, projection="3d")
    self.euler_axes = self.fig.add_subplot(2, 2, 3, projection="3d")
    self.stormer_axes = self.fig.add_subplot(2, 2, 4, projection="3d")

    # Heun (RK2)
    #self.axes.plot(xs=q_heun[:,0,0], ys=q_heun[:,0,1], zs=q_heun[:,0,2])
    self.heun_axes.plot(xs=q_heun[:,1,0], ys=q_heun[:,1,1], zs=q_heun[:,1,2])
    self.heun_axes.set_title("Heun (RK2)")

    # RK4
    self.rk4_axes.plot(xs=q_rk4[:,0,0], ys=q_rk4[:,0,1], zs=q_rk4[:,0,2])
    self.rk4_axes.plot(xs=q_rk4[:,1,0], ys=q_rk4[:,1,1], zs=q_rk4[:,1,2])
    self.rk4_axes.set_title("RK4")

    # Euler Symplectic
    self.euler_axes.plot(xs=q_euler[:,0,0], ys=q_euler[:,0,1], zs=q_euler[:,0,2])
    self.euler_axes.plot(xs=q_euler[:,1,0], ys=q_euler[:,1,1], zs=q_euler[:,1,2])
    self.euler_axes.set_title("Euler Symplectique")

    # Stormer Verlet (symplectic)
    self.stormer_axes.plot(xs=q_stormer[:,0,0], ys=q_stormer[:,0,1], zs=q_stormer[:,0,2])
    self.stormer_axes.plot(xs=q_stormer[:,1,0], ys=q_stormer[:,1,1], zs=q_stormer[:,1,2])
    self.stormer_axes.set_title("Stormer Verlet (Symplectique)")

    plt.tight_layout()
    plt.show()

  def animate(self):
    self.animation = animation.FuncAnimation(
      fig=self.fig,
      func=self.live_orbit,
      #frames=,
      repeat=True
    )

  def live_orbit(self):
    pass
