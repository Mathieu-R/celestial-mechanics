import numpy as np
import matplotlib.pyplot as plt

from .edo import (two_body_dhdr, two_body_dhdp)
from .solvers import (heun, rk4, euler_symp, stormer_verlet_it)

from consts import (sun_position0, sun_impulsion0, jupiter_position0, jupiter_impulsion0)
from consts import (t0, tN, dt)

np.set_printoptions(precision=3, suppress=True)

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

    q_heun, p_heun = heun(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q, p=p, dt=self.dt, nt=nt)
    #q_rk4, p_rk4 = rk4(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q, p=p, dt=self.dt, nt=nt)
    #q_euler, p_euler = euler_symp(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q, p=p, dt=self.dt, nt=nt)
    #q_stormer, p_stormer = stormer_verlet_it(dhdr=two_body_dhdr, dhdp=two_body_dhdp, q=q, p=p, dt=self.dt, nt=nt)

    return q_heun#, q_rk4, q_euler, q_stormer

  def plot(self):
    q_heun = self.solve()

    fig, ax = plt.subplots(2,2, figsize=(8,8))

    #Heun (RK2)
    ax[0][0].plot(xs=q_heun[:,0,0], ys=q_heun[:,0,1], zs=q_heun[:,0,2], projection="3d")
    ax[0][0].plot(xs=q_heun[:,1,0], ys=q_heun[:,1,1], zs=q_heun[:,1,2], projection="3d")

    plt.tight_layout()
    plt.show()
