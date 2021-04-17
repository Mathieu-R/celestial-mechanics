import numpy as np
from consts import (M_sun as m1, M_jup as m2, M_sat as m3, G)

def n_body_dqdt(qk, pk, bodies):
  """
  :param qk: n-body position state vector => qk: [x1, y1, z1, x2, y2, z2,..., xN, yN, zN]
  :param pk: n-body impulsion state vector => pk: [px1, py1, pz1,..., pxN, pyN, pzN]
  :param bodies: set of bodies (Sun, Jupiter,...)
  :type q_state: ndarray
  :type p_state: ndarray
  :type bodies: ndarray

  :return: \dot{q} = dh/dp
  :rtype: ndarray
  """
  # reshape impulsions vector as a list of 3D vectors
  pk = pk.reshape([-1, 3])
  dqdt = np.zeros(pk.shape)
  for i, pi in enumerate(pk):
    dqdt[i] = pi / bodies[i].mass

  # return the unstructured array of the same shape as before
  return dqdt.flatten()

def n_body_dpdt(qk, pk, bodies):
  """
  :param qk: n-body position state vector at k * dt time => qk: [x1, y1, z1, x2, y2, z2,..., xN, yN, zN]
  :param pk: n-body impulsion state vector at k * dt time => pk: [px1, py1, pz1,..., pxN, pyN, pzN]
  :param bodies: set of bodies (Sun, Jupiter,...)
  :type qk: ndarray
  :type pk: ndarray
  :type bodies: ndarray

  :return: \dot{p} = - dh/dr => \dot{p} = [\dot{px1}, \dot{py2}, \dot{pz1},..., \dot{pxN}, \dot{pyN}, \dot{pzN}]
  :rtype: ndarray
  """
  # reshape positions vector as a list of 3D vectors
  qk = qk.reshape([-1, 3])
  dpdt = np.zeros(qk.shape)
  # loop each body
  for i in range(len(bodies)):
    for j in range(len(bodies)):
      if i != j:
        # compute distance between mass i and mass j
        r_ij = r_dist(qk[i], qk[j])
        dpdt[i] -= ((G * bodies[i].mass * bodies[j].mass) / r_ij ** 3) * (qk[i] - qk[j])

  # return unstructured array
  return dpdt.flatten()

def r_dist(ri, rj):
  """
  Calculate distance between two points in \R^3
  :param ri: point of \R^3: [xi, yi, zi]
  :param rj: point of \R^3: [xj, yj, zj]

  :return: distance betwen the two points
  :rtype: float
  """
  return np.sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2 + (ri[2] - rj[2]) ** 2)
