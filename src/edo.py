import numpy as np
from numpy import linalg
from numpy.linalg.linalg import norm
from consts import (M_sun as m1, M_jup as m2, M_sat as m3, G)

"""
Compute total energy of the N-Body system.
"""
def hamiltonian(qk, pk, bodies):
  qk = qk.reshape([-1, 3])
  pk = pk.reshape([-1, 3])

  p_sum = 0
  for i, pi in enumerate(pk):
    p_sum += (pi ** 2) / 2 * bodies[i].mass

  q_sum = 0
  for i in range(len(bodies)):
    for j in range(len(bodies)):
      if i != j:
        # compute distance between mass i and mass j
        r_ij = r_dist(qk[i], qk[j])
        q_sum += ((G * bodies[i].mass * bodies[j].mass) / r_ij)

  # return hamiltonian (= energy)
  return np.linalg.norm(p_sum + q_sum)

"""
Compute total angular momentum of the N-Body system.
"""
def compute_angular_momentum(qk, pk, bodies):
  qk = qk.reshape([-1, 3])
  pk = pk.reshape([-1, 3])

  Lk = np.zeros(3)

  L_sum = 0
  # L(t_j) = \sum_{i=1}^{N} l_i(t_j)
  for i in range(len(bodies)):
    pi, qi = pk[i], qk[i]
    vi = pi / bodies[i].mass

    rdotv = np.dot(qi, vi)
    rnorm = np.linalg.norm(qi)
    vnorm = np.linalg.norm(vi)

    # r \cdot v = ||r|| ||v|| \cos(\theta)
    theta = np.arccos(rdotv / (rnorm * vnorm))
    # L = r x mv = m ||r|| ||v|| \sin(\theta)
    l = bodies[i].mass * rnorm * vnorm * np.sin(theta)
    L_sum += l
    #Lk[i] = L

  return L_sum #Lk

def compute_area_swept(qk, pk, qk_next, pk_next, bodies):
  qk = qk.reshape([-1, 3])
  pk = pk.reshape([-1, 3])

  qk_next = qk_next.reshape([-1, 3])
  pk_next = pk_next.reshape([-1, 3])

  darea_sum = 0
  for i in range(len(bodies)):
    pi, qi = pk[i], qk[i]
    pi_next, qi_next = pk_next[i], qk_next[i]

    # \theta = tan^{-1}(y / x)
    theta = np.arctan(abs(qi[1] / qi[0]))
    theta_next = np.arctan(abs(qi_next[1] / qi_next[0]))

    # r_moyen = (1/2) * (ri_norm + rinext_norm)
    #rmoy = 0.5 * (np.linalg.norm(qi) + np.linalg.norm(qi_next))

    dtheta = np.abs(theta - theta_next)

    # compute delta area between two successive points q[i] and q[i+1]
    # \Delta A = (1/2) ||r_1|| ||r_2|| \Delta theta
    darea = 0.5 * dtheta * (np.linalg.norm(qi) * np.linalg.norm(qi_next))
    darea_sum += darea

  return darea_sum


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
