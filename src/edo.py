import numpy as np
from consts import (M_sun as m1, M_jup as m2, M_sat as m3, G)

def two_body_dqdt(q_state, p_state):
  """

  """
  p1, p2 = p_state[0], p_state[1]
  #print(f"p2/m2: {p2[0] / m2}, {p2[1] / m2}, {p2[2] / m2}")

  return np.array([
    [
      p1[0] / m1,
      p1[1] / m1,
      p1[2] / m1
    ],
    [
      p2[0] / m2,
      p2[1] / m2,
      p2[2] / m2
    ]
  ])

def two_body_dpdt(q_state, p_state):
  r1, r2 = q_state[0], q_state[1]
  #print(r1, r2)

  r12 = np.sqrt((r1[0] - r2[0]) ** 2 + (r1[1] - r2[1]) ** 2 + (r1[2] - r2[2]) ** 2)
  #print(f"r12: {r12}")
  #print(f"G: {G}")
  #print(f"m1: {m1}, m2: {m2}, x1 - x2: {r1[0] - r2[0]}")
  #print(f"x1: {((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0])}")
  return np.array([
    [
      ((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0]),
      ((G * m1 * m2) / r12 ** 3) * (r1[1] - r2[1]),
      ((G * m1 * m2) / r12 ** 3) * (r1[2] - r2[2])
    ],
    [
      ((G * m1 * m2) / r12 ** 3) * (r2[0] - r1[0]),
      ((G * m1 * m2) / r12 ** 3) * (r2[1] - r1[1]),
      ((G * m1 * m2) / r12 ** 3) * (r2[2] - r1[2])
    ]
  ])

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
  #print(f"dq/dt -- Jupiter: x:{qk[3]}, y:{qk[4]}, z:{qk[5]}")
  #print(f"dq/dt -- Jupiter: px:{pk[3]}, py:{pk[4]}, pz:{pk[5]}")
  #print(f"Jupiter Mass: {bodies[1].mass}")
  dqdt = np.zeros(len(qk))
  for i in range(len(bodies)):
    offset = i * 3
    # px_i / m_i
    dqdt[offset] = pk[offset] / bodies[i].mass
    # py_i / m_i
    dqdt[offset + 1] = pk[offset + 1] / bodies[i].mass
    # pz_i / m_i
    dqdt[offset + 2] = pk[offset + 2] / bodies[i].mass

  return dqdt

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
  #print(f"dp/dt -- Jupiter: x:{qk[3]}, y:{qk[4]}, z:{qk[5]}")
  #print(f"dq/dt -- Jupiter: px:{pk[3]}, py:{pk[4]}, pz:{pk[5]}")
  dpdt = np.zeros(len(qk))
  # loop each body
  for i in range(len(bodies)):
    i_offset = i * 3
    for j in range(len(bodies)):
      if i != j:
        j_offset = j * 3

        #print(f"i: {i}, i_offset: {i_offset}, j: {j}, j_offset: {j_offset}")
        #print(f"G: {G}, Mass 1: {bodies[i].mass}, Mass 2: {bodies[j].mass}")

        r_ij = r_dist(ri=[qk[i_offset], qk[i_offset + 1], qk[i_offset + 2]], rj=[qk[j_offset], qk[j_offset + 1], qk[j_offset + 2]])

        r_ij2 = r_dist(ri=qk[i_offset:i_offset+3], rj=qk[j_offset:j_offset+3])

        r_ij3 = np.linalg.norm([qk[i_offset:i_offset+3], qk[j_offset:j_offset+3]])

        dpdt[i_offset] -= ((G * bodies[i].mass * bodies[j].mass) / r_ij ** 3) * (qk[i_offset] - qk[j_offset])

        dpdt[i_offset + 1] -= ((G * bodies[i].mass * bodies[j].mass) / r_ij ** 3) * (qk[i_offset + 1] - qk[j_offset + 1])

        dpdt[i_offset + 2] -= ((G * bodies[i].mass * bodies[j].mass) / r_ij ** 3) * (qk[i_offset + 2] - qk[j_offset + 2])

  return dpdt

def r_dist(ri, rj):
  """
  Calculate distance between two points in \R^3
  :param ri: point of \R^3: [xi, yi, zi]
  :param rj: point of \R^3: [xj, yj, zj]

  :return: distance betwen the two points
  :rtype: float
  """
  return np.sqrt((ri[0] - rj[0]) ** 2 + (ri[1] - rj[1]) ** 2 + (ri[2] - rj[2]) ** 2)
