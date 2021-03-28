import numpy as np

from consts import (M_sun as m1, M_jup as m2, M_sat as m3, G)

# two-body problem
# H = p_1^2 / 2m_1 + p_2^2 / 2m_2 + Gm_1m_2 / |r_1 - r_2|
# v1 = [[x1, y1, z1], [p1x, p1y, p1z]]
"""
  :state: state vector => state = [p1, r1, p2, r2]
  :p1: 3D impulsion of Sun => p1 = [p1x, p1y, p1z]
  :r1: 3D position of Sun => r1 = [x1, y1, z1]
  :p2: 3D impulsion of Jupiter => p2 = [p2x, p2y, p2z]
  :r2: 3D position of Jupiter => r2 = [x2, y2, z2]
"""
# def two_body(state):
#   p1, r1, p2, r2 = state[0], state[1], state[2], state[3]

#   r12 = np.sqrt(((r1[0] - r2[0]) ** 2 + (r1[1] - r2[1]) ** 2 + (r1[2] - r2[2]) ** 2))
#   return np.array([
#     [
#       ((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0]),
#       ((G * m1 * m2) / r12 ** 3) * (r1[1] - r2[1]),
#       ((G * m1 * m2) / r12 ** 3) * (r1[2] - r2[2])
#     ],
#     [
#       p1[0] / m1,
#       p1[1] / m1,
#       p1[2] / m1
#     ],
#     [
#       - ((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0]),
#       - ((G * m1 * m2) / r12 ** 3) * (r1[1] - r2[1]),
#       - ((G * m1 * m2) / r12 ** 3) * (r1[2] - r2[2])
#     ],
#     [
#       p2[0] / m2,
#       p2[1] / m2,
#       p2[2] / m2
#     ]
#   ])

def two_body_dhdr(q_state, p_state):
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


def n_body_dhdr(q_state, p_state, bodies):
  """
  :param q_state: n-body position state vector => q_state: [x1, y1, z1, x2, y2, z2,..., xN, yN, zN]
  :param p_state: n-body impulsion state vector => p_state: [px1, py1, pz1,..., pxN, pyN, pzN]
  :param bodies: set of bodies (Sun, Jupiter,...)
  :type q_state: ndarray
  :type p_state: ndarray
  :type bodies: ndarray

  :return: dh/dr
  :rtype: ndarray
  """
  new_state = np.zeros(len(p_state))
  for i in range(len(bodies)):
    # px_i / m_i
    new_state[i] = q_state[i] / bodies[i].mass
    # py_i / m_i
    new_state[i + 1] = q_state[i + 1] / bodies[i].mass
    # pz_i / m_i
    new_state[i + 2] = q_state[i + 1] / bodies[i].mass

  return new_state


def two_body_dhdp(q_state, p_state):
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
      - ((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0]),
      - ((G * m1 * m2) / r12 ** 3) * (r1[1] - r2[1]),
      - ((G * m1 * m2) / r12 ** 3) * (r1[2] - r2[2])
    ]
  ])

def n_body_dhdp(q_state, p_state, bodies):
  """
  :param q_state: n-body position state vector => q_state: [x1, y1, z1, x2, y2, z2,..., xN, yN, zN]
  :param p_state: n-body impulsion state vector => p_state: [px1, py1, pz1,..., pxN, pyN, pzN]
  :param bodies: set of bodies (Sun, Jupiter,...)
  :type q_state: ndarray
  :type p_state: ndarray
  :type bodies: ndarray

  :return: dh/dr
  :rtype: ndarray
  """
  new_state = np.zeros(len(q_state))
  for i in range(len(bodies)):
    pass

  return new_state

# three-body problem
# H = p_1^2 / 2m_1 + p_2^2 / 2m_2 + p_3^2 / 2m_3 + Gm_1m_2 / |r_1 - r_2| + Gm_2m_3 / |r_2 - r_3| + Gm_1m_3 / |r_1 - r_2|
def three_body():
  pass
