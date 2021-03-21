import numpy as np

from astropy.constants import (M_sun as m1, M_jup as m2, G)
from consts import (M_sat as m3)

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
def two_body(state):
  p1, r1, p2, r2 = state[0], state[1], state[2], state[3]
  
  r12 = ((r1[0] - r2[0]) ** 2 + (r1[1] - r2[1]) ** 2 + (r1[2] - r2[2]) ** 2)
  return np.array([
    [
      p1[0] / m1,
      p1[1] / m1,
      p1[2] / m1
    ],
    [
      ((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0]),
      ((G * m1 * m2) / r12 ** 3) * (r1[1] - r2[1])
      ((G * m1 * m2) / r12 ** 3) * (r1[2] - r2[2])
    ],
    [
      p2[0] / m2,
      p2[1] / m2,
      p2[2] / m2
    ],
    [
      - ((G * m1 * m2) / r12 ** 3) * (r1[0] - r2[0]),
      - ((G * m1 * m2) / r12 ** 3) * (r1[1] - r2[1])
      - ((G * m1 * m2) / r12 ** 3) * (r1[2] - r2[2])
    ]
  ])
  
  
# three-body problem
# H = p_1^2 / 2m_1 + p_2^2 / 2m_2 + p_3^2 / 2m_3 + Gm_1m_2 / |r_1 - r_2| + Gm_2m_3 / |r_2 - r_3| + Gm_1m_3 / |r_1 - r_2|
def three_body():
  pass