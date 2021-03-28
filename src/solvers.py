import numpy as np
import tqdm

def rk2_derivatives(edo, qk, pk, dt):
  k1 = dt * edo(qk, pk)
  k2 = dt * edo(qk + (dt * k1), pk + (dt * k1))
  return (k1 + k2) / 2

def heun(dhdr, dhdp, q, p, dt, nt):
  for k in range(0, nt - 1):
    pk, qk = p[k], q[k]
    q[k + 1] = qk + rk2_derivatives(dhdr, qk, pk, dt)
    p[k + 1] = pk + rk2_derivatives(dhdp, qk, pk, dt)

  #print(q)
  return q, p


def rk4_derivatives(edo, qk, pk, dt):
  k1 = edo(qk, pk)
  k2 = edo(qk + ((dt / 2) * k1), pk + ((dt / 2) * k1))
  k3 = edo(qk + ((dt / 2) * k2), pk + ((dt / 2) * k2))
  k4 = edo(qk + (dt * k3), pk + (dt * k3))
  return (k1 + 2*k2 + 2*k3 + k4)

def rk4(dhdr, dhdp, q, p, dt, nt):
  for k in range(0, nt - 1):
    q[k + 1] = q[k] + (dt / 6) * rk4_derivatives(dhdr, q[k], p[k], dt)
    p[k + 1] = p[k] + (dt / 6) * rk4_derivatives(dhdp, q[k], p[k], dt)

  return q, p


def euler_symp(dhdr, dhdp, q, p, dt, nt):
  for k in range(0, nt - 1):
    q[k + 1] = q[k] + dt * dhdp(q[k], p[k + 1])
    p[k + 1] = p[k] - dt * dhdr(q[k], p[k])

  return q, p


def stormer_verlet(dhdr, dhdp, q, p, dt, nt):
  for k in range(0, nt - 1):
    p_half = p[k] - (dt / 2) * dhdr(q[k], p[k])
    q[k+1] = q[k] + dt * dhdp(q[k], p_half)
    p[k + 1] = p_half - (dt / 2) * dhdr(q[k + 1], p_half)

  return q, p
