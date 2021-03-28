import numpy as np
import tqdm

def rk2_derivatives(edo, qk, pk, dt, args):
  k1 = dt * edo(qk, pk, args)
  k2 = dt * edo(qk + (dt * k1), pk + (dt * k1), args)
  return (k1 + k2) / 2

def heun(dqdt, dpdt, q, p, dt, nt, args):
  for k in range(0, nt - 1):
    pk, qk = p[k], q[k]
    q[k + 1] = qk + rk2_derivatives(dqdt, qk, pk, dt, args)
    p[k + 1] = pk + rk2_derivatives(dpdt, qk, pk, dt, args)

  #print(q)
  return q, p


def rk4_derivatives(edo, qk, pk, dt, args):
  k1 = edo(qk, pk, args)
  k2 = edo(qk + ((dt / 2) * k1), pk + ((dt / 2) * k1), args)
  k3 = edo(qk + ((dt / 2) * k2), pk + ((dt / 2) * k2), args)
  k4 = edo(qk + (dt * k3), pk + (dt * k3), args)
  return (k1 + 2*k2 + 2*k3 + k4)

def rk4(dqdt, dpdt, q, p, dt, nt, args):
  for k in range(0, nt - 1):
    q[k + 1] = q[k] + (dt / 6) * rk4_derivatives(dqdt, q[k], p[k], dt, args)
    p[k + 1] = p[k] + (dt / 6) * rk4_derivatives(dqdt, q[k], p[k], dt, args)

  return q, p


def euler_symp(dqdt, dpdt, q, p, dt, nt, args):
  for k in range(0, nt - 1):
    q[k + 1] = q[k] + dt * dqdt(q[k], p[k + 1], args)
    p[k + 1] = p[k] + dt * dpdt(q[k], p[k], args)

  return q, p


def stormer_verlet(dqdt, dpdt, q, p, dt, nt, args):
  for k in range(0, nt - 1):
    p_half = p[k] + (dt / 2) * dpdt(q[k], p[k], args)
    q[k+1] = q[k] + dt * dqdt(q[k], p_half, args)
    p[k + 1] = p_half + (dt / 2) * dpdt(q[k + 1], p_half, args)

  return q, p
