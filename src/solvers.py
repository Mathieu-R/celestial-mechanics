import numpy as np
import tqdm

def rk2_derivatives_dqdt(edo, qk, pk, dt, bodies):
  k1 = dt * edo(qk, pk, bodies)
  k2 = dt * edo(qk + (dt * k1), pk, bodies)
  return (k1 + k2) / 2.

def rk2_derivatives_dpdt(edo, qk, pk, dt, bodies):
  k1 = dt * edo(qk, pk, bodies)
  k2 = dt * edo(qk, pk + (dt * k1), bodies)
  return (k1 + k2) / 2.

def heun(dqdt, dpdt, q, p, dt, nt, bodies):
  for k in range(0, nt - 1):
    qk, pk = q[k], p[k]
    print("rk2 -- iteration:", k, qk)
    q[k + 1] = qk + rk2_derivatives_dqdt(dqdt, qk, pk, dt, bodies)
    p[k + 1] = pk + rk2_derivatives_dpdt(dpdt, qk, pk, dt, bodies)

  return q, p


def rk4_derivatives_dqdt(edo, qk, pk, dt, bodies):
  k1 = dt * edo(qk, pk, bodies)
  k2 = dt * edo(qk + ((dt / 2) * k1), pk, bodies)
  k3 = dt * edo(qk + ((dt / 2) * k2), pk, bodies)
  k4 = dt * edo(qk + (dt * k3), pk, bodies)
  return (k1 + 2*k2 + 2*k3 + k4) / 6.

def rk4_derivatives_dpdt(edo, qk, pk, dt, bodies):
  k1 = dt * edo(qk, pk, bodies)
  k2 = dt * edo(qk, pk + ((dt / 2) * k1), bodies)
  k3 = dt * edo(qk, pk + ((dt / 2) * k2), bodies)
  k4 = dt * edo(qk, pk + (dt * k3), bodies)
  return (k1 + 2*k2 + 2*k3 + k4) / 6.

def rk4(dqdt, dpdt, q, p, dt, nt, bodies):
  for k in range(0, nt - 1):
    qk, pk = q[k], p[k]
    print("rk4 -- iteration:", k, qk)
    q[k + 1] = qk + rk4_derivatives_dqdt(dqdt, qk, pk, dt, bodies)
    p[k + 1] = pk + rk4_derivatives_dpdt(dpdt, qk, pk, dt, bodies)

  return q, p


def euler_symp(dqdt, dpdt, q, p, dt, nt, bodies):
  for k in range(0, nt - 1):
    qk, pk = q[k], p[k]
    print("euler -- iteration:", k, qk)
    p[k + 1] = pk + dt * dpdt(qk, pk, bodies)
    q[k + 1] = qk + dt * dqdt(qk, p[k + 1], bodies)

  return q, p


def stormer_verlet(dqdt, dpdt, q, p, dt, nt, bodies):
  for k in range(0, nt - 1):
    qk, pk = q[k], p[k]
    #print("stormer-verlet -- iteration:", k, qk)
    p_half = pk + ((dt / 2) * dpdt(qk, pk, bodies))
    q[k + 1] = qk + (dt * dqdt(qk, p_half, bodies))
    p[k + 1] = p_half + ((dt / 2) * dpdt(q[k + 1], p_half, bodies))

  return q, p
