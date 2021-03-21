import numpy as np
import tqdm

def rk2_derivatives(edo, v, dt):
  k1 = edo(v)
  k2 = edo(v + (dt * k1))
  return (k1 + k2)

def heun(dhdr, dhdp, q, p, dt, nt):
  for t in range(0, nt - 1):
    q[t + 1] = q[t] + (dt / 2) * rk2_derivatives(dhdr, q[t], dt)
    p[t + 1] = p[t] + (dt / 1) * rk2_derivatives(dhdp, p[t], dt)

  return q, p


def rk4_derivatives(edo, v, dt):
  k1 = edo(v)
  k2 = edo(v + ((dt / 2) * k1))
  k3 = edo(v + ((dt / 2) * k2))
  k4 = edo(v + (dt * k3))
  return (k1 + 2*k2 + 2*k3 + k4)

def rk4(dhdr, dhdp, q, p, dt, nt):
  for t in tqdm(range(0, nt - 1)):
    q[t + 1] = q[t] + (dt / 6) * rk4_derivatives(dhdr, q[t], dt)
    p[t + 1] = p[t] + (dt / 6) * rk4_derivatives(dhdp, p[t], dt)

  return q, p


def euler_symp(dhdr, dhdp, q, p, dt, nt):
  for t in tqdm(range(0, nt - 1)):
    p[t+1] = p[t] + dt * dhdr(q[t])
    q[t+1] = q[t] + dt * dhdp(p[t+1])

  return q, p


def stormer_verlet_it(dhdr, dhdp, q, p, dt, nt):
  for t in tqdm(range(0, nt - 1)):
    pdemi = - (dt / 2) * dhdr(q[t])
    q[t+1] = dt * dhdp(pdemi)
    p[t + 1] = pdemi - (dt / 2) * dhdr(q[t + 1])

  return q, p
