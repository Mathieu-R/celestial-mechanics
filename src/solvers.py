def heun():
  pass


def rk4_derivatives(edo, v, dt, *args):
  k1 = edo(v, *args)
  k2 = edo(v + ((dt / 2) * k1), *args)
  k3 = edo(v + ((dt / 2) * k2), *args)
  k4 = edo(v + (dt * k3), *args)
  return (k1 + 2*k2 + 2*k3 + k4)

def rk4(edo, v, dt, *args):
  return (dt / 6) * rk4_derivatives(edo, v, dt, *args)

def euler_symp():
  pass

def stormer_verlet():
  pass
  