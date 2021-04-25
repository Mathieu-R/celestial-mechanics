#!/usr/bin/env python3.9
import click
from click.decorators import option
import matplotlib.pyplot as plt
import pyfiglet # ascii art

from src.nbody import NBodySimulation
from src.body import Body

from consts import (sun_position0, sun_impulsion0, jupiter_position0, jupiter_impulsion0, saturn_position0, saturn_impulsion0, M_sun, M_jup, M_sat)
from consts import (t0, tN, dt)

plt.style.use("science")

# prepare bodies
Sun = Body(name="Sun", initial_positions=sun_position0, initial_impulsions=sun_impulsion0, mass=M_sun, color='y', markersize=9)
Jupiter = Body(name="Jupiter", initial_positions=jupiter_position0, initial_impulsions=jupiter_impulsion0, mass=M_jup, color='r', markersize=8)
Saturn = Body(name="Saturn", initial_positions=saturn_position0, initial_impulsions=saturn_impulsion0, mass=M_sat, color='b', markersize=7)

bodies = {
  "Sun": Sun,
  "Jupiter": Jupiter,
  "Saturn": Saturn
}

@click.command()
@click.option(
  "--body", "-cb",
  type=click.Choice(["Sun", "Jupiter", "Saturn"], case_sensitive=False),
  multiple=True,
  default=["Sun", "Jupiter"],
  show_default=True,
  help="Celestial body to simulate. You can add multiple bodies typing multiple times -cb"
)
@click.option(
  "--plot", "-p",
  type=click.Choice(["static", "animated"], case_sensitive=False),
  default="static",
  show_default=True,
  help="(S)tatic or (A)nimated plot"
)
@click.option(
  "--solver", "-s",
  type=click.Choice(["heun", "euler-symplectic", "stormer-verlet"]),
  default="stormer-verlet",
  show_default=True,
  help="Type of numerical scheme. Euler symplectic and Stormer-Verlet are symplectic schemes. Only needed for animation. For static plot, all the scheme are used at once to compute 4 subplots."
)
@click.option(
  "--save", "-sa",
  is_flag=True,
  help="Save plot or animation"
)
def main(body, plot, solver, save):
  # clear terminal (even history)
  print('\033c', end=None)
  # ascii art - for fun.
  print(pyfiglet.print_figlet("CELESTIAL"))

  options = {"save": save}

  nbody = NBodySimulation(
    bodies=[bodies[b] for b in body],
    t0=t0,
    tN=tN,
    dt=dt,
    options=options
  )

  if plot == "static":
    nbody.simulate()
    nbody.plot3D()
    nbody.plot_energy()
    nbody.plot_angular_momentum()
  elif plot == "animated":
    nbody.animate(solver, save)

if __name__ == "__main__":
  main()
