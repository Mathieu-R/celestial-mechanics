#!/usr/bin/env python3.9
import click
from click.decorators import option
from click.types import Choice
import matplotlib.pyplot as plt
import pyfiglet # ascii art

from src.nbody import NBodySimulation
from src.body import Body

from consts import (sun_position0, sun_impulsion0, jupiter_position0, earth_position0, jupiter_impulsion0, saturn_position0, saturn_impulsion0, earth_impulsion0, M_sun, M_jup, M_sat, M_earth)
from consts import (t0, tN, dt)

plt.style.use("science")

# prepare bodies
Sun = Body(name="Sun", initial_positions=sun_position0, initial_impulsions=sun_impulsion0, mass=M_sun, color='gold', marker="o", marker_anim="o", markersize=9)
Jupiter = Body(name="Jupiter", initial_positions=jupiter_position0, initial_impulsions=jupiter_impulsion0, mass=M_jup, color='r', marker=",", marker_anim="o-", markersize=8)
Saturn = Body(name="Saturn", initial_positions=saturn_position0, initial_impulsions=saturn_impulsion0, mass=M_sat, color='sandybrown', marker=",", marker_anim="o-", markersize=7)
Earth = Body(name="Earth", initial_positions=earth_position0, initial_impulsions=earth_impulsion0, mass=M_earth, color='dodgerblue', marker=",", marker_anim="o-", markersize=7)

bodies = {
  "Sun": Sun,
  "Jupiter": Jupiter,
  "Saturn": Saturn,
  "Earth": Earth
}

@click.command()
@click.option(
  "--body", "-cb",
  type=click.Choice(["Sun", "Jupiter", "Saturn", "Earth"], case_sensitive=False),
  multiple=True,
  default=["Sun", "Jupiter"],
  show_default=True,
  help="Celestial body to simulate. You can add multiple bodies typing multiples -cb"
)
@click.option(
  "--dimensions", "-d",
  type=int,
  default=2,
  show_default=True,
  help="Only for static plot. Dimensions of the orbital evolution static plot (2D / 3D)"
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
  help="Only needed for animation. Type of numerical scheme. Euler symplectic and Stormer-Verlet are symplectic schemes. For static plot, all the scheme are used at once to compute 4 subplots."
)
@click.option(
  "--time", "-t",
  type=int,
  multiple=True,
  default=[5000],
  show_default=True,
  help="Time of integration (in years). You can add multiples times typing multiples -t"
)
@click.option(
  "--save", "-sa",
  is_flag=True,
  help="Save plot or animation"
)
def main(body, dimensions, plot, solver, time, save):
  # clear terminal (even history)
  print('\033c', end=None)
  # ascii art - for fun.
  print(pyfiglet.print_figlet("CELESTIAL"))

  options = {"save": save}

  # possibility to perform multiple simulations
  # if multiples times given
  for t in time:
    nbody = NBodySimulation(
      bodies=[bodies[b] for b in body],
      t0=t0,
      tN=t * 365.25,
      dt=dt,
      options=options
    )

    if plot == "static":
      nbody.simulate()

      if dimensions == 2:
        nbody.plot2D()
      elif dimensions == 3:
        nbody.plot3D()

      nbody.plot_energy()
      nbody.plot_angular_momentum()
      nbody.plot_area_swept()
    elif plot == "animated":
      nbody.animate(solver)

if __name__ == "__main__":
  main()
