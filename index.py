#!/usr/bin/env python3.9
import click
import matplotlib.pyplot as plt

from src.two_body import TwoBody
from src.nbody import NBodySimulation
from src.body import Body

from consts import (sun_position0, sun_impulsion0, jupiter_position0, jupiter_impulsion0, saturn_position0, saturn_impulsion0, M_sun, M_jup, M_sat)
from consts import (t0, tN, dt)

plt.style.use("science")

# prepare bodies
Sun = Body("Sun", sun_impulsion0, sun_impulsion0, M_sun)
Jupiter = Body("Jupiter", jupiter_position0, jupiter_impulsion0, M_jup)
Saturn = Body("Saturn", saturn_position0, saturn_impulsion0, M_sat)

@click.command()
@click.option("--stype", default="twobody", help="2-Body problem or 3-Body problem")
@click.option("--plot", default="static", help="static or animated plot")
def main(stype, plot):
  if stype == "twobody":
    twobody = TwoBody()
    if plot == "static":
      twobody.plot()
    elif plot == "animated":
      pass
  elif stype == "threebody":
    ThreeBody = NBodySimulation(bodies=[Sun, Jupiter, Saturn], t0=t0, tN=tN, dt=dt)
    ThreeBody.simulate()
    if plot == "static":
      ThreeBody.plot()
    elif plot == "animated":
      ThreeBody.animate()

if __name__ == "__main__":
  main()
