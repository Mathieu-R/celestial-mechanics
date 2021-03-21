#!/usr/bin/env python3.9
import click

from src.two_body import TwoBody

twobody = TwoBody()

@click.command()
@click.option("--stype", default="twobody", help="2-Body problem or 3-Body problem")
@click.option("--plot", default="static", help="static or animated plot")
def main(stype, plot):
  if stype=="twobody":
    if plot=="static":
      twobody.plot()

if __name__ == "__main__":
  main()
