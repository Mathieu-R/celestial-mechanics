import click

@click.command()
@click.option("--type", default="twobody", help="2-Body problem or 3-Body problem")
@click.option("--plot", default="animated", help="static or animated plot")
def main():
  pass

if __name__ == "__main__":
  main()
