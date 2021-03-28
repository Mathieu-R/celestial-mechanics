class Body():
  """
  Body (Sun, Jupiter,...) in cartesian coordinates
  """
  def __init__(self, name, initial_positions, initial_impulsions, mass):
    self.name = name
    self.initial_positions = initial_positions
    self.initial_impulsions = initial_impulsions
    self.mass = mass
