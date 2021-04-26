import numpy as np

class Body():
  """
  Body (Sun, Jupiter,...) in cartesian coordinates
  """
  def __init__(self, name, initial_positions, initial_impulsions, mass, color, markersize):
    self.name = name
    self.initial_positions = np.array(initial_positions)
    self.initial_impulsions = np.array(initial_impulsions)
    self.mass = mass
    self.markersize = markersize

    self.color = color
    #self.line = None
