import numpy as np

class Body():
  """
  Body (Sun, Jupiter,...) in cartesian coordinates
  """
  def __init__(self, name, initial_positions, initial_impulsions, mass, color, marker, marker_anim, markersize):
    self.name = name
    self.initial_positions = np.array(initial_positions)
    self.initial_impulsions = np.array(initial_impulsions)
    self.mass = mass
    self.marker = marker # type of marker for static plot
    self.marker_anim = marker_anim # type of marker for animated plot
    self.markersize = markersize

    self.color = color
