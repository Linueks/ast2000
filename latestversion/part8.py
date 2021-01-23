from __future__ import division, print_function
from ast2000solarsystem_27_v6 import AST2000SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style

seed = 11466
driddu = 1

my_SS = AST2000SolarSystem(seed)
my_SS.part2C_5(number_of_light_signals = 50)
