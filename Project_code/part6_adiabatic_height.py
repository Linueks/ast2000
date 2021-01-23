from __future__ import division, print_function
from ast2000solarsystem_27_v6 import AST2000SolarSystem
from scipy.constants import Boltzmann, m_p, gravitational_constant
import numpy as np

solmass_to_kg = 1.9891e30
mu = 17.03
temperature = 271

kg_to_earthmass = 1.673360107e-25


myss = AST2000SolarSystem(11466)
driddu_mass = myss.mass[1] * solmass_to_kg
driddu_radius = myss.radius[1] * 1000
gravity = gravitational_constant * driddu_mass / (driddu_radius)**2
print(gravity)
rho_0 = myss.rho0[1]
scale_height = Boltzmann * temperature / (m_p * mu * gravity)
gamma = 1.4


r_t0 = scale_height * gamma / (2*(gamma-1))


print(r_t0)
