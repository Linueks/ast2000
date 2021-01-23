from __future__ import division, print_function
import numpy as np
from scipy.constants import speed_of_light, Boltzmann, m_p, gravitational_constant
from ast2000solarsystem_27_v5 import AST2000SolarSystem
import matplotlib.pyplot as plt
import matplotlib.style as style
style.use('ggplot')


seed = 11466
star_system = AST2000SolarSystem(seed)


gravity = gravitational_constant * star_system.mass[1] * 1.9891e30 / (star_system.radius[1] * 1000)**2

scale_height = Boltzmann * 271 / (gravity * 18 * 2 * m_p)
gamma = 1.4

print(scale_height)
print(star_system.mass[1])


def temp_profile(r):
    if r <= scale_height * gamma / (2*(gamma-1)):
        temp = 271 * (1 - (gamma - 1) / gamma * (r/scale_height))
        return temp
    else:
        temp = 271 / 2
        return temp

rho_0 = 7.9

def dens_profile(r):
    if r <= scale_height * gamma / (2 * (gamma - 1)):
        density = rho_0 * (1 - (gamma - 1) / gamma * (r/scale_height))
        return density
    else:
        density = (rho_0 / 2)**(1/(gamma+1))*np.exp(-r/scale_height)
        return density

r = np.linspace(0, 80000, 2000000)
#plt.plot(r, map(temp_profile, r))
plt.plot(r, map(dens_profile, r))
plt.show()




"""
print(star_system.temperature)

lambdas = np.array([719.9975117, 1659.95523846, 2869.91143779])
sigmas = np.array([0.00085463, 0.00154593, 0.00338083])
weights = np.array([18.01528, 16.04, 44.013]) * m_p


temps = weights / Boltzmann * (speed_of_light * sigmas / lambdas)**2

#print(temps)
"""
