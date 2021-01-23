from __future__ import division
from ast2000solarsystem_27_v4 import AST2000SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
style.use('fivethirtyeight')
seed = 11466

star_system = AST2000SolarSystem(seed)
gravitational_constant = 4 * np.pi**2                                            # [AU**3 * yr**-2 * M_solar**-1]
planet_mass = star_system.mass                                                   # [Solar masses]
star_mass = star_system.star_mass                                                # [Solar masses]
planet_semimajor_axes = star_system.a                                            # [AU]
planet_number = 6


planet_orbital_period = np.sqrt(planet_semimajor_axes[planet_number]**3 / (star_mass+planet_mass[planet_number]))
time_steps_per_year = 100000
simulation_years = 4
total_time = planet_orbital_period * simulation_years
total_time_steps = time_steps_per_year * simulation_years
dt = total_time / total_time_steps


positions = np.zeros((2, 2, total_time_steps))
velocities = np.zeros((2, 2, total_time_steps))
acceleration = np.zeros((2, 2, total_time_steps))
system_energy = np.zeros(total_time_steps)


positions[:, 0, 0] = [0, 0]
positions[:, 1, 0] = [star_system.x0[planet_number], star_system.y0[planet_number]]
velocities[:, 0, 0] = [0, 0]
velocities[:, 1, 0] = [star_system.vx0[planet_number], star_system.vy0[planet_number]]
acceleration[:, 0, 0] = -gravitational_constant * planet_mass[planet_number] * (positions[:, 0, 0] - positions[:, 1, 0]) / np.linalg.norm(positions[:, 0, 0] - positions[:, 1, 0])**3
acceleration[:, 1, 0] = -gravitational_constant * star_mass * (positions[:, 1, 0] - positions[:, 0, 0]) / np.linalg.norm(positions[:, 1, 0] - positions[:, 0, 0])**3



for t in xrange(total_time_steps-1):
    r_vector = positions[:, 0, t] - positions[:, 1, t]
    r = np.linalg.norm(r_vector)

    #Star
    positions[:, 0, t+1] = positions[:, 0, t] + velocities[:, 0, t] * dt + 0.5 * acceleration[:, 0, t] * dt**2
    acceleration[:, 0, t+1] = -gravitational_constant * planet_mass[planet_number] * r_vector / r**3
    velocities[:, 0, t+1] = velocities[:, 0, t] + 0.5*(acceleration[:, 0, t] + acceleration[:, 0, t+1]) * dt

    #Planet
    positions[:, 1, t+1] = positions[:, 1, t] + velocities[:, 1, t] * dt + 0.5 * acceleration[:, 1, t] * dt ** 2
    acceleration[:, 1, t+1] = -gravitational_constant * star_mass * (-r_vector) / r ** 3
    velocities[:, 1, t+1] = velocities[:, 1, t] + 0.5* (acceleration[:, 1, t] + acceleration[:, 1, t + 1]) * dt




"""
plt.plot(
    positions[0, 0, :], positions[1, 0, :]
)
plt.plot(
    positions[0, 1, :], positions[1, 1, :]
)
"""


# Choosing a reference point far away on the x-axis, the radial velocity just becomes the x-velocity
time = np.linspace(0, total_time, total_time_steps)
y = velocities[0, 0, :]
# using this to simulate "real" values using a Gaussian distribution function to emulate the observational noise
y = y + np.random.normal(loc=0.0, scale=0.2 * np.amax(velocities[0, 0, :]), size=total_time_steps)
plt.xlabel('[Yr]')
plt.ylabel('[AU/yr]')
plt.plot(time, y, linewidth=0.3)
plt.show()
