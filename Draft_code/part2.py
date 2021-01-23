from __future__ import division
from ast2000solarsystem_27_v4 import AST2000SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
style.use('ggplot')
seed = 11466


gravitational_constant = 4 * np.pi**2                                            # [AU**3 * yr**-2 * M_solar**-1]

star_system = AST2000SolarSystem(seed)
star_mass = star_system.star_mass                                                # [Solar mass]
star_radius = star_system.star_radius                                            # [km]


planet_masses = star_system.mass                                                 # [Solar mass]
planet_radii = star_system.radius                                                # [km]
number_of_planets = star_system.number_of_planets
planet_semimajor_axes = star_system.a                                            # [AU]
planet_eccentricities = star_system.e
planet_orbit_angle = star_system.psi


time_steps_per_year = 1000
home_planet_orbit_test = np.sqrt(planet_semimajor_axes[0]**3)

def get_planet_orbital_period(planet):
    planet_orbital_period = np.sqrt(planet_semimajor_axes[planet]**3 / (planet_masses[planet] + star_mass))
    return planet_orbital_period



years_home_planet = 1
total_time = years_home_planet * get_planet_orbital_period(0)
total_time_steps = int(round(time_steps_per_year * total_time))
dt = get_planet_orbital_period(0) / 1000

times = np.linspace(0, total_time, total_time_steps)


planet_positions = np.zeros((2, number_of_planets, total_time_steps))                                                                           # [AU]
planet_velocities = np.zeros((2, number_of_planets, total_time_steps))                                                                          # [AU/Yr]
planet_accelerations = np.zeros((2, number_of_planets, total_time_steps))                                                                       # [AU/Yr^2]
distance_to_star = np.zeros(number_of_planets)

for p in xrange(number_of_planets):
    planet_positions[:, p, 0] = [star_system.x0[p], star_system.y0[p]]                                                                          # [AU]
    planet_velocities[:, p, 0] = [star_system.vx0[p], star_system.vy0[p]]                                                                       # [AU/Yr]
    distance_to_star[p] = np.sqrt(planet_positions[0, p, 0]**2 + planet_positions[1, p, 0]**2)                                                  # [AU]
    planet_accelerations[:, p, 0] = -gravitational_constant * star_mass * planet_positions[:, 0, 0] / distance_to_star[p]**3


for t in xrange(total_time_steps-1):
    planet_positions[:, :, t+1] = planet_positions[:, :, t] + planet_velocities[:, :, t] * dt + 0.5 * planet_accelerations[:, :, t] * dt**2
    r = np.sqrt(planet_positions[0, :, t]**2 + planet_positions[1, :, t]**2)
    planet_accelerations[:, :, t+1] = -gravitational_constant * star_mass * planet_positions[:, :, t+1] / r**3
    planet_velocities[:, :, t+1] = planet_velocities[:, :, t] + 0.5 * (planet_accelerations[:, :, t+1] + planet_accelerations[:, :, t]) * dt



def analytic_orbit(planet, resolution):
    #print planet_eccentricities, planet index 5 has the largest eccentricity

    analytical_position = np.zeros(resolution)
    analytical_position[0] = np.sqrt(planet_positions[0, planet, 0]**2 + planet_positions[1, planet, 0]**2)
    polar_angle = np.linspace(0, 2 * np.pi, resolution)

    a = planet_semimajor_axes[planet]
    e = planet_eccentricities[planet]
    f = polar_angle - planet_orbit_angle[planet]

    for i in xrange(resolution-1):
        analytical_position[i+1] = a * (1 - e**2) / (1 - e*np.cos(f[i+1]))

    x = np.cos(polar_angle) * analytical_position
    y = np.sin(polar_angle) * analytical_position

    return x, y



def check_kepler_third():
    """What i need to do is to integrate the positions of a planet
     over some time at one point in the orbit, and then integrate
     the positions over the same amount of time another point in the orbit
     Checking the law for planet index 5 because the task said so..."""
    orbital_period = get_planet_orbital_period(5)

    placeholder_positions = 0.0
    placeholder_times = np.linspace(0, orbital_period, total_time_steps)

    orbit_x_integral = np.trapz(planet_positions[0, 5, :], axis=0)
    orbit_y_integral = np.trapz(planet_positions[1, 5, :], axis=0)
    #integral_aphelion = np.trapz(placeholder_positions, placeholder_times)
    #integral_perihelion = np.trapz(placeholder_positions, placeholder_times)
    return orbit_x_integral, orbit_y_integral

#x, y = analytic_orbit(5, 1000)
#plt.plot(planet_positions[0, 5, :], planet_positions[1, 5, :])
#plt.plot(x, y)
#plt.show()

print check_kepler_third()












"""dritt"""
