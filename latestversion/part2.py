from __future__ import division
import matplotlib.pyplot as plt
from ast2000solarsystem_27_v4 import AST2000SolarSystem
import numpy as np

solar_system = AST2000SolarSystem(11466)

gravitational_constant = 4 * np.pi**2                               #[AU^3 * Yr^(-2) * solar masses^(-1)]
number_of_planets = solar_system.number_of_planets
semi_major_axes = solar_system.a                                    #[AU]
star_mass = solar_system.star_mass                                  #[solar masses]
star_radius = solar_system.star_radius                              #[km]
planet_masses = solar_system.mass                                   #[solar masses]
planet_radii = solar_system.radius                                  #[km]

home_planet_orbital_period = np.sqrt((semi_major_axes[0]**3) / (star_mass + planet_masses[0]))
total_time = 50 #Years
time_steps_per_year = 70000
total_time_steps = time_steps_per_year*total_time
dt = total_time/total_time_steps

times = np.linspace(0, total_time, total_time_steps)

planet_positions = np.zeros((2, number_of_planets, total_time_steps))         #[AU]
planet_velocities = np.zeros((2, number_of_planets, total_time_steps))        #[AU/Yr]
planet_accelerations = np.zeros((2, number_of_planets, total_time_steps))
distance_to_star = np.zeros(number_of_planets)

for p in range(number_of_planets):
    planet_positions[:, p, 0] = [solar_system.x0[p], solar_system.y0[p]]     #[AU]
    planet_velocities[:, p, 0] = [solar_system.vx0[p], solar_system.vy0[p]]  #[AU/Yr]
    distance_to_star[p] = np.sqrt(planet_positions[0, p, 0]**2 + planet_positions[1, p, 0]**2)
    planet_accelerations[:, p, 0] = -gravitational_constant*star_mass * planet_positions[:, 0, 0]/distance_to_star[p]**3


for t in range(total_time_steps-1):
        planet_positions[:, :, t+1] = planet_positions[:, :, t] + planet_velocities[:, :, t]*dt + 0.5*planet_accelerations[:, :, t]*dt**2
        r = np.sqrt(planet_positions[0, :, t] ** 2 + planet_positions[1, :, t] ** 2)
        planet_accelerations[:, :, t+1] = -gravitational_constant * star_mass * planet_positions[:, :, t+1] / r ** 3
        planet_velocities[:, :, t+1] = planet_velocities[:, :, t] + 0.5*(planet_accelerations[:, :, t] + planet_accelerations[:, :, t+1])*dt
        if t%70000 == 0:
            print t/70000

def analytical_orbit(planet):
    analytical_position = np.zeros(total_time_steps)
    analytical_position[0] = np.sqrt(solar_system.x0[planet]**2 + solar_system.y0[planet]**2)
    theta = np.linspace(0, 2*np.pi, total_time_steps)

    a = solar_system.a[planet]
    e = solar_system.e[planet]
    f = theta - solar_system.psi[planet]

    for t in range(total_time_steps-1):
        analytical_position[t+1] = a*(1-e**2)/(1+e*np.cos(f[t+1]))

    return analytical_position*np.cos(theta), analytical_position*np.sin(theta)

def integrate_area(planet, theta1, theta2):
    delta_t = total_time_steps/100  #1/100th of the orbit
    #theta1 = solar_system.psi[planet] -

#solar_system.check_planet_positions(planet_positions, total_time, time_steps_per_year)
#solar_system.orbit_xml(planet_positions, times)
plt.plot(planet_positions[0, 0, :], planet_positions[1, 0, :])
#plt.plot(planet_positions[0, 1, :], planet_positions[1, 1, :])
#plt.plot(planet_positions[0, 2, :], planet_positions[1, 2, :])
#plt.plot(planet_positions[0, 3, :], planet_positions[1, 3, :])
#plt.plot(planet_positions[0, 4, :], planet_positions[1, 4, :])
plt.plot(planet_positions[0, 5, :], planet_positions[1, 5, :])
#plt.plot(planet_positions[0, 6, :], planet_positions[1, 6, :])
#plt.legend(["Home Planet", "Target Planet", "Gas Giant" ])
plt.legend(["Home Planet", "Planet 5"])
plt.xlabel("AU")
plt.ylabel("AU")
plt.plot(analytical_orbit(0)[0], analytical_orbit(0)[1])
plt.plot(analytical_orbit(5)[0], analytical_orbit(5)[1])
#plt.plot(analytical_orbit(6)[0], analytical_orbit(6)[1])
plt.savefig("planet_orbits_with_analytical.png")
plt.show()