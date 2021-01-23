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
star_radius = star_system.star_radius                                            # [Km]
planet_semimajor_axes = star_system.a                                            # [AU]
star_temperature = star_system.temperature


stefan_boltzmann_sigma = 5.670373e-8                                             # [W * m**-2 * K**-4]
astrounits_to_meters = 1.496e11



def calculate_irradiance(planet):
    star_luminosity = 4 * np.pi * (star_radius * 1000)**2 * stefan_boltzmann_sigma * star_temperature**4
    irradiance = star_luminosity / (4 * np.pi * (planet_semimajor_axes[planet] * astrounits_to_meters)**2)

    return irradiance


def calculate_solar_cell_area(irradiance, power_output=40, efficiency=0.12):
    """Lander needs 40W during daytime, panels have efficiency 12%"""
    area_needed = power_output / (irradiance * efficiency)

    return area_needed


def find_planet_temperature(planet):
    planet_temperature = ((star_radius**2 * star_temperature**4) / 4 *(planet_semimajor_axes[planet] * astrounits_to_meters)**2)**0.25

    return planet_temperature
