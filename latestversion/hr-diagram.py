from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
import astropy.units as u
from astropy.constants import k_B, M_sun, L_sun, sigma_sb, R_sun



star_temps = np.linspace(2000, 40000, 100) * u.Kelvin
star_radiuses = np.linspace(0.1, 100, 100) * R_sun


def luminosity(star_temp, star_radius):
    lum = 4 * np.pi * star_radius**2 * sigma_sb * star_temp**4
    return lum

luminosity = luminosity(star_temps, star_radiuses)


plt.scatter(star_temps[::-1], (luminosity / L_sun))
plt.show()
