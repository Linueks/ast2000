from __future__ import division, print_function
from scipy.constants import c, Boltzmann, m_p, m_n
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
import math as math

style.use('ggplot')
tolerance = 1e-8

planet_spectrum_data = np.load('planet_atmosphere_spectrum_data.npy')
planet_spectrum_standard_deviation = np.load('standard_deviation_planet_spectrum.npy')

boltzmann_constant = Boltzmann
mass_proton = m_p
mass_neutron = m_n
speed_of_light = c


wavelengths, flux = planet_spectrum_data.T
standard_deviation = planet_spectrum_standard_deviation.T[1]
noisy_wavelengths = wavelengths + standard_deviation

gases= {
        'O2':{'lines': [630, 690, 760], 'weight': 32},
        'CO':{'lines': [2340], 'weight': 28},
        'H2O':{'lines': [720, 820, 940], 'weight': 18},
        'CO2':{'lines': [1400, 1600], 'weight': 44},
        'CH4':{'lines': [1660, 2200], 'weight': 16},
        'N2O':{'lines': [2870], 'weight':30}
}


def get_standard_deviation_range(gas, spectral_line, T_min=150, T_max=450, resolution=40):
    k = boltzmann_constant
    m = gases[gas]['weight'] * (mass_proton + mass_neutron) / 2
    lambd = spectral_line
    temps = np.linspace(T_min, T_max, resolution)
    standard_deviation = lambd / c * np.sqrt((k*temps/m))

    return standard_deviation


def get_flux_range(resolution=40):
    flux_min_range = np.linspace(flux.min() / flux.max(), flux.max() / flux.max(), resolution)
    ## flux_min_range = np.linspace(flux.min(), flux.max(), resolution)


    return flux_min_range


def get_wavenlength_center_range(gas, spectral_line, width=None, res=40):
    satellite_velocity = 10000  # [m/s]
    max_doppler_shift = spectral_line * satellite_velocity / c

    if not width:
        width = max_doppler_shift

    mask = distance_from_center_wavelength_mask(wavelenghts, spectral_line, width)
    wavelength = wavelengths[mask]
    wavelength_center_range = np.linspace(wavelength[0], wavelength[-1], resolution)

    return wavelength_center_range


def distance_from_center_wavelength_mask(wavelenghts, center, width):
    mask = np.abs(wavelenghts - center) < width
    return mask


def find_nearest_index(array,value):
    index = np.searchsorted(array, value, side="left")
    if index > 0 and (index == len(array) or math.fabs(value - array[index-1]) < math.fabs(value - array[index])):
        return index-1
    else:
        return index


def f_model(wavelength, sigma, flux_min, spectral_line):
    flux_max = 1
    return flux_max - (flux_max, flux_min) * np.exp(-(wavelength - spectral_line)**2 / (2 * sigma**2))

"""
for gas in gases:
    for x_coord in gases[gas]['lines']:
        plt.axvline(x = x_coord, color='b')
"""

def plot_spectrum():
    for gas in sorted(gases.keys()):
        info = gases[gas]
        lines = info['lines']

        fig, axes = plt.subplots(len(lines))

        try:
            axes[0]
        except TypeError:
            axes = [axes]

        for i, line in enumerate(lines):
            mask = distance_from_center_wavelength_mask(wavelengths, center=line, width=0.1)
            axes[i].plot(wavelengths[mask], flux[mask])
            plt.show()


plt.xlabel('wavelengths [nm]')
plt.ylabel('Flux')
plot_spectrum()
#plt.plot(wavelengths, flux)
#plt.savefig('wavelenghts_flux.png')
