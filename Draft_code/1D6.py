from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from scipy.constants import speed_of_light, Boltzmann, Avogadro
style.use('ggplot')


wavelengths, flux = np.load('planet_atmosphere_spectrum_data.npy').T
_, flux_standard_deviation = np.load('standard_deviation_planet_spectrum.npy').T


###PHYSICAL CONSTANTS
speed_of_light = speed_of_light
boltzmann_constant = Boltzmann
avogadro_constant = Avogadro


###UNIT CONVERSION
amu_to_kg = 1.66054e-27



gases= {
        'O2':{'lines': [630, 690, 760], 'weight': 32},
        'CO':{'lines': [2340], 'weight': 28},
        'H2O':{'lines': [720, 820, 940], 'weight': 18},
        'CO2':{'lines': [1400, 1600], 'weight': 44},
        'CH4':{'lines': [1660, 2200], 'weight': 16},
        'N2O':{'lines': [2870], 'weight':30}
}


def sigma_range(name_of_gas, T_lower=150, T_upper=450):
    """
    This is the standard deviation of the spectral line, due to Doppler effect.
    Don't confuse this with sigma_n, the noise in the measurement of the flux. Loaded from file.
    """
    T = np.linspace(T_lower, T_upper, 30)
    k = boltzmann_constant
    c = speed_of_light
    m = gases[name_of_gas]['weight'] * amu_to_kg
    lambda_center = np.reshape(np.array(gases[name_of_gas]['lines']), (3,1))

    sigma_range = lambda_center / c * np.sqrt((k*T/m))
    return sigma_range


def flux_min_range():
    flux_min_range = np.linspace(0.7, 1, 30)
    return flux_min_range


def lambda_center_range(name_of_gas):
    """
    This is the range of wavelenghts due to satellite velocity causing doppler shift
    """
    v_max = 10000
    c = speed_of_light
    lambda_center = gases[name_of_gas]['lines']
    lambda_center_range = np.zeros((len(lambda_center), 300))
    for line in xrange(len(lambda_center))
        max_doppler_shift = lambda_center * v_max / c
        wl_within_max_shift = np.abs(wavelengths - lambda_center) < max_doppler_shift

        lambda_center_range = np.linspace(wavelengths[wl_within_max_shift][0], wl[wavelengths_within_max_shift][-1], 300)
    return lambda_center_range


def flux_model(sigma_range, flux_min_range, lambda_center_range):
    flux_model = flux_min_range[-1] + (flux_min_range - flux_min_range[-1]) * np.exp(-(wavelengths - lambda_center_range)**2 / (2 * sigma_range**2))

    return flux_model

print(flux_model(sigma_range('O2'), flux_min_range(), lambda_center_range('O2')))
