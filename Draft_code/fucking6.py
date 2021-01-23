from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from scipy.constants import speed_of_light, Boltzmann, Avogadro
style.use('ggplot')

flux_wavelength_data = np.load('planet_atmosphere_spectrum_data.npy')
wavelengths, flux = flux_wavelength_data.T
_, flux_standard_deviation = np.load('standard_deviation_planet_spectrum.npy').T
print(flux_standard_deviation)

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



## FUNCTIONS FOR RANGES
def sigma_range(name_of_gas, T_lower=150, T_upper=450):
    """
    This is the standard deviation of the spectral line, due to Doppler effect.
    Don't confuse this with sigma_n, the noise in the measurement of the flux. Loaded from file.

    Returns sigma_range, shape = (len(lines), 30)
    """
    T = np.linspace(T_lower, T_upper, 30)
    k = boltzmann_constant
    c = speed_of_light
    m = gases[name_of_gas]['weight'] * amu_to_kg
    lambda_center = np.array(gases[name_of_gas]['lines'])
    sigma_range = np.zeros((len(lambda_center), 30))
    for line in xrange(len(lambda_center)):
        sigma_range[line] = lambda_center[line] * np.sqrt(k * T / m) / c

    return sigma_range


def flux_min_range():
    """
    Range of flux if no absorption

    Returns flux_min_range, shape = (len(lines), 30)

    """

    flux_min_range = np.linspace(0.7, 1, 30)

    return flux_min_range


def lambda_center_range(name_of_gas):
    """
    This is the range of wavelenghts due to satellite velocity causing doppler shift

    Returns lambda_center_range
    """
    v_max = 10000
    c = speed_of_light
    lambda_center = gases[name_of_gas]['lines']
    lambda_center_range = np.zeros((len(lambda_center), 300))
    for line in xrange(len(lambda_center)):
        max_doppler_shift = lambda_center[line] * v_max / c
        wl_within_max_shift = np.abs(wavelengths - lambda_center[line]) < max_doppler_shift
        lambda_center_range[line] = np.linspace(wavelengths[wl_within_max_shift][0], wavelengths[wl_within_max_shift][-1], 300)

    return lambda_center_range


## MASKING INTERESTING AREA
def get_mask(lambda_center_range):
    mask1 = wavelengths < lambda_center_range[-1]
    mask2 = wavelengths > lambda_center_range[0]
    return np.logical_and(mask1, mask2)






## CALCULATE F_MODEL TERM IN XHI SQUARE
def flux_model(lambda_values, sigma, lambda_center, flux_min):
    flux_max = 1
    flux_model = flux_max + (flux_min - flux_max) * np.exp(-(lambda_values - lambda_center)**2 / (2 * sigma**2))

    return flux_model




## XHI SQUARE MINIMIZATION
def chi_square(sigma_range, lambda_center_range, flux_min_range):
    chi_square = np.zeros((len(sigma_range), len(sigma_range[0]), len(lambda_center_range[0]), len(flux_min_range)))

    for line in xrange(len(sigma_range)):
        print('line:', line)
        mask = get_mask(lambda_center_range[line])
        lambda_values = wavelengths[mask]
        for sigma in xrange(len(sigma_range[line])):
            print('sigma:', sigma)
            for lambda_center in xrange(len(lambda_center_range[line])):
                if lambda_center % 20 == 0:
                    print('lambda center:', lambda_center)
                for f_min in xrange(len(flux_min_range)):
                    f_model = flux_model(lambda_values, sigma_range[line, sigma], lambda_center_range[line, lambda_center], f_min)
                    chi_values = (flux[mask] - f_model)**2 / flux_standard_deviation[mask]**2
                    chi_square[line, sigma, lambda_center, f_min] = np.sum(chi_values)
    """shape = lines, sigma, lambda_center"""
    return chi_square


"""
chi_square = chi_square(
                sigma_range('CO'),
                lambda_center_range('CO'),
                flux_min_range()
            )
"""
#np.save('chi_CO2', chi_square(sigma_range('CO2'), lambda_center_range('CO2'), flux_min_range()))
