-from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
import scipy as sp
from scipy.constants import speed_of_light, Boltzmann, Avogadro

style.use('ggplot')


## IMPORTING DATA
"""The file contains flux measurements in the range 600nm
to 3000nm normalized such that the background flux is 1"""
planet_flux_measurements = np.load('planet_atmosphere_spectrum_data.npy')

"""A file called sigma_noise.txt with the standard deviation sigma_n of the random
noise fluctuations for each observed wavelength is found in the same directory. """
planet_spectrum_standard_deviation = np.load('standard_deviation_planet_spectrum.npy')

_, standard_deviation = planet_spectrum_standard_deviation.T
wavelengths, flux = planet_flux_measurements.T



## PHYSICAL CONSTANTS
speed_of_light = speed_of_light
boltzmann_constant = Boltzmann
avogadro_constant = Avogadro


## UNIT CONVERSION
gas_weight_to_kg = 1.674e-27


## Dictionary of gases' spectral lines we are trying to identify
gases= {
        'O2':{'lines': [630, 690, 760], 'weight': 32},
        'CO':{'lines': [2340], 'weight': 28},
        'H2O':{'lines': [720, 820, 940], 'weight': 18},
        'CO2':{'lines': [1400, 1600], 'weight': 44},
        'CH4':{'lines': [1660, 2200], 'weight': 16},
        'N2O':{'lines': [2870], 'weight':30}
}


## PROBLEM SPECIFIC CONSTANTS
observation_speed_max = 10000   # [m/s]


def difference_from_spectral_line(wavelengths, name, line, width=0.01):
    mask = np.abs(wavelengths - gases[name][line]) < width
    return mask


def wavelengths_inside_wavelength_range(wavelengths, min_wavelength, max_wavelength):
    mask1 = wavelengths < max_wavelength
    mask2 = wavelengths > min_wavelength
    wavelengths_inside_wavelength_range = np.logical_and(mask1, mask2)

    return wavelengths_inside_wavelength_range


def standard_deviation_of_gaussian_absorption_line(FWHM):
    standard_deviation = FWHM / np.sqrt(8*np.log(2))

    return standard_deviation


def max_deviation_in_wavelength(name, line):
    c = speed_of_light
    for line in gases[name][line]:
        max_change_in_wavelength = observation_speed_max * line / c

    return max_change_in_wavelength


def find_wavelength_range(name, line, width=None, res=300):
    if not width:
        width = max_deviation_in_wavelength(name, line)
    mask = difference_from_spectral_line(wavelengths, name, line, width)
    wavelengths_within_max_shift = wavelengths[mask]
    wavelength_range = np.linspace(wavelengths_within_max_shift[0], wavelengths_within_max_shift[-1], res)
    corresponding_fluxes = flux[mask]
    corresponding_flux_range = np.linspace(corresponding_fluxes[0], corresponding_fluxes[1], res)
    corresponding_noise = standard_deviation[mask]
    return wavelength_range


def flux_range(res=30):
    F_min_min = 0.7
    F_min_max = 1
    F_min = np.linspace(F_min_min, F_min_max, res)
    return F_min


def calculate_sigma(name, line, T_lower = 150, T_upper = 450, res=30):
    N_a = avogadro_constant
    k = boltzmann_constant
    c = speed_of_light

    def compute_FWHM(name, line, T):
        FWHM = (2*gases[name][line])/c * np.sqrt(2*k*T*np.log(2) / (gases[name][weight] * gas_weight_to_kg))
        return FWHM
    return np.linspace(compute_FWHM(name, line, T_lower), compute_FWHM(name, line, T_upper), res)

def f_model(interesting_wavelengths, F_min, sigma, center_wavelength_range):
    F_max = 1
    return F_max - (F_max - F_min) * np.exp(-(interesting_wavelengths - center_wavelength_range)**2 / (2 * sigma**2))


def plot_spectra():
    """
    Looking at plots possible absorption around spectral lines:
    CO2 - 1600nm
    H2O - 720nm
    N2O - 2870nm
    """
    for gas in sorted(gases.keys()):
        gas_info = gases[gas]
        spectral_lines = gas_info['lines']
        print('-----------', gas, spectral_lines, '-----------')

        fig, axes = plt.subplots(len(spectral_lines))

        try:
            axes[0]
        except TypeError:
            axes = [axes]

        for i, line in enumerate(spectral_lines):
            mask = difference_from_spectral_line(wavelengths, gases[], line)
            axes[i].plot(wavelengths[mask], flux[mask])
            axes[i].axvline(x=line,
                            ymin = 0,
                            ymax = 1,
                            c='b',
                            linewidth = 2)
            axes[i].legend(['spectra around %dnm' % line,
                            '%dnm' % line])
            axes[i].set_title(gas+','+str(line))
            axes[i].set_xlabel('wavelenght [nm]')
            axes[i].set_ylabel('Normalized Flux [lumen]')

        plt.show()


def find_lowest_flux_at_lines(gas1='CO2', gas2='H2O', gas3='N2O'):
    return


"""
def lets_try():

    #polyfit(x, y, deg)
    #trenger a bruke polyfit pa interessante bolgelengder med tilhorende flux, per linje.. ma loope over gases og lines i guess

    wavelength_range, corresponding_flux_range = find_wavelength_range('CO', 'lines')
    x = np.linspace(wavelength_range[0], wavelength_range[-1], 10000)
    third, cubic, lin, constant = sp.polyfit(wavelength_range, corresponding_flux_range, 3)
    plt.plot(x, third * x**3 + cubic * x**2 + lin * x + constant)
    plt.plot(wavelength_range, corresponding_flux_range)
    plt.show()
lets_try()
"""


def search_for_light(name = "CO2", line = 1400):
    gas_names = sorted(gases.keys())
    m = [gases[name]['weight'] for name in gas_names]

    F_values =  flux_range()
    l_c_values = find_wavelength_range(name, line)
    sigma_values = calculate_sigma(name, line)

    sigma_max = np.max(sigma_values)
    max_wl = np.max(l_c_values) + 4*sigma_max
    min_wl = np.min(l_c_values) - 4*sigma_max
    mask = get_mask2(wavelengths, min_wl, max_wl)
    lambda_values = wavelengths[mask]
    flux_values = flux[mask]
    noise_values = noise[mask]

    xhi = np.zeros((300, 30, 30))
    data = np.array([l_c_values, F_values, sigma_values])

    for k,sigma in enumerate(sigma_values):
        print(k)
        for j,F_min in enumerate(F_values):
            for i,l_c in enumerate(l_c_values):
                mask = get_mask2( lambda_values, l_c - 4*sigma_max,
                                  l_c + 4*sigma_max)
                lambs = lambda_values[mask]
                fmodl = f_model(lambs, F_min, sigma, l_c)
                val=(spectrum_values[mask]-fmodl)**2/noise_values[mask]**2
                xhi[i,j,k] = np.sum(val)
    return xhi, data




#search_for_light('N2O', 2870)

plot_spectra()





#dritt
