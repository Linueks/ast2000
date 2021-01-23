from __future__ import division, print_function
from scipy.constants import Boltzmann, Avogadro, speed_of_light, m_p, gravitational_constant
from ast2000solarsystem_27_v6 import AST2000SolarSystem
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
style.use('ggplot')

seed = 11466
star_system = AST2000SolarSystem(seed)

# CONVERSIONS
solar_to_kg = 1.9884e30
km_to_m = 1000


k = Boltzmann
n_a = Avogadro
c = speed_of_light
m_p = m_p
G = gravitational_constant



gravity = G * star_system.mass[1] * solar_to_kg / (star_system.radius[1] * km_to_m)**2
surface_temp = 271                                                               #[K]
gamma = 1.4
rho_0 = star_system.rho0[1]
mu = 17.02
scale_height = Boltzmann * surface_temp / (2 * m_p * mu * gravity)
planet_mass = star_system.mass[1] * 1.989e30

print('star mass: ', 2.44*solar_to_kg)
print('star radius: ', (star_system.star_radius*1000))
print('star radius^2: ', (star_system.star_radius*1000)**2)
print('star temp: ', star_system.temperature)
print('star volume: ', 4/3 * np.pi * (star_system.star_radius*1000)**3)
print('gravitational constant: ', gravitational_constant)
print('boltzmann constant: ', Boltzmann)
print('mass hydrogen: ', m_p)
print('mu: ', 1.757)



def compute_adiabatic_analytic(r):
    prefac = ((gamma-1) / (gamma)) * ((G * planet_mass * mu * m_p) / (k))
    T = prefac * ((1/r) - (1/star_system.radius[1]*km_to_m)) + surface_temp
    rho = rho_0 * (T / surface_temp)**(1/(gamma-1))
    return rho


def compute_isothermal_analytic(r_0, r):
    rho_0 = compute_adiabatic_analytic(r_0)
    exponent = ((2 * G * planet_mass * mu * m_p) / (k * surface_temp))
    rho = rho_0 * np.exp(exponent * (1/r - 1/r_0))
    return rho





def compute_r():
    r_inv = (1.0/star_system.radius[1] * km_to_m) - (surface_temp/2) * (gamma / (gamma-1)) * (k / (G * planet_mass * mu * m_p))
    return (r_inv)**(-1)

print('her er rrrrr:', compute_r())


def compute_atmos_numeric():
    P = k * rho_0 * surface_temp / (mu * m_p)
    adiabatic_constant = (P**(1-gamma)) * (surface_temp**gamma)
    r = star_system.radius[1] * km_to_m

    T = surface_temp
    rho_list = [rho_0]
    ideal_gas_constant = mu * m_p / k
    exponent = (gamma-1)/gamma
    mu_gravity = -G * planet_mass
    dr = 1

    while T > surface_temp / 2:
        P += mu_gravity * rho_list[-1] / (r**2) * dr
        T = (adiabatic_constant / (P**(1-gamma)))**(1/gamma)
        rho_list.append(ideal_gas_constant*(P/T))
        r+=dr

    no_points_adiabatic = len(rho_list)
    max_r_adiabatic = r
    max_r_analytic = compute_r()


    while rho_list[-1] > 1e-10:
        P += mu_gravity * rho_list[-1] / r**2 * dr
        rho_list.append(ideal_gas_constant * (P/T))
        r += dr

    r_tot = np.linspace(star_system.radius[1] * km_to_m, r, len(rho_list))
    no_points_iso = len(rho_list) - no_points_adiabatic
    r_iso = np.linspace(max_r_adiabatic, r, no_points_iso)
    r_adiabatic = np.linspace(star_system.radius[1]*km_to_m, max_r_adiabatic, no_points_adiabatic)



    rho_adiabatic_analytic = compute_adiabatic_analytic(r_adiabatic)
    rho_isothermal_analytic = compute_isothermal_analytic(max_r_adiabatic, r_iso)
    rho_tot_analytic = np.append(rho_adiabatic_analytic, rho_isothermal_analytic)

    plt.plot(r_tot - (star_system.radius[1] * km_to_m), rho_list)
    #plt.plot(r_tot - (star_system.radius[1] * km_to_m), np.log(rho_list))
    plt.xlabel('Distance from planet surface [m]')
    plt.ylabel('Density [kg/m**3]')
    plt.legend(['Logarithmic Atmosphere Density Profile'])
    plt.hold('on')

    #plt.plot(r_adiabatic, rho_adiabatic_analytic)
    #plt.plot(r_iso, rho_isothermal_analytic)
    #np.save('rho_values.npy', (rho_list, r_tot))
    plt.show()


#compute_atmos_numeric()


"""
r_1 = np.linspace(star_system.radius[1], compute_r(), 1000000)
r_2 = np.linspace(compute_r(), 300000, 1000000)
print(compute_adiabatic_analytic(r_1))
"""

"""
plt.plot(r_1, compute_adiabatic_analytic(r_1))
plt.hold('on')
plt.plot(r_2, compute_isothermal_analytic(compute_r(), r_2))
plt.show()
"""
