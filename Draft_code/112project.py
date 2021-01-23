from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numerical_integration import trapezoidal
from matplotlib import style
style.use('ggplot')


m = 3.34746e-27         #[kg]
T = 10000               #[K]
k = 1.38064852e-23      #[m**2 * kg * s**-2 * K**-1]
sigma = np.sqrt(k*T/m)
velocity_values = np.linspace(-25000, 25000, 100000)
print(sigma)


def probability_distribution(vx):
    return (1 / (np.sqrt(2*np.pi) * sigma)) * np.exp((-1/2) * (vx**2 / sigma**2))


integral_right_side = trapezoidal(probability_distribution, sigma, 25000, 10000)
integral_left_side = trapezoidal(probability_distribution, -25000, -sigma, 10000)

total_integral = integral_left_side + integral_right_side
print(total_integral)

plt.plot(velocity_values, probability_distribution(velocity_values))
plt.show()
