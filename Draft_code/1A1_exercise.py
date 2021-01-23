from __future__ import division
from numerical_integration import trapezoidal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style

style.use('ggplot')


def prob_density(T):
    standard_deviation = 2
    mean_temp = 29

    density = (1 / (standard_deviation * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((T - mean_temp)**2) / standard_deviation**2)
    return density


print(trapezoidal(prob_density, 30, 100, 100000))
print(trapezoidal(prob_density, -1000, 30, 100000))
