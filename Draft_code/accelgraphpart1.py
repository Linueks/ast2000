from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

force = 1
mass = np.linspace(0, 100000)

accel = force / mass

plt.plot(mass, accel)
plt.show()
