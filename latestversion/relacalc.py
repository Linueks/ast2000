import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from scipy.constants import c
style.use('ggplot')

xs = np.array([0.0, 3.64*10**5, 2.32*10**5, 0.0])
ts = np.array([0.0, 1.34*1e-3, 1.34*1e-3, 2.67*1e-3])
vel = 0.65
gamma = 1 / np.sqrt(1 - vel**2)



def lorentz(x, t):
    x_lor = gamma * x + -vel * gamma * t
    t_lor = -vel * gamma * x + gamma * t

    return x_lor, t_lor

x_lor, t_lor = lorentz(xs, ts)
print('initial x: ', xs)
print('initial t: ', ts)

print('transformed x:', x_lor)
print('transformed t:', t_lor)



plt.subplot(212)
plt.scatter(x_lor, np.abs(t_lor))
plt.title('transformed')

plt.subplot(211)
plt.scatter(xs, ts)
plt.title('before')
plt.show()
