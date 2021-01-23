from __future__ import division
import numpy as np

def trapezoidal(f, a, b, n):
    h = float(b - a) / n
    s = 0.0
    s += f(a) / 2.0
    for i in xrange(1, n):
        s += f(a + i*h)
    s += f(b) / 2.0
    return s * h

def midpoint(f, a, b, n):
    h = (b-a)/n
    Int = 0.0
    x = a + 0.5 * h
    for k in range(1, n+1):
        Int += f(x)
        x += h
    return (h * Int)
