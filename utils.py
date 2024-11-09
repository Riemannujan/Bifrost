import random
import numpy as np

def rand_pol(t, q):
    v = [random.randint(0, q - 1) for _ in range(t)]
    return np.polynomial.Polynomial(v)

def rand_bd_pol(t, B):
    v = [random.randint(-B, B) for _ in range(t)]
    return np.polynomial.Polynomial(v)

def pol_unit(t):
    return np.polynomial.Polynomial([1] * t)

def vect_ip(v, w, q):
    n = len(v)
    return sum(v[i] * w[i] for i in range(n)) % q
