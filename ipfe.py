from utils import *
import random
import numpy as np

q = 2**607 - 1
B = 10 #Bound of the values of x and y entries
t = 4 #Size of the polynomial n sur le papier
n = 5 #Size of the vector l sur le papier
K = n*(B*n)*(B*n) + 1
delta = q // K
P = np.polynomial.Polynomial([1] + [0] * (t-2) + [1])

def setup():
    a = rand_pol(t, q)
    s = [rand_bd_pol(t, 1) for _ in range(n)]
    e = [rand_bd_pol(t, 1) for _ in range(n)]
    b = [a * s[i] + e[i] % P for i in range(n)]
    mpk = (a, b)
    return (s, mpk)

def encrypt(mpk, x):
    (a, b) = mpk
    u = rand_bd_pol(t, 1)
    e0 = rand_bd_pol(t, 1)
    e1 = [rand_bd_pol(t, 1) for _ in range(n)]
    c0 = a * u + e0 % P
    mu = [delta * x[i] * pol_unit(t) for i in range(n)]
    c1 = [u * b[i] + e1[i] + mu[i] % P for i in range(n)]
    return (c0, c1)

def keygen(msk, y):
    dk = sum(y[i] * msk[i] for i in range(n)) % P
    return dk

def decrypt(dk, y, c):
    (c0, c1) = c
    res1 = sum(y[i] * c1[i] for i in range(n)) % P
    res0 = c0 * dk % P
    res = list(res1 - res0)[0]
    return round(float(res) / delta)

y = [random.randint(0, B - 1) for _ in range(n)]
x = [random.randint(0, B - 1) for _ in range(n)]
(msk, mpk) = setup()
c = encrypt(mpk, x)
dk = keygen(msk, y)
dec = decrypt(dk, y, c)
res = vect_ip(x, y, q)
print(res)
print(dec)
