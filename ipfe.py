from utils import *

q = 28871271685163
B = 100 #Bound of the values of x and y
t = 1 #Size of the polynomial
n = 2 #Size of the vector
Bx = 10 #Bound of x
By = 10 #Bound of y
K = 2*B*n*B*n + 1
delta = q // K

def setup():
    #a = sample_uniform_vector(t, q)
    a = random.randint(0, q - 1)
    s = sample_bounded_vector(n, 1)
    e = sample_bounded_vector(n, 1)
    b = vector_vector_addition(scalar_vector_multiplication(a, s, q), e, q)
    mpk = (a, b)
    return (s, mpk)

def encrypt(mpk, x):
    (a, b) = mpk
    #u = sample_bounded_vector(t, 1)
    u = random.randint(-1, 1)
    c0 = a*u % q
    mu = scalar_vector_multiplication(delta, x, q)
    c1 = vector_vector_addition(scalar_vector_multiplication(u, b, q), mu, q)
    return (c0, c1)

def keygen(msk, y):
    dk = vector_vector_inner_product(msk, y, q)
    return dk

def decrypt(dk, y, c):
    (c0, c1) = c
    res1 = vector_vector_inner_product(y, c1, q)
    res0 = dk*c0 % q
    res = res1 - res0 % q
    return res // delta

y = sample_uniform_vector(n, B)
x = sample_uniform_vector(n, B)
(msk, mpk) = setup()
c = encrypt(mpk, x)
dk = keygen(msk, y)
dec = decrypt(dk, y, c)
res = vector_vector_inner_product(x, y, q)
print(res == dec)