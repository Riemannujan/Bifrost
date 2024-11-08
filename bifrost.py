from utils import *

def key_gen():
    A = sample_uniform_vector(n, q)
    sk = sample_bounded_vector(n, B)
    As = vector_vector_multiplication(A, sk, q)
    e = scalar_vector_multiplication(t, sample_bounded_vector(n, B), q)
    b = vector_vector_addition(As, e, q)
    pk = (A, b)
    return (sk, pk)

def encrypt(pk, m):
    (A, b) = pk
    u = sample_bounded_vector(n, B)
    e1 = scalar_vector_multiplication(t, sample_bounded_vector(n, B), q)
    e2 = scalar_vector_multiplication(t, sample_bounded_vector(n, B), q)
    Aalpha = vector_vector_multiplication(A, alpha, q)
    c1 = vector_vector_addition(Aalpha, vector_vector_addition(e1, u, q), q)
    balpha = vector_vector_multiplication(b, alpha, q)
    mu = vector_vector_multiplication(m, u, q)
    c2 = vector_vector_addition(balpha, vector_vector_addition(e2, mu, q), q)
    return (c1, c2)

def keyder(pk, sk, v):
    (A, b) = pk
    vs = vector_vector_addition(v, scalar_vector_multiplication(-1, sk, q), q)
    dk = vector_vector_multiplication(A, vector_vector_multiplication(alpha, vs, q), q)
    return dk

def decrypt(dk, v, ct):
    (c1, c2) = ct
    vc = vector_vector_multiplication(v, c1, q)
    dkvc = vector_vector_addition(dk, scalar_vector_multiplication(-1, vc, q), q)
    res = vector_vector_addition(c2, dkvc, q)
    res = [res[i] % q for i in range(n)]
    return [res[i] % t for i in range(n)]

(sk, pk) = key_gen()
m = sample_uniform_vector(n, t)
v = m
ct = encrypt(pk, m)
dk = keyder(pk, sk, v)
res = decrypt(dk, v, ct)
print(res)
