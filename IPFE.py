from utils import *

def IPFE_setup(
    ciphertext_modulus: int,
    poly_modulus: np.ndarray,
    size: int,
):
    # Draw the secret
    s = [random_ternary_poly(ciphertext_modulus, poly_modulus) for _ in range(size)]
    # Generate noise e
    e = [random_normal_poly(ciphertext_modulus, poly_modulus) for _ in range(size)]
    # Generate uniform poly a
    a = random_uniform_poly(ciphertext_modulus, poly_modulus)
    # RLWE instance b
    b = [a * s[i] + e[i] for i in range(size)]
    return s, a, b

def IPFE_encrypt(
    msg: list,
    pk0: QuotientRingPoly,
    pk1: list,
    ciphertext_modulus: int,
    poly_modulus: np.ndarray,
    size: int,
    delta: int,
):
    # Generate random polynomial u
    u = random_ternary_poly(ciphertext_modulus, poly_modulus)
    # Generate noises e0, e1
    e0 = random_normal_poly(ciphertext_modulus, poly_modulus)
    e1 = [random_normal_poly(ciphertext_modulus, poly_modulus) for _ in range(size)]
    # Mask the randomness u
    c0 = pk0 * u + e0
    # Mask the message with a rlwe instance (b * u + e + m)
    unit = unit_poly(ciphertext_modulus, poly_modulus)
    c1 = [pk1[i] * u + e1[i] + cst_poly(delta * msg[i], unit, ciphertext_modulus, poly_modulus) for i in range(size)]
    return c0, c1

def IPFE_keyder(
    sk: list,
    y: list,
    ciphertext_modulus: int,
    poly_modulus: np.ndarray,
    size: int,
):
    dk = cst_poly(y[0], sk[0], ciphertext_modulus, poly_modulus)
    for i in range(1, size):
        dk += cst_poly(y[i], sk[i], ciphertext_modulus, poly_modulus)
    return dk

def IPFE_decrypt(
    dk: QuotientRingPoly,
    y: list,
    c0: QuotientRingPoly,
    c1: list,
    ciphertext_modulus: int,
    poly_modulus: np.ndarray,
    delta: int,
    size: int,
):
    msg = cst_poly(y[0], c1[0], ciphertext_modulus, poly_modulus)
    for i in range(1, size):
        msg += cst_poly(y[i], c1[i], ciphertext_modulus, poly_modulus)
    msg = round((msg - c0 * dk).coef[0] / delta)
    return msg

#===================================================================
# PARAMETERS =======================================================
#===================================================================

# Ciphertext space: Z_q[X]/X^n+1
# ciphertext_modulus = q
ciphertext_modulus = 2 ** 521 - 1
# Size of the polynomials
n = 2 ** 4
# poly_modulus = X^n + 1
poly_modulus = init_poly_modulus(n)

# Size (number of entries) of the plaintext vector
size = 50

# Bounds on the entries of vectors x and y
Bx = 2**31-1
By = 100

# Scaling factor
K = size * Bx * By + 1
delta = ciphertext_modulus // K

#=======================================
# IPFE experiments =====================
#=======================================

"""
# Generate FE keys
msk, mpk0, mpk1 = IPFE_setup(ciphertext_modulus, poly_modulus, size)

# Generate random small message and encrypt it
msg = random_uniform_vector(size, Bx)
c0, c1 = IPFE_encrypt(msg, mpk0, mpk1, ciphertext_modulus, poly_modulus, size, delta)

# Generate random weight vector and encrypt it for function-hiding
y = random_uniform_vector(size, By)

# Derivation of the key
dk = IPFE_keyder(msk, y, ciphertext_modulus, poly_modulus, size)

# Retrieve the result
res = IPFE_decrypt(dk, y, c0, c1, ciphertext_modulus, poly_modulus, delta, size)

# Experiments
print("")
print("Expected result:" + str(sum(msg[i] * y[i] for i in range(size))))
print("")
print("Decryption result:" + str(res))
print("")
"""
