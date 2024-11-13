from utils import *


def IPFE_setup(
    coef_modulus: int,
    poly_modulus: np.ndarray,
    size: int,
):
    # Draw the secret
    s = [random_ternary_poly(coef_modulus, poly_modulus) for _ in range(size)]
    # Generate noise e
    e = [random_normal_poly(coef_modulus, poly_modulus) for _ in range(size)]
    # Generate uniform poly a
    a = random_uniform_poly(coef_modulus, poly_modulus)
    # RLWE instance b
    b = [a * s[i] + e[i] for i in range(size)]
    return s, a, b

def IPFE_encrypt(
    msg: list,
    pk0: QuotientRingPoly,
    pk1: list,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    size: int,
    delta: int,
):
    # Generate random polynomial u
    u = random_ternary_poly(coef_modulus, poly_modulus)
    # Generate noises e0, e1
    e0 = random_normal_poly(coef_modulus, poly_modulus)
    e1 = [random_normal_poly(coef_modulus, poly_modulus) for _ in range(size)]
    # Mask the randomness u
    c0 = pk0 * u + e0
    # Mask the message with a rlwe instance (b * u + e + m)
    unit = unit_poly(coef_modulus, poly_modulus)
    c1 = [pk1[i] * u + e1[i] + cst_poly(delta * msg[i], unit, coef_modulus, poly_modulus) for i in range(size)]
    return c0, c1

def IPFE_keyder(
    sk: list,
    y: list,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    size: int,
):
    dk = cst_poly(y[0], sk[0], coef_modulus, poly_modulus)
    for i in range(1, size):
        dk += cst_poly(y[i], sk[i], coef_modulus, poly_modulus)
    return dk

def IPFE_decrypt(
    dk: QuotientRingPoly,
    y: list,
    c0: QuotientRingPoly,
    c1: list,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    delta: int,
    size: int,
):
    msg = cst_poly(y[0], c1[0], coef_modulus, poly_modulus)
    for i in range(1, size):
        msg += cst_poly(y[i], c1[i], coef_modulus, poly_modulus)
    msg = round((msg - c0 * dk).coef[0] / delta)
    return msg

