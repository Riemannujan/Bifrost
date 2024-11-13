from utils import *


def setup(
    coef_modulus: int,
    poly_modulus: np.ndarray,
    plaintext_modulus: int,
):
    # Draw the secret
    sk = random_ternary_poly(coef_modulus, poly_modulus)
    # Generate the public key
    e = random_normal_poly(coef_modulus, poly_modulus)
    a = random_uniform_poly(coef_modulus, poly_modulus)
    b = a * sk + e * plaintext_modulus
    return sk, a, b

def encrypt(
    msg: QuotientRingPoly,
    pk0: QuotientRingPoly,
    pk1: QuotientRingPoly,
    alpha: QuotientRingPoly,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    plaintext_modulus: int,
):
    u = random_ternary_poly(coef_modulus, poly_modulus)
    e0 = random_normal_poly(coef_modulus, poly_modulus)
    e1 = random_normal_poly(coef_modulus, poly_modulus)

    c0 = pk0 * alpha + e0 * plaintext_modulus + u
    c1 = pk1 * alpha + e1 * plaintext_modulus + u * msg

    return (c0, c1)

def keyder(
    sk: QuotientRingPoly,
    pk0: QuotientRingPoly,
    alpha: QuotientRingPoly,
    v: QuotientRingPoly,
):
    dk = pk0 * alpha * (sk - v)
    return dk

def decrypt(
    dk: QuotientRingPoly,
    c0: QuotientRingPoly,
    c1: QuotientRingPoly,
    v: QuotientRingPoly,
    plaintext_modulus: int,
):
    res = c1 - dk - v * c0
    return res % plaintext_modulus
