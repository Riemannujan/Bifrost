from utils import *

#==============================================================#
# NOTATIONS ===================================================#
#==============================================================#

# coef_modulus (q): modulus of the finite group Z_q where the ciphertext lies
# poly_modulus: the polynomial X^n + 1 that divides the polynomial field
# plaintext_modulus (t): maximum size of the plaintext entries, modulus of Z_t 
# msg lies in Z_t[X] / X^n + 1
# sk, pk0, pk1 lie in Z_q

#==============================================================#
# BGV PKE scheme ==============================================#
#==============================================================#

def BGV_keygen(
    coef_modulus: int,
    poly_modulus: np.ndarray,
    plaintext_modulus,
    base: int,
):
    # Draw the secret key
    sk = random_ternary_poly(coef_modulus, poly_modulus)
    # Generate the public key
    a = random_uniform_poly(coef_modulus, poly_modulus)
    e = random_normal_poly(coef_modulus, poly_modulus)
    b = a * sk + e * plaintext_modulus

    # Generate the evaluation (relinearization) key
    n_terms = math.ceil(math.log(coef_modulus, base))
    evk = []
    for i in range(n_terms):
        ai = random_uniform_poly(coef_modulus, poly_modulus)
        ei = random_normal_poly(coef_modulus, poly_modulus)
        evk0 = ai * sk + ei * plaintext_modulus + (sk * sk) * (base ** i)
        evk1 = -ai
        evk.append((evk0, evk1))

    return sk, b, -a, evk

def BGV_encrypt(
    msg: QuotientRingPoly,
    pk0: QuotientRingPoly,
    pk1: QuotientRingPoly,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    plaintext_modulus: int,
):
    u = random_ternary_poly(coef_modulus, poly_modulus)
    e0 = random_normal_poly(coef_modulus, poly_modulus)
    e1 = random_normal_poly(coef_modulus, poly_modulus)
    # Mask the message with a rlwe instance (b * r + te)
    c0 = pk0 * u + e0 * plaintext_modulus + msg
    c1 = pk1 * u + e1 * plaintext_modulus
    return c0, c1

def BGV_decrypt(
    c0: QuotientRingPoly,
    c1: QuotientRingPoly,
    sk: QuotientRingPoly,
    plaintext_modulus: int,
):
    msg = c0 + c1 * sk
    return msg % plaintext_modulus


#==============================================================#
# EVALUATION ==================================================#
#==============================================================#


# Performs the homomorphic addition of 2 ciphertexts
def add(
    c0: QuotientRingPoly,
    c1: QuotientRingPoly,
    d0: QuotientRingPoly,
    d1: QuotientRingPoly
):
    return c0 + d0, c1 + d1


# Multiplies a ciphertext by a constant
def cst_x_ct(
    pt: int,
    c0: QuotientRingPoly,
    c1: QuotientRingPoly,
    zero: list,
):
    if pt == 0:
        return zero
    res0, res1 = c0, c1
    for _ in range(1, pt):
        res0, res1 = add(res0, res1, c0, c1)
    return res0, res1


def mul(c00, c01, c10, c11, evk, base, coef_modulus, poly_modulus):

    # Preliminary multiplication
    c0 = c00 * c10
    c1 = c00 * c11 + c01 * c10
    c2 = c01 * c11

    # Decompose c2
    n_terms = math.ceil(math.log(coef_modulus, base))
    degree = len(poly_modulus) - 1
    coeffs = np.zeros((degree, n_terms), dtype=object)

    for i, coef in enumerate(c2.coef):
        digits = int2base(coef % coef_modulus, base)
        # Pad with 0s
        pad_len = n_terms - len(digits)
        digits = np.pad(digits, (0, pad_len))
        coeffs[i] = digits

    # Make list of polynomials
    c2_polys = []
    for i in range(n_terms):
        p = QuotientRingPoly(
            coeffs[:, i], coef_modulus=coef_modulus, poly_modulus=poly_modulus
        )
        c2_polys.append(p)

    assert len(c2_polys) == len(evk)

    # Construct c0_hat, c1_hat
    c0_hat = c0
    c1_hat = c1
    for c2_i, (ek0, ek1) in zip(c2_polys, evk):
        c0_hat += c2_i * ek0
        c1_hat += c2_i * ek1

    return c0_hat, c1_hat
