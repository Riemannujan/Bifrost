from utils import *

n = 2**4
# q
# coef_modulus = getRandomNBitInteger(32)
coef_modulus = 874
plaintext_modulus = 7
poly_modulus = init_poly_modulus(n)

def gen_secret_key(coef_modulus: int, poly_modulus: np.ndarray):
    # Draw the secret
    s = random_ternary_poly(coef_modulus, poly_modulus)
    return s

def gen_public_key(
    sk: QuotientRingPoly,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    plaintext_modulus: int,
):
    # Generate noise e
    e = random_normal_poly(coef_modulus, poly_modulus)
    # Generate unifrom poly a
    a = random_uniform_poly(coef_modulus, poly_modulus)
    # RLWE instance b
    b = a * sk + e * plaintext_modulus
    return b, -a

def bgv_encrypt(
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


def bgv_decrypt(
    c0: QuotientRingPoly,
    c1: QuotientRingPoly,
    sk: QuotientRingPoly,
    plaintext_modulus: int,
    return_noise: bool = False,
):
    msg = c0 + c1 * sk
    noise = np.max(np.abs(msg.coef))
    msg = msg % plaintext_modulus

    if return_noise:
        msg_not_reduced = np.max((c0 + c1 * sk).coef)
        # t = mod_center(polymul(c0.coef, sk.coef), c0.coef_modulus)
        # t = mod_center(polyadd(c1.coef, t), c0.coef_modulus)
        # _, t = polydiv(t, c0.poly_modulus)
        # msg_not_reduced = mod_center(t, c0.coef_modulus)
        noise = np.max(np.abs(msg_not_reduced))

        return msg, noise
    else:
        return msg


# Generate secret key
# sk = gen_secret_key(coef_modulus, poly_modulus)

# Generate public key pair
# pk0, pk1 = gen_public_key(sk, coef_modulus, poly_modulus, plaintext_modulus)

# Generate random small message in R_t
# msg = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)

# Encrypt the message into two ciphertexts
# c0, c1 = bgv_encrypt(msg, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)

# Decrypt the message and return the noise
# msg_decr, noise = bgv_decrypt(c0, c1, sk, plaintext_modulus, return_noise=True)

def add(c0, c1, c0_, c1_):
    return c0 + c0_, c1 + c1_

def mul(c0, c1, c0_, c1_):
    return c0 * c0_, c0 * c1_ + c1 * c0_, c1 * c1_

def decrypt_quad(c0, c1, c2, sk, plaintext_modulus, return_noise: bool = False):
    # Evaluate the quadratic equation
    msg = c0 + c1 * sk + c2 * sk * sk
    noise = np.max(np.abs(msg.coef))
    msg = msg % plaintext_modulus

    if return_noise:
        return msg, noise
    else:
        return msg

def pt_x_ct(
    pt: int,
    c0: QuotientRingPoly,
    c1: QuotientRingPoly,
    zero: list,
) -> QuotientRingPoly:
    if pt == 0:
        return zero
    res0, res1 = c0, c1
    for _ in range(1, pt):
        res0, res1 = add(res0, res1, c0, c1)
    return res0, res1

# sk = gen_secret_key(coef_modulus, poly_modulus)
# pk0, pk1 = gen_public_key(sk, coef_modulus, poly_modulus, plaintext_modulus)
# msg = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)
# msg_ = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)

# (c0, c1) = bgv_encrypt(msg, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)

# (c0_, c1_) = bgv_encrypt(msg_, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)
# msg_res = (msg + msg_) % plaintext_modulus
# c0_res, c1_res = add(c0, c1, c0_, c1_)
# msg_res_decr, noise = bgv_decrypt(c0_res, c1_res, sk, plaintext_modulus, return_noise=True)

# msg_res = (msg * msg_) % plaintext_modulus
# c0_res, c1_res, c2_res = mul(c0, c1, c0_, c1_)
# msg_res_decr, noise = decrypt_quad(
#     c0_res, c1_res, c2_res, sk, plaintext_modulus, return_noise=True
# )

def poly2base(poly: QuotientRingPoly, base: int) -> List[QuotientRingPoly]:
    """Converts a polynomial to a list of polynomials that represent the
    polynomial's coefficients in the given base.
    """

    coef_modulus = poly.coef_modulus
    poly_modulus = poly.poly_modulus
    n_terms = math.ceil(math.log(coef_modulus, base))
    degree = len(poly.poly_modulus) - 1
    # Make a matrix to store the coefficients
    # the rows will represent each coefficient's decomposition in `base`
    # c^(i) will be formed by taking the columns (each coefficient at step i).
    coeffs = np.zeros((degree, n_terms), dtype=object)

    for i, coef in enumerate(poly.coef):
        digits = int2base(coef % coef_modulus, base)
        # Pad with 0s
        pad_len = n_terms - len(digits)
        digits = np.pad(digits, (0, pad_len))
        coeffs[i] = digits

    # Make list of polynomials
    res = []
    for i in range(n_terms):
        p = QuotientRingPoly(
            coeffs[:, i], coef_modulus=coef_modulus, poly_modulus=poly_modulus
        )
        res.append(p)
    return res

def reconstruct_poly(poly_list: List[QuotientRingPoly], base: int) -> QuotientRingPoly:
    """Reconstructs a polynomial from a list of polynomials that represent the
    polynomial's coefficients in the given base.
    """
    if len(poly_list) == 0:
        raise ValueError("Can't construct poly from empty list")

    poly = QuotientRingPoly(
        np.array([0], dtype=object),
        poly_list[0].coef_modulus,
        poly_list[0].poly_modulus,
    )
    for i, poly_ in enumerate(poly_list):
        poly = poly + poly_ * base**i
    return poly

# t = random_uniform_poly(12345, np.array([1, 0, 0, 0, 0, 1]))
# base = 17
# res = poly2base(t, base)

# t_ = reconstruct_poly(res, base)

def gen_relinearization_key(sk, base, coef_modulus, poly_modulus, plaintext_modulus):
    n_terms = math.ceil(math.log(coef_modulus, base))

    eks = []
    for i in range(n_terms):
        # ai = random_uniform_poly(coef_modulus, poly_modulus)
        # ei = random_normal_poly(coef_modulus, poly_modulus)
        # ek0 = ai * sk + ei * plaintext_modulus + (sk * sk) * (base ** i)
        # ek1 = -ai
        b, ai = gen_public_key(sk, coef_modulus, poly_modulus, plaintext_modulus)
        ek0 = b + (sk * sk) * (base**i)
        ek1 = ai
        eks.append((ek0, ek1))
    return eks

def relinearize(c0, c1, c2, eks, base, coef_modulus, poly_modulus):
    # Decompose c2
    c2_polys = poly2base(c2, base)
    assert len(c2_polys) == len(eks)

    # Construct c0_hat, c1_hat
    c0_hat = c0
    c1_hat = c1
    for c2_i, (ek0, ek1) in zip(c2_polys, eks):
        c0_hat += c2_i * ek0
        c1_hat += c2_i * ek1

    return c0_hat, c1_hat

# max_noise = coef_modulus // 2

# base = 5
# sk = gen_secret_key(coef_modulus, poly_modulus)
# pk0, pk1 = gen_public_key(sk, coef_modulus, poly_modulus, plaintext_modulus)
# eks = gen_relinearization_key(sk, base, coef_modulus, poly_modulus, plaintext_modulus)

# msg = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)
# msg_ = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)

# (c0, c1) = bgv_encrypt(msg, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)

# (c0_, c1_) = bgv_encrypt(msg_, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)

# msg_res = (msg * msg_) % plaintext_modulus

# c0_res, c1_res, c2_res = mul(c0, c1, c0_, c1_)

# msg_quad_decr = decrypt_quad(c0_res, c1_res, c2_res, sk, plaintext_modulus)
# c0_hat, c1_hat = relinearize(
#     c0_res, c1_res, c2_res, eks, base, coef_modulus, poly_modulus
# )

# msg_res_decr, noise = bgv_decrypt(c0_hat, c1_hat, sk, plaintext_modulus, return_noise=True)
 
def scale(
    x: QuotientRingPoly, big_mod: int, small_mod: int, plaintext_modulus: int
) -> QuotientRingPoly:
    """Scales a polynomial x from a big coefficient modulus to a smaller one.
    Adjust post scaling."""

    # Scale x by small_mod / big_mod
    x_ = x * (small_mod / big_mod)

    # Adjust x_ to be closest vector such that x_ = x mod plaintext_modulus
    # We calculate the difference and add it.
    diff = (x.coef % plaintext_modulus) - (x_.coef % plaintext_modulus)
    # print(diff % plaintext_modulus)
    x_.coef += diff
    assert all(x.coef % plaintext_modulus == x_.coef % plaintext_modulus)

    return x_

def scale2(
    x: QuotientRingPoly, big_mod: int, small_mod: int, plaintext_modulus: int
) -> QuotientRingPoly:
    """Scales a polynomial x from a big coefficient modulus to a smaller one.
    Adjust pre scaling."""

    assert big_mod % small_mod == 0
    delta = big_mod // small_mod

    # Adjustment term (small delta)
    d = ((-x * pow(plaintext_modulus, -1, delta)) % delta) * plaintext_modulus

    assert d % delta == (-x % delta)
    assert all((d % plaintext_modulus).coef == 0)

    # Scale x by small_mod / big_mod
    x_ = (x + d) * (small_mod / big_mod)
    assert all(x.coef % plaintext_modulus == x_.coef % plaintext_modulus)

    return x_

# We now consider Q - big_mod, q - small_mod and q | Q
# while True:
#     small_mod = getRandomNBitInteger(32)
#     delta = getPrime(96)
#     big_mod = small_mod * delta
#     if (
#         delta % plaintext_modulus == 1
#     ):  # and small_mod % plaintext_modulus == 1 and big_mod % plaintext_modulus == 1:
#         break

# Generate keys
# sk = gen_secret_key(big_mod, poly_modulus)
# pk0, pk1 = gen_public_key(sk, big_mod, poly_modulus, plaintext_modulus)

# Generate message
# msg = random_uniform_poly(big_mod, poly_modulus, high=plaintext_modulus)

# Encrypt message
# (c0, c1) = bgv_encrypt(msg, pk0, pk1, big_mod, poly_modulus, plaintext_modulus)
# c0_scaled = scale2(c0, big_mod, small_mod, plaintext_modulus)
# c1_scaled = scale2(c1, big_mod, small_mod, plaintext_modulus)

# Change ring
# c0_scaled.coef_modulus = small_mod
# c1_scaled.coef_modulus = small_mod
# sk_ = sk.copy()
# sk_.coef_modulus = small_mod
# msg_decr_before, noise_before = bgv_decrypt(c0, c1, sk, plaintext_modulus, return_noise=True)
# msg_decr_after, noise_after = bgv_decrypt(c0_scaled, c1_scaled, sk_, plaintext_modulus, return_noise=True)

def gen_relinearization_key_ms(sk, big_mod, small_mod, poly_modulus, plaintext_modulus):
    sk = sk.copy()
    sk.coef_modulus = big_mod

    a = random_uniform_poly(coef_modulus=big_mod, poly_modulus=poly_modulus)
    e = random_normal_poly(big_mod, poly_modulus)

    ek0 = a * sk + e * plaintext_modulus + (sk * sk) * (big_mod // small_mod)
    ek1 = -a

    return ek0, ek1

def relinearize_ms(
    c0, c1, c2, ek0, ek1, big_mod, small_mod, poly_modulus, plaintext_modulus
):
    # Copy ciphertexts and change ring to big_mod.
    c0 = c0.copy()
    c0.coef_modulus = big_mod
    c1 = c1.copy()
    c1.coef_modulus = big_mod
    c2 = c2.copy()
    c2.coef_modulus = big_mod

    delta = big_mod // small_mod

    c0_hat = c0 * delta + c2 * ek0
    c1_hat = c1 * delta + c2 * ek1

    c0_hat = scale2(c0_hat, big_mod, small_mod, plaintext_modulus)
    c1_hat = scale2(c1_hat, big_mod, small_mod, plaintext_modulus)

    return c0_hat, c1_hat

# We now consider Q - big_mod, q - small_mod and q | Q
"""
while True:
    small_mod = getRandomNBitInteger(32)
    delta = getPrime(96)
    big_mod = small_mod * delta
    if (
        delta % plaintext_modulus == 1
    ):  # and small_mod % plaintext_modulus == 1 and big_mod % plaintext_modulus == 1:
        break
"""
# sk = gen_secret_key(small_mod, poly_modulus)
# pk0, pk1 = gen_public_key(sk, small_mod, poly_modulus, plaintext_modulus)

# ek0, ek1 = gen_relinearization_key_ms(sk, big_mod, small_mod, poly_modulus, plaintext_modulus)

# msg = random_uniform_poly(small_mod, poly_modulus, high=plaintext_modulus)
# msg_ = random_uniform_poly(small_mod, poly_modulus, high=plaintext_modulus)

# (c0, c1) = bgv_encrypt(msg, pk0, pk1, small_mod, poly_modulus, plaintext_modulus)

# (c0_, c1_) = bgv_encrypt(msg_, pk0, pk1, small_mod, poly_modulus, plaintext_modulus)

# msg_res = (msg * msg_) % plaintext_modulus

# c0_res, c1_res, c2_res = mul(c0, c1, c0_, c1_)

# msg_quad_decr = decrypt_quad(c0_res, c1_res, c2_res, sk, plaintext_modulus)
"""
c0_hat, c1_hat = relinearize_ms(
    c0_res,
    c1_res,
    c2_res,
    ek0,
    ek1,
    big_mod,
    small_mod,
    poly_modulus,
    plaintext_modulus,
)
"""

# msg_res_decr, noise = bgv_decrypt(c0_hat, c1_hat, sk, plaintext_modulus, return_noise=True)
 
# Change ring
# c0_hat.coef_modulus = small_mod
# c1_hat.coef_modulus = small_mod
# sk_ = sk.copy()
# sk_.coef_modulus = small_mod

# msg_res_decr, noise = bgv_decrypt(c0_hat, c1_hat, sk_, plaintext_modulus, return_noise=True)

"""
n = 2**4
plaintext_modulus = 7
poly_modulus = init_poly_modulus(n)
coef_modulus = getRandomNBitInteger(32)
base = 5

msg0 = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)
msg1 = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)
msg2 = random_uniform_poly(coef_modulus, poly_modulus, high=plaintext_modulus)

msg_res = (msg0 * msg1 + msg2) % plaintext_modulus

sk = gen_secret_key(coef_modulus, poly_modulus)
pk0, pk1 = gen_public_key(sk, coef_modulus, poly_modulus, plaintext_modulus)
eks = gen_relinearization_key(sk, base, coef_modulus, poly_modulus, plaintext_modulus)
c00, c01 = bgv_encrypt(msg0, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)
c10, c11 = bgv_encrypt(msg1, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)
c20, c21 = bgv_encrypt(msg2, pk0, pk1, coef_modulus, poly_modulus, plaintext_modulus)
c_star0, c_star1, c_star2 = mul(c00, c01, c10, c11)
c_hat0, c_hat1 = relinearize(c_star0, c_star1, c_star2, eks, base, coef_modulus, poly_modulus)
c_add0, c_add1 = add(c_hat0, c_hat1, c20, c21)
msg_decr = bgv_decrypt(c_add0, c_add1, sk, plaintext_modulus)

print(msg_res)
print(msg_decr)
"""
