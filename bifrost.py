from utils import *
from BGV import *
from IPFE import *
from SPADE import *


#===================================================================
# PARAMETERS =======================================================
#===================================================================

# Ciphertext space
# q
coef_modulus = 2 ** 31 - 1
# Size of the polynomials
n = 2 ** 4
poly_modulus = init_poly_modulus(n)

# BGV base
base = 5

# Plaintext space of SPADE
SPD_modulus = coef_modulus

# Plaintext space of IPFE
IPFE_modulus = 97

# Size of the plaintext vector
size = 5

# Bounds on x and y
Bx = 10
By = 10
K = size * Bx * By + 1

# Scaling factor
delta = coef_modulus // K


#===================================================================
# EVAL KEYDER ======================================================
#===================================================================


# Key derivation with a, y as a ciphertext for function hiding
def SPD_ev_fh(
    c_msk0: QuotientRingPoly,
    c_msk1: QuotientRingPoly,
    c_v0: QuotientRingPoly,
    c_v1: QuotientRingPoly,
    c_mpk00: QuotientRingPoly,
    c_mpk01: QuotientRingPoly,
    c_alpha0: QuotientRingPoly,
    c_alpha1: QuotientRingPoly,
    coef_modulus: int,
    poly_modulus: int,
    base: int,
):
    # Computes HE evaluation of msk - v0 encryption
    res0, res1 = add(c_msk0, c_msk1, c_v0, c_v1)
    # Computes HE evaluation of alpha * (msk - v0) encryption
    res0, res1 = mul(res0, res1, c_alpha0, c_alpha1, evk, base, coef_modulus, poly_modulus)
    # Computes HE evaluation of mpk0 * alpha * (msk - v0) encryption
    res0, res1 = mul(res0, res1, c_mpk00, c_mpk01, evk, base, coef_modulus, poly_modulus)
    return res0, res1


# Key derivation with y as a plaintext
def IPFE_ev(
    C_msk: list,
    y: list,
    zero: list,
    size: int,
):
    dk0, dk1 = cst_x_ct(y[0], (C_msk[0])[0], (C_msk[0])[1], zero)
    for i in range(1, size):
        c0, c1 = cst_x_ct(y[i], (C_msk[i])[0], (C_msk[i])[1], zero)
        dk0, dk1 = add(dk0, dk1, c0, c1)
    return dk0, dk1


# Key derivation with y as a ciphertext for function hiding
def IPFE_ev_fh(
    C_msk: list,
    C_y: list,
    evk: QuotientRingPoly,
    coef_modulus: int,
    poly_modulus: np.ndarray,
    base: int,
    size: int,
):
    dk0, dk1 = mul((C_msk[0])[0], (C_msk[0])[1], (C_y[0])[0], (C_y[0])[1], evk, base, coef_modulus, poly_modulus)
    for i in range(1, size):
        c_hat0, c_hat1 = mul((C_msk[i])[0], (C_msk[i])[1], (C_y[i])[0], (C_y[i])[1], evk, base, coef_modulus, poly_modulus)
        dk0, dk1 = add(c_hat0, c_hat1, dk0, dk1)
    return dk0, dk1


#==============================================================================================
# SPADE EXPERIMENT ============================================================================
#==============================================================================================

# Generate FE and HE keys
msk, mpk0, mpk1 = SPD_setup(coef_modulus, poly_modulus, SPD_modulus)
sk, pk0, pk1, evk = BGV_keygen(coef_modulus, poly_modulus, SPD_modulus, base)

# Encrypt the master secret key and public key
c_msk0, c_msk1 = BGV_encrypt(msk, pk0, pk1, coef_modulus, poly_modulus, SPD_modulus)
c_mpk00, c_mpk01 = BGV_encrypt(mpk0, pk0, pk1, coef_modulus, poly_modulus, SPD_modulus)

# Generate user's identifiers
alpha = random_ternary_poly(coef_modulus, poly_modulus)
c_alpha0, c_alpha1 = BGV_encrypt(alpha, pk0, pk1, coef_modulus, poly_modulus, SPD_modulus)

# Generate random small message and encrypt it
msg = random_uniform_poly(coef_modulus, poly_modulus, high=SPD_modulus)
c0, c1 = SPD_encrypt(msg, mpk0, mpk1, alpha, coef_modulus, poly_modulus, SPD_modulus)

# Generate the seeked value v
v = random_uniform_poly(coef_modulus, poly_modulus, high=SPD_modulus)
v = msg
# Encrypt -v
c_v0, c_v1 = BGV_encrypt(-v, pk0, pk1, coef_modulus, poly_modulus, SPD_modulus)

# Evaluation of the key derivation
dk0, dk1 = SPD_ev_fh(
    c_msk0,
    c_msk1,
    c_v0,
    c_v1,
    c_mpk00,
    c_mpk01,
    c_alpha0,
    c_alpha1,
    coef_modulus,
    poly_modulus,
    base,
)

# Decryption of the decryption key
dk_ev = BGV_decrypt(dk0, dk1, sk, SPD_modulus)
res_ev = SPD_decrypt(dk_ev, c0, c1, v, SPD_modulus)

# Traditional SPADE
dk = SPD_keyder(msk, mpk0, alpha, v)
res = SPD_decrypt(dk, c0, c1, v, SPD_modulus)

# Experiments
print("Transciphering of SPADE:" + str(res_ev == res))
print("")


#=======================================
# IPFE experiments =====================
#=======================================

# Generate FE and HE keys
msk, mpk0, mpk1 = IPFE_setup(coef_modulus, poly_modulus, size)
sk, pk0, pk1, evk = BGV_keygen(coef_modulus, poly_modulus, IPFE_modulus, base)

# Encrypt the master secret key
C_msk = [BGV_encrypt(msk[i], pk0, pk1, coef_modulus, poly_modulus, IPFE_modulus) for i in range(size)]

# Generate random small message and encrypt it
msg = random_uniform_vector(size, Bx)
c0, c1 = IPFE_encrypt(msg, mpk0, mpk1, coef_modulus, poly_modulus, size, delta)

# Generate random weight vector and encrypt it for function-hiding
y = random_uniform_vector(size, By)

# Encrypt the function
Pt_y = msg2pt(y, coef_modulus, poly_modulus)
C_y = [BGV_encrypt((Pt_y[i])[0], pk0, pk1, coef_modulus, poly_modulus, IPFE_modulus) for i in range(size)]

# Evaluation of the key derivation in the function-hiding setup
dk0, dk1 = IPFE_ev_fh(C_msk, C_y, evk, coef_modulus, poly_modulus, base, size)
dk_ev_fh = mod_center(BGV_decrypt(dk0, dk1, sk, IPFE_modulus), IPFE_modulus, True)

# Evaluation of the key derivation in the normal setup
zero = BGV_encrypt(QuotientRingPoly([0], coef_modulus, poly_modulus), pk0, pk1, coef_modulus, poly_modulus, IPFE_modulus)
dk0, dk1 = IPFE_ev(C_msk, y, zero, size)
dk_ev = mod_center(BGV_decrypt(dk0, dk1, sk, IPFE_modulus), IPFE_modulus, True)

# Retrieve the result
res_ev = IPFE_decrypt(dk_ev, y, c0, c1, coef_modulus, poly_modulus, delta, size)
res_ev_fh = IPFE_decrypt(dk_ev_fh, y, c0, c1, coef_modulus, poly_modulus, delta, size)
res = sum(msg[i] * y[i] for i in range(size))

# Experiments
print("Transciphering of IPFE:" + str(res_ev == res))
print("")
print("Transciphering of FH-IPFE:" + str(res_ev_fh == res))
print("")
