from utils import *
from BGV import *
from IPFE import *

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

# Plaintext space of IPFE
IPFE_modulus = 2**31 - 1

# Size (number of entries) of the plaintext vector
size = 50

# Bounds on the entries of vectors x and y
Bx = IPFE_modulus - 1
By = 100

# Scaling factor
K = size * Bx * By + 1
delta = ciphertext_modulus // K

# BGV base
base = 5

#===================================================================
# EVAL KEYDER ======================================================
#===================================================================

# Key derivation performed homomorphically
def Bifrost(
    C_msk: list, # Homomorphically encrypted master secret key
    C_y: list, # Homomorphically encrypted function (f = inner product for vector y)
    evk: QuotientRingPoly, # Evaluation key
    ciphertext_modulus: int,
    poly_modulus: np.ndarray,
    base: int,
    size: int,
):
    # Same as the IPFE keyder algorithm but with homomorphic evaluation
    dk0, dk1 = mul((C_msk[0])[0], (C_msk[0])[1], (C_y[0])[0], (C_y[0])[1], evk, base, ciphertext_modulus, poly_modulus)
    for i in range(1, size):
        c_hat0, c_hat1 = mul((C_msk[i])[0], (C_msk[i])[1], (C_y[i])[0], (C_y[i])[1], evk, base, ciphertext_modulus, poly_modulus)
        dk0, dk1 = add(c_hat0, c_hat1, dk0, dk1)
    return dk0, dk1

#=======================================
# IPFE experiments =====================
#=======================================

# Generate FE and HE keys
msk, mpk0, mpk1 = IPFE_setup(ciphertext_modulus, poly_modulus, size)
sk, pk0, pk1, evk = BGV_keygen(ciphertext_modulus, poly_modulus, IPFE_modulus, base)

# Encrypt the master secret key
C_msk = [BGV_encrypt(msk[i], pk0, pk1, ciphertext_modulus, poly_modulus, IPFE_modulus) for i in range(size)]

# Generate random small message and encrypt it under IPFE
msg = random_uniform_vector(size, Bx)
c0, c1 = IPFE_encrypt(msg, mpk0, mpk1, ciphertext_modulus, poly_modulus, size, delta)

# Generate random weight vector and encrypt it for function-hiding
y = random_uniform_vector(size, By)
# First represent it as a plaintext, then encrypt it
Pt_y = msg2pt(y, ciphertext_modulus, poly_modulus)
C_y = [BGV_encrypt((Pt_y[i])[0], pk0, pk1, ciphertext_modulus, poly_modulus, IPFE_modulus) for i in range(size)]

# Evaluation of the key derivation in the function-hiding setup
dk0, dk1 = Bifrost(C_msk, C_y, evk, ciphertext_modulus, poly_modulus, base, size)
dk_ev_fh = mod_center(BGV_decrypt(dk0, dk1, sk, IPFE_modulus), IPFE_modulus, True)

# Retrieve the result
res = IPFE_decrypt(dk_ev_fh, y, c0, c1, ciphertext_modulus, poly_modulus, delta, size)

# Experiments
print("Transciphering of FH-IPFE:" + str(res == sum(msg[i] * y[i] for i in range(size))))
print("")
