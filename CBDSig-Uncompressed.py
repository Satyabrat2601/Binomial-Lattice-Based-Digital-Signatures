import os
import hashlib
from typing import List, Tuple, Optional

# ----------------------------
# Parameters
# ----------------------------
Q = 7681    # Prime modulus
N = 256      # Polynomial degree

ML_DSA_PARAMS = {
    44: {'k': 2, 'l': 3, 'eta': 1, 'gamma1': 339, 'tau': 39, 'beta': 39},
    65: {'k': 3, 'l': 4, 'eta': 1, 'gamma1': 469, 'tau': 49, 'beta': 49},
    87: {'k': 4, 'l': 6, 'eta': 1, 'gamma1': 510, 'tau': 60, 'beta': 60},
}

# ----------------------------
# Basic modular helpers
# ----------------------------
def mod_exp(base, exp, mod):
    r = 1
    base %= mod
    while exp > 0:
        if exp & 1:
            r = (r * base) % mod
        base = (base * base) % mod
        exp >>= 1
    return r

def sym_mod(a, mod):
    """Symmetric modulus"""
    r = a % mod
    if mod % 2 == 0:
        if r > mod/2:
            r = r - mod
        if r <= - mod/2:
            r = r + mod
    else:
        if r > (mod-1)/2:
            r = r - mod
        if r < - mod/2:
            r = r + mod
    return r

# ----------------------------
# NTT (proper twiddle scheduling, no bit-reversal)
# ----------------------------
def _powmod(a, e, m):
    r = 1
    a %= m
    while e:
        if e & 1:
            r = (r * a) % m
        a = (a * a) % m
        e >>= 1
    return r

def factorize(n):
    factors = set()

    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1

    if n > 1:
        factors.add(n)

    return factors


def find_primitive_root(q):
    """
    Find primitive generator of F_q*
    """

    phi = q - 1
    factors = factorize(phi)

    for g in range(2, q):

        ok = True

        for p in factors:
            if _powmod(g, phi // p, q) == 1:
                ok = False
                break

        if ok:
            return g

    raise ValueError("No primitive root found")


def find_principal_nth_root(q, n):
    """
    Find primitive n-th root of unity modulo q.

    Requires:
        n | (q-1)
    """

    assert (q - 1) % n == 0

    g = find_primitive_root(q)

    omega = _powmod(g, (q - 1) // n, q)

    # Verify exact order n
    assert _powmod(omega, n, q) == 1

    tmp = n

    for p in factorize(n):
        assert _powmod(omega, n // p, q) != 1

    return omega

PSI = find_principal_nth_root(Q, N)

INV_PSI = _powmod(PSI, Q - 2, Q)

N_INV = _powmod(N, Q - 2, Q)

def ntt(a: List[int]) -> List[int]:
    """Iterative Cooley–Tukey DIT NTT (length N)."""
    if len(a) != N:
        raise ValueError(f"Polynomial must have exactly {N} coefficients")
    A = [x % Q for x in a]
    m = 1
    while m < N:
        w_m = _powmod(PSI, N // (2*m), Q)
        for k in range(0, N, 2*m):
            w = 1
            for j in range(m):
                t = (w * A[k + j + m]) % Q
                u = A[k + j]
                A[k + j]     = (u + t) % Q
                A[k + j + m] = (u - t) % Q
                w = (w * w_m) % Q
        m <<= 1
    return A

def ntt_inv(A: List[int]) -> List[int]:
    """Inverse of the above NTT (Gentleman–Sande style) with 1/N scaling."""
    if len(A) != N:
        raise ValueError(f"Polynomial must have exactly {N} coefficients")
    a = [x % Q for x in A]
    m = N >> 1
    while m >= 1:
        w_m = _powmod(INV_PSI, N // (2*m), Q)
        for k in range(0, N, 2*m):
            w = 1
            for j in range(m):
                u = a[k + j]
                v = a[k + j + m]
                a[k + j]     = (u + v) % Q
                a[k + j + m] = ((u - v) * w) % Q
                w = (w * w_m) % Q
        m >>= 1
    for i in range(N):
        a[i] = (a[i] * N_INV) % Q
    return a

# ----------------------------
# Poly ops
# ----------------------------
def poly_add(a: List[int], b: List[int], mod: int = Q) -> List[int]:
    if len(a) != N or len(b) != N:
        raise ValueError("Polynomials must have exactly N coefficients")
    return [(a[i] + b[i]) % mod for i in range(N)]

def poly_sub(a: List[int], b: List[int], mod: int = Q) -> List[int]:
    if len(a) != N or len(b) != N:
        raise ValueError("Polynomials must have exactly N coefficients")
    return [(a[i] - b[i]) % mod for i in range(N)]

def poly_scalar_mult(poly: List[int], scalar: int, mod: int = Q) -> List[int]:
    return [(coeff * scalar) % mod for coeff in poly]

def poly_mult_ntt(a: List[int], b: List[int]) -> List[int]:
    """Polynomial multiplication modulo Q using NTT."""
    if len(a) != N or len(b) != N:
        raise ValueError(f"Both polynomials must have exactly {N} coefficients")
    a_ntt = ntt(a)
    b_ntt = ntt(b)
    c_ntt = [(a_ntt[i] * b_ntt[i]) % Q for i in range(N)]
    return ntt_inv(c_ntt)

# ----------------------------
# SHAKE helpers
# ----------------------------
def shake256_absorb_finalize(data: bytes):
    shake = hashlib.shake_256()
    shake.update(data)
    return shake

def shake256_squeeze(shake_ctx, length: int) -> bytes:
    return shake_ctx.digest(length)

# ----------------------------
# Expand A[i,j] via rejection sampling (q < 2^16)
# ----------------------------
def expand_a_element(rho: bytes, i: int, j: int) -> List[int]:
    """
    Generate A[i,j] coefficients uniformly in [0, Q) using SHAKE128.
    rho must be 32 bytes. Input domain separation: rho || j || i (FIPS 204 style).
    """
    if len(rho) != 32:
        raise ValueError("Seed rho must be exactly 32 bytes")

    xof = hashlib.shake_128(rho + bytes([j, i]))
    # Expect ~ N * (65536/Q) candidates; take a generous buffer
    buf = xof.digest(2 * N * 3)
    coeffs, idx = [], 0
    while len(coeffs) < N:
        if idx + 2 > len(buf):
            more = xof.digest(len(buf) + 1024)[len(buf):]
            buf += more
        val = buf[idx] | (buf[idx+1] << 8)
        idx += 2
        if val < Q:
            coeffs.append(val)
    return coeffs

# ----------------------------
# CBD sampler (general eta) – returns coefficients in [0,Q)
# ----------------------------
def sample_poly_cbd(B: bytes, eta: int) -> List[int]:
    """Center Binomial Distribution; keeps coefficients modulo Q in [0,Q)."""
    if len(B) * 8 < 2 * eta * N:
        raise ValueError(f"Byte array too short for eta={eta} and N={N}. "
                         f"Need {2*eta*N} bits, got {len(B)*8}")
    f = [0] * N
    bit_index = 0
    for i in range(N):
        a = 0
        b = 0
        for _ in range(eta):
            byte_pos = (bit_index // 8)
            bit_offset = bit_index % 8
            a += (B[byte_pos] >> bit_offset) & 1
            bit_index += 1
        for _ in range(eta):
            byte_pos = (bit_index // 8)
            bit_offset = bit_index % 8
            b += (B[byte_pos] >> bit_offset) & 1
            bit_index += 1
        f[i] = sym_mod((a - b), Q)
    return f

# ----------------------------
# Matrix-vector (over R_q[x]/(x^N+1)) using NTT
# ----------------------------
def matrix_vector_mult_ntt(A: List[List[List[int]]], s: List[List[int]]) -> List[List[int]]:
    """
    Compute (k x l) polynomial matrix times (l) polynomial vector using NTT.
    All arithmetic modulo Q.
    """
    k = len(A)
    l = len(A[0]) if A else 0
    if len(s) != l:
        raise ValueError(f"Vector s must have length {l}")
    result = []
    for i in range(k):
        acc = [0] * N
        for j in range(l):
            prod = poly_mult_ntt(A[i][j], s[j])
            acc = poly_add(acc, prod, mod=Q)
        result.append(acc)
    return result

# ----------------------------
# Demo / Test pipeline
# ----------------------------
if __name__ == "__main__":
    k = ML_DSA_PARAMS[65]['k']
    l = ML_DSA_PARAMS[65]['l']
    eta = ML_DSA_PARAMS[65]['eta']

    print(f"== Modified ML-DSA-44 keygen test ==")
    print(f"Parameters: Q={Q}, N={N}, k={k}, l={l}, eta={eta}")
    print(f"NTT sanity: PSI^N mod Q = {_powmod(PSI, N, Q)} (should be 1)")

    # Seeds
    seed = b'\x42' * 32
    zeta = seed[:32] if seed is not None else os.urandom(32)

    shake_ctx = shake256_absorb_finalize(zeta + b'\x00')
    rho = shake256_squeeze(shake_ctx, 32)

    shake_ctx = shake256_absorb_finalize(zeta + b'\x01')
    rho_prime = shake256_squeeze(shake_ctx, 64)

    print("\n[Step 1] Generate A1 in R_q^{k x l} ...")
    A1: List[List[List[int]]] = []
    for i in range(k):
        row = []
        for j in range(l):
            row.append(expand_a_element(rho, i, j))
        A1.append(row)
    print("A1 generated.",A1)

    print("\n[Step 2] Sample secrets s1 (length l) and s2 (length k) via CBD_eta ...")
    s1: List[List[int]] = []
    s2: List[List[int]] = []

    for i in range(l):
        nonce = i.to_bytes(2, 'little')
        xof = shake256_absorb_finalize(rho_prime + nonce)
        B = shake256_squeeze(xof, 64 * eta)
        s1.append(sample_poly_cbd(B, eta))
    
    print("s1 sampled", s1)

    for i in range(k):
        nonce = (l + i).to_bytes(2, 'little')
        xof = shake256_absorb_finalize(rho_prime + nonce)
        B = shake256_squeeze(xof, 64 * eta)
        s2.append(sample_poly_cbd(B, eta))

    print("s2 sampled.", s2)

    print("\n[Step 3] Compute b = A1*s1 + s2 over R_q ...")
    A1_s1 = matrix_vector_mult_ntt(A1, s1)       # k polys
    b: List[List[int]] = []
    for i in range(k):
        b.append(poly_add(A1_s1[i], s2[i], mod=Q))
    print("b computed.", b)

    print("\n[Step 4] Build column (-2b + q*e0) (integer domain) ...")
    neg_2b_plus_qe0: List[List[int]] = []
    for i in range(k):
        poly = [(-2 * b[i][t]) for t in range(N)]   # -2b[i] as polynomial
        if i == 0: 
            poly[0] += Q                               # add q to constant coefficient (q*e0)
        neg_2b_plus_qe0.append(poly)
    print("(-2b + q*e0) column ready.", neg_2b_plus_qe0)

    # Instead of forming the huge block-matrix A and multiplying it by s,
    # compute A*s mod 2Q directly:
    # A*s = (-2b + q*e0) + 2*(A1*s1) + 2*s2  (mod 2Q)
    print("\n[Step 5] Compute A·s mod 2Q (without forming A explicitly) ...")
    twoQ = 2 * Q
    result_mod_2Q: List[List[int]] = []
    for i in range(k):
        term = [0] * N
        for t in range(N):
            val = neg_2b_plus_qe0[i][t] + 2 * A1_s1[i][t] + 2 * s2[i][t]
            term[t] = val % twoQ
        result_mod_2Q.append(term)
    print("A·s mod 2Q computed.", result_mod_2Q)

    # Expected: q*e0 (i.e., constant coeff == Q, others == 0) for each row
    expected = [[0]*N for _ in range(k)]
    for i in range(k):
        expected[0][0] = Q

    ok = (result_mod_2Q == expected)
    print("\n[Step 6] Check: A·s ≡ q·e0 (mod 2Q) ? ->", "PASS" if ok else "FAIL")

    # Optional: show a small slice to avoid huge prints
    print("\nResult[0][:8] =", result_mod_2Q[0][:8], "(first 8 coeffs)")
    print("Expected   [:8] =", expected[0][:8], "(first 8 coeffs)")

    # Strong assertion
    assert ok, "A·s mod 2Q does not equal q·e0 — something is off."
    print("\nAll good.")

def build_big_matrix(b, A1, Q):
    k = len(b)       # number of rows
    l = len(A1[0])   # number of columns in A1
    N = len(b[0])    # polynomial length
    
    # Step 1: -2b + q*e0
    neg2b_qe0 = []
    for i in range(k):
        poly = [(-2 * b[i][t]) for t in range(N)]
        if i == 0:
            poly[0] += Q   # only first row gets +Q
        neg2b_qe0.append(poly)
    
    # Step 2: 2*A1
    twoA1 = [[[2 * coeff for coeff in A1[i][j]] for j in range(l)] for i in range(k)]
    
    # Step 3: 2*I_k (as polynomials)
    twoIk = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == j:
                row.append([2] + [0]*(N-1))  # polynomial "2"
            else:
                row.append([0]*N)            # zero polynomial
        twoIk.append(row)
    
    # Step 4: Concatenate blocks horizontally
    big_matrix = []
    for i in range(k):
        row = [neg2b_qe0[i]] + twoA1[i] + twoIk[i]
        big_matrix.append(row)
    
    return big_matrix

A = build_big_matrix(b, A1, Q)
print("Shape:", len(A), "x", len(A[0]))
print("First row of A:", A[0])

s = [[1] + [0] * (N - 1)] + s1 + s2
#print(s)
############################################################################################################################
####################################                  SIGNING                 ##############################################
############################################################################################################################

def hash_to_ball(message: bytes, w1_bytes: bytes, tau: int) -> List[int]:
    """IMPROVED: Hash function H(M || w₁) ∈ B_τ using FIPS 204 SHAKE256"""
    # Use SHAKE256 as per FIPS 204
    shake_input = message + w1_bytes
    shake_ctx = shake256_absorb_finalize(shake_input)
    hash_bytes = shake256_squeeze(shake_ctx, 256)  # Get more bytes for better distribution
    
    # Convert to polynomial in B_τ (τ coefficients are ±1, rest are 0)
    c = [0] * N
    
    # Use rejection sampling to get τ distinct positions
    positions = set()
    byte_idx = 0
    
    # Get τ distinct positions
    while len(positions) < tau and byte_idx < len(hash_bytes) - 1:
        pos = (hash_bytes[byte_idx] | (hash_bytes[byte_idx + 1] << 8)) % N
        positions.add(pos)
        byte_idx += 2
        
        # Get more bytes if needed
        if byte_idx >= len(hash_bytes) - 1:
            hash_bytes += shake256_squeeze(shake_ctx, 256)
    
    # Set signs for the chosen positions
    positions_list = list(positions)[:tau]
    sign_bytes = shake256_squeeze(shake_ctx, tau)
    
    for i, pos in enumerate(positions_list):
        # Set coefficient to ±1 based on sign bit
        sign = 1 if (sign_bytes[i] & 1) == 0 else -1
        c[pos] = sign
    
    return c


# Signing loop
z = None
attempt = 0
max_attempts = 5  # Reduced for faster testing

print("Starting signing loop...")
message = b"Hello, Modified ML-DSA with Verification!"
while z is None and attempt < max_attempts:
    attempt += 1
    print(f"  Attempt {attempt}/{max_attempts}", end="", flush=True)

    y = []
    cbd_eta = ML_DSA_PARAMS[65]['gamma1']
    bytes_per_poly = (2 * cbd_eta * N + 7) // 8

    # we'll derive unique nonces for each poly using a master XOF rather than one huge slice
    master_nonce = (1000 + attempt).to_bytes(4, 'little')
    shake_input = rho_prime + master_nonce
    shake_ctx = shake256_absorb_finalize(shake_input)

    # generate y polynomials one-by-one (l + k total)
    total_polys = l + k + 1  # if you actually need l+k+1, change accordingly
    for i in range(total_polys):
        # domain-separate per-poly using an index-based nonce appended to XOF input
        per_poly_xof = hashlib.shake_128()
        per_poly_xof.update(shake_input + i.to_bytes(2, 'little'))
        B_slice = per_poly_xof.digest(bytes_per_poly)
        poly = sample_poly_cbd(B_slice, cbd_eta)
        y.append(poly)
    print("The vector y:", y)

    # compute w = A * y
    w = matrix_vector_mult_ntt(A, y)

    import struct

    def poly_to_bytes(poly):
        return b"".join(struct.pack("<H", coeff % Q) for coeff in poly)

    def vec_to_bytes(vec):
        return b"".join(poly_to_bytes(poly) for poly in vec)

    # c = H(w,M)
    c = hash_to_ball(message, vec_to_bytes(w), 60)
    #print(c)

    def sample_uniform_poly_fast() -> List[int]:
        """Fast uniform polynomial sampling from {0, 1}^n"""
        random_bytes = os.urandom(N // 8)
        poly = []
        for i in range(N):
            byte_idx = i // 8
            bit_idx = i % 8
            bit_val = (random_bytes[byte_idx] >> bit_idx) & 1
            poly.append(bit_val)
        return poly
    
    # Step 6: u ←$ CBD_κ, b ←$ {0, 1}^n (optimized)
    kappa = 1
    nonce_u = (2000 + attempt).to_bytes(2, 'little')
    shake_input = rho_prime + nonce_u
    shake_ctx = shake256_absorb_finalize(shake_input)
    B_u = shake256_squeeze(shake_ctx, 64 * kappa)
    u = sample_poly_cbd(B_u, kappa)
            
    b = sample_uniform_poly_fast()  # Use faster sampling
            
    # Step 7: ξⱼ = 2uⱼ if cⱼ = 0, 2uⱼ + 1 - 2bⱼ if cⱼ = 1
    xi = []
    for j in range(N):
        if c[j] == 0:
            xi_j = 0
        else:
            xi_j = sym_mod(2 * u[j] + 1 - 2 * b[j],4) 
        xi.append(xi_j)
            
    def poly_sub_sym(a: List[int], b: List[int], q: int = Q) -> List[int]:
        """Polynomial subtraction with symmetric modulus"""
        return [sym_mod(ai - bi, q) for ai, bi in zip(a, b)]
    # Step 8: h = (ξ - c)s using NTT multiplication (simplified)
    xi_minus_c = poly_sub_sym(xi, c, Q)
    #print(xi_minus_c)


    #The term h


    h=[]
    for i in range(l+k+1):
        h.append(poly_mult_ntt(xi_minus_c, s[i]))
    #print(h)

    import numpy as np

    #print("\n A·h mod 2Q=0")
    twoQ = 2 * Q
    term = matrix_vector_mult_ntt(A, h)
    #print(term)
    term = np.array(term) % twoQ
    #print("A·h mod 2Q computed.", term)

    #The term cs

    cs=[]
    for i in range(l+k+1):
        cs.append(poly_mult_ntt(c, s[i]))

    def poly_add_sym(a: List[int], b: List[int], q: int = Q) -> List[int]:
        """Add two polynomials with symmetric modulus reduction."""
        if len(a) != N or len(b) != N:
            raise ValueError("Polynomials must have exactly N coefficients")
        return [sym_mod(ai + bi, q) for ai, bi in zip(a, b)]


    def vec_add_sym(vec1: List[List[int]], vec2: List[List[int]], q: int = Q) -> List[List[int]]:
        """Add two vectors of polynomials with symmetric modulus reduction."""
        if len(vec1) != len(vec2):
            raise ValueError("Polynomial vectors must have the same length")
        return [poly_add_sym(p1, p2, q) for p1, p2 in zip(vec1, vec2)]
    
    #Signature z=y+cs+h

    z = vec_add_sym(vec_add_sym(y,cs,Q),h,Q)
    # print(z)
    

    def poly_norm_infinity(poly: List[int]) -> int:
        """Calculate the infinity norm of a polynomial"""
        max_coeff = 0
        for coeff in poly:
            # Convert to signed representation
            signed_coeff = coeff if coeff <= Q // 2 else coeff - Q
            max_coeff = max(max_coeff, abs(signed_coeff))
        return max_coeff

    def vector_norm_infinity(vector: List[List[int]]) -> int:
        """Calculate the infinity norm of a vector of polynomials"""
        max_norm = 0
        for poly in vector:
            poly_norm = poly_norm_infinity(poly)
            max_norm = max(max_norm, poly_norm)
        return max_norm
    
    #Norm of z
    z_norm = vector_norm_infinity(z)
    print("Norm of z", z_norm)


############################################################################################################################
####################################                  Verify                  ##############################################
############################################################################################################################

Az = []
for i in range(k):
    # Initialize result polynomial for row i
    az_i = [0] * N

    for j in range(len(A[i])):  # A[i] has (1 + l + k) polynomials
        # Multiply the j-th polynomial of A[i] with the j-th polynomial of z
        product = poly_mult_ntt(A[i][j], z[j])
        az_i = poly_add(az_i, product)

    Az.append(az_i)

print(Az)

e_0 = [[0] * N for _ in range(k)]
e_0[0][0] = 1  # [1, 0, ..., 0] in first polynomial

# Compute ce_0 only for the first polynomial (others are zero)
ce_0 = poly_mult_ntt(c, e_0[0])

# Scale by Q
qce_0 = poly_scalar_mult(ce_0, Q)

# Now subtraction works since both are polynomials of length N
Az_minus_qcj = []
for i in range(k):
    diff = poly_sub(Az[i], qce_0)   # both are length-N polynomials
    diff_mod_2q = [coeff % (2 * Q) for coeff in diff]
    Az_minus_qcj.append(diff_mod_2q)

# print(Az_minus_qcj)

c_prime = hash_to_ball(message, vec_to_bytes(Az_minus_qcj), 60)

c_match = all(c[i] == c_prime[i] for i in range(N))
print(f"Challenge comparison: c == c' ? {c_match}")
print(f"  - Norm bound: {z_norm} < {300} = {z_norm < 300}")
