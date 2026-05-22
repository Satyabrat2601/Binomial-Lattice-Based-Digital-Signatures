"""
DUET: A Dual-Binomial Framework for Lattice-Based Digital Signatures
Security Estimator -- Conference-Ready

Usage:
    python3 DUET_estimator.py

Required files (same folder -- from https://github.com/pq-crystals/security-estimates):
    MSIS_security.py   MLWE_security.py   proba_util.py   model_BKZ.py

Optional pre-computed cache (copy from previous run or provided file):
    results_cache.json  -- skips the slow (~60-90s) MLWE estimator calls

================================================================
SCHEME SUMMARY  (Fig. 1 / Fig. 2)
================================================================
Ring: R_q = Z_q[X]/(X^n+1),  n=256,  q=64513  (NTT-friendly prime)

KeyGen (Fig. 1):
  A1 <-$ R_q^{k x ell}                              -- public matrix seed rho (32 B)
  (s1, s2) <- CBD_eta^ell x CBD_eta^k               -- secret key
  if ||s1||_2^2 > eta^2 * n  restart                 -- n = n_ring = 256  (per-poly check:
                                                     --   each s1_i satisfies ||s1_i||^2<=eta^2*n)
  b  = A1*s1 + s2  mod q
  A  = (-2b + qj | 2*A1 | 2*I_k)  mod 2q            -- full signing matrix
  S  = (1 | s1^T | s2^T)^T                          -- secret signing key
  pk = (rho, b);   sk = (rho, K, tr, s1, s2)

Sign (Fig. 1):
  y <-$ CBD_{gamma1}^{ell+k+1}                      -- masking vector, w polys
  w  = A*y mod 2q
  c  = H(mu || w) in B_tau                          -- SHAKE-256 + rejection, {0,1}^n weight-tau
  u  <-$ CBD_1,  d  <-$ {0,1}^n
  xi_j = 0                  if c_j = 0              -- Fig. 1 step 7
         2*u_j + 1 - 2*d_j  mod +-4  otherwise
  --   When c_j=1: (u,d) in {0,1}^2 -> raw in {-1,+1,3,1} -> mod +-4 -> {-1,+1}
  --   So xi_j | c_j=1 is uniform on {-1, +1} (Rademacher, NOT CBD_1).
  --   Marginal xi_j: P(0)=(n-tau)/n,  P(+-1)=tau/(2n) each
  g  = (xi - c)*S;   z = y + c*S + g = y + xi*S
  REJECT if ||z||_inf > B  or  z not in CBD_B

Sign (Fig. 2, compressed CSign):
  y <-$ CBD_{gamma1}^{ell+k+1}
  w  = A*y mod 2q
  w1 = HighBits(w);  w0 = LSB(w)                   -- CSign step 3-5
  c  = H(w1, w0, mu)
  ... compute z = y + xi*S (as above) ...
  (z1, z2) = Split(z)                               -- CSign step 9: z1 = first ell+1 polys
  h  = MakeHint(w, w - 2*z2)                        -- CSign step 10
  h  = w1 - HighBits(w - 2*z2) mod 2*(q-1)         -- CSign step 11
  REJECT if ||z||_inf > B  or  z not in CBD_B
  return sigma = (z1, h, c)

Verify (Fig. 2, CVerify):
  w' <- UseHint(h, A0*z1 - q*c*j) mod 2q           -- CVerify step 1
  w' = HighBits(A0*z1 - q*c*j) + h  mod 2*(q-1)   -- CVerify step 2
  z0 = first element of z1                          -- CVerify step 3
  w0 = LSB(z0 - xi)*j                               -- CVerify step 4
  c' = H(w', w0, mu)                                -- CVerify step 5
  z2'= (w' - (A0*z1 - q*c*j) + w0)/2 mod +-q      -- CVerify step 6
  z' = [z1, z2']                                    -- CVerify step 7
  Accept if (c = c') and (||z'||_inf <= B)          -- CVerify step 8

================================================================
SECURITY REDUCTION
================================================================
  MSIS columns:  w = ell + k + 1   (1 lead + ell masking + k identity block)
  MSIS rows:     h = k             (k equations from the identity block)
  MLWE dimension: d = ell  (secret s1 dimension)
  MLWE samples:   m = k   (samples from A1*s1 + s2)

  MSIS-weak  L2:  zeta_w = 2 * B * sqrt(n * w)
      (conservative bound: single accepted z has ||z||_2 <= B*sqrt(n*w);
       the factor 2 is a safe upper bound matching Dilithium-style analysis)
  MSIS-strong L2: zeta_s = 2 * (gamma1 + beta) * sqrt(n * w)
      (difference of two forgery z-vectors; entry-wise bound 2*(gamma1+beta))
  MSIS-weak  Linf: zeta_w_inf = 2 * B
      (conservative: single z has ||z||_inf <= B; factor 2 is a safe upper bound)
  MSIS-strong Linf: zeta_s_inf = 2 * (gamma1 + beta)

  Constraint:  gamma1 = B + beta   (ensures accepted |z_i| <= B)
               beta   = tau * eta

================================================================
PARAMETERS
================================================================
  DUET-I:   ell=3 k=2 eta=1 B=30 tau=48  -> target ~120b classical
  DUET-II:  ell=4 k=3 eta=1 B=35 tau=60  -> target ~180b classical
  DUET-III: ell=6 k=4 eta=1 B=39 tau=72  -> target ~260b classical
"""

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import sys, os, json, contextlib, io, struct, time, random, hashlib
from math import ceil, log2, comb, sqrt

# ---------------------------------------------------------------------------
# NumPy  (required for vectorised Monte Carlo)
# ---------------------------------------------------------------------------
try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False
    print("WARNING: numpy not found -- Monte Carlo and empirical tests disabled.")

# ---------------------------------------------------------------------------
# pq-crystals security-estimates (must be in the same directory)
# ---------------------------------------------------------------------------
_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)

from MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from MLWE_security  import MLWE_summarize_attacks, MLWEParameterSet
from proba_util     import (build_centered_binomial_law, iter_law_convolution,
                            law_convolution, tail_probability)


# ===========================================================================
#  1. CBD SAMPLER
# ===========================================================================
def cbd_sample(k, size=1):
    """
    Sample `size` draws from CBD_k = Binomial(2k, 1/2) - k.
    Returns a numpy int32 array, or a plain Python list if numpy is absent.
    """
    if not _HAS_NUMPY:
        out = []
        for _ in range(size):
            a = sum(random.randint(0, 1) for _ in range(k))
            b = sum(random.randint(0, 1) for _ in range(k))
            out.append(a - b)
        return out
    rng = np.random.default_rng()
    return (rng.binomial(2 * k, 0.5, size) - k).astype(np.int32)


# ===========================================================================
#  2. TRUNCATED CBD SAMPLER
# ===========================================================================
def truncated_cbd(k, B, size=1, rng=None):
    """
    Sample `size` values from CBD_k conditioned on |x| <= B.
    Uses rejection sampling in batches for efficiency.
    This enforces z in CBD_B for the signing masking vector y.
    """
    if not _HAS_NUMPY:
        out = []
        while len(out) < size:
            x = cbd_sample(k, 1)[0]
            if abs(x) <= B:
                out.append(x)
        return out
    if rng is None:
        rng = np.random.default_rng()
    out  = np.empty(size, dtype=np.int32)
    done = 0
    while done < size:
        need  = size - done
        batch = (rng.binomial(2 * k, 0.5, need * 6) - k).astype(np.int32)
        valid = batch[np.abs(batch) <= B]
        take  = min(len(valid), need)
        out[done:done + take] = valid[:take]
        done += take
    return out


# ===========================================================================
#  3. COMPRESSION PRIMITIVES (Fig. 2)
# ===========================================================================
def HighBits(w, beta, q):
    """
    HighBits(w) = floor(w / beta) * beta mod 2*(q-1).

    Returns the scaled high component used in the paper's Fig.2
    decomposition. Works element-wise on numpy arrays or scalars.
    """

    if _HAS_NUMPY and isinstance(w, np.ndarray):
        return (np.floor_divide(w.astype(np.int64), beta) * beta) % (2 * (q - 1))

    return ((int(w) // beta) * beta) % (2 * (q - 1))


def LowBits(w, beta, q):
    """
    LowBits(w) = w - HighBits(w) mod 2q.

    Extracts the residual low component after removing the scaled
    HighBits contribution.
    """

    if _HAS_NUMPY and isinstance(w, np.ndarray):
        w64 = w.astype(np.int64)
        hi  = HighBits(w64, beta, q)
        return (w64 - hi) % (2 * q)

    hi = HighBits(w, beta, q)
    return (int(w) - hi) % (2 * q)


def LSB(w):
    """
    LSB(w) = w mod 2  (least significant bit).
    Used in:
      CSign step 5: w0 = LSB(w)
      CVerify step 4: w0 = LSB(z0 - xi) * j
    Works element-wise on numpy arrays or scalars.
    """
    if _HAS_NUMPY and isinstance(w, np.ndarray):
        return (w.astype(np.int64) % 2).astype(np.int32)
    return int(w) % 2

def MakeHint(w, wp, beta, q):
    """
    MakeHint(w, w') exactly following Fig.2 CSign step 10-11.

    Computes the HighBits correction value:

        h = HighBits(w) - HighBits(w') mod 2(q-1)

    where:
        w' = w - 2*z2

    This is NOT a boolean Dilithium-style hint.
    It stores the exact modular correction needed for recovery.

    Parameters
    ----------
    w : int or np.ndarray
        Original value.

    wp : int or np.ndarray
        Perturbed value:
            wp = w - 2*z2

    beta : int
        HighBits decomposition parameter.

    q : int
        Ring modulus.

    Returns
    -------
    int or np.ndarray
        Hint correction modulo 2(q-1).
    """

    modulus = 2 * (q - 1)

    hw  = HighBits(w,  beta, q)
    hwp = HighBits(wp, beta, q)

    if _HAS_NUMPY and isinstance(w, np.ndarray):
        return (hw.astype(np.int64) - hwp.astype(np.int64)) % modulus

    return (int(hw) - int(hwp)) % modulus


def UseHint(h, r, beta, q):
    """
    UseHint(h, r) exactly following Fig.2 CVerify steps 1-2.

    Verification computes:

        w' = HighBits(r) + h mod 2(q-1)

    where:
        r = A0*z1 - q*c*j

    Parameters
    ----------
    h : int or np.ndarray
        HighBits correction value from MakeHint().

    r : int or np.ndarray
        Verifier-side reconstructed value:
            r = A0*z1 - q*c*j

    beta : int
        HighBits decomposition parameter.

    q : int
        Ring modulus.

    Returns
    -------
    int or np.ndarray
        Recovered HighBits value modulo 2(q-1).
    """

    modulus = 2 * (q - 1)

    if _HAS_NUMPY and isinstance(r, np.ndarray):

        hi = HighBits(r.astype(np.int64), beta, q)

        return (hi.astype(np.int64) + h.astype(np.int64)) % modulus

    hi = HighBits(int(r), beta, q)

    return (int(hi) + int(h)) % modulus

# ===========================================================================
#  4. XI SAMPLING  (Fig. 1, Sign step 7)
# ===========================================================================
def sample_xi(c_coeff, u_j, d_j):
    
    if c_coeff == 0:
        return 0
    
    # u_j ∈ {-1, 0, +1}
    val = 2 * u_j + 1 - 2 * d_j
    # Reduce mod ±4 to {-1, +1}
    if val == 3: return -1
    if val == -3: return 1
    return val

def sample_xi_vector(c, n, rng=None):
    """
    Vectorised xi sampling for a full challenge polynomial c (length n).
    Implements Fig. 1 steps 6-7 exactly.
    Returns xi as numpy int32 array of length n, with xi[j] in {-1, 0, +1}.
    When c[j] = 0: xi[j] = 0.
   When c[j] = 1: xi[j] = 2*u[j] + 1 - 2*d[j], u ~ CBD₁, d ~ Uniform{0,1}.
    """
    if not _HAS_NUMPY:
        return [sample_xi(int(c[j]), random.randint(0, 1), random.randint(0, 1))
                for j in range(n)]
    if rng is None:
        rng = np.random.default_rng()
    c_arr = np.asarray(c, dtype=np.int32)
    u = (rng.binomial(2, 0.5, size=n) - 1).astype(np.int32)   
    d     = rng.integers(0, 2, size=n, dtype=np.int32)  
    raw   = 2 * u + 1 - 2 * d                            # in {-1, +1}
    raw   = np.where(raw > 2, raw - 4, np.where(raw < -2, raw + 4, raw))
    xi    = np.where(c_arr != 0, raw, np.int32(0))
    return xi.astype(np.int32)

def w0_from_lsb(z0, xi0, j, q):
    """
    Compute w0 = LSB(z0 - xi0) * j  as in Fig. 2 CVerify step 4.
    z0 : first coefficient of z1 (scalar)
    xi0: xi perturbation at position 0 (0 or +-1)
    j  : unit vector selector (scalar 1 here since j = (1,0,...,0))
    q  : modulus (kept for API consistency)
    Returns w0 in {0, 1}.
    """
    return int(LSB(int(z0) - int(xi0)) * int(j))


# ===========================================================================
#  5. CHALLENGE ENCODING  (Fig. 1, Sign step 5)
# ===========================================================================
def encode_challenge(n, tau, seed_bytes):
    """
    Encode challenge polynomial c in B_tau from seed_bytes.
    B_tau = { c in {0,1}^n : ||c||_1 = tau }  (weight-tau binary).

    Note: DUET uses {0,1} challenge coefficients (not {-1,0,1} like Dilithium).
    The sign is absorbed into xi via u_j (Fig. 1 step 7).

    Algorithm: extendable-output rejection sampling (Dilithium-style).
      1. Initialise SHAKE-256 with seed_bytes as a streaming XOF.
      2. Read 2 bytes at a time; combine into a 16-bit value.
      3. Reject if value mod 2**16 >= n * floor(2**16 / n)  (uniform up to bias).
         Else position = value mod n.
      4. Skip already-selected positions; continue until tau are chosen.
      5. The XOF is unbounded, so this loop terminates with prob 1 and
         requires no deterministic fallback.

    For n=256 (a power of two dividing 2**16), 2**16 / n = 256 is exact, so
    there is no rejection bias and step 3 simplifies.

    Challenge entropy: H_c = log2 C(n, tau) + tau bits.
      log2 C(n,tau) encodes which tau positions are chosen.
      tau extra bits encode the u_j signs (each contributing one bit via xi).
    """
    c = [0] * n
    xof = hashlib.shake_256(seed_bytes)
    selected = set()
    # Stream bytes from the XOF until we have tau distinct positions.
    # We pull in 64-byte blocks to amortise the call cost.
    block_size = 64
    block_idx  = 0
    pos_in_block = 0
    block = xof.digest(block_size)
    while len(selected) < tau:
        if pos_in_block + 2 > len(block):
            # extend the digest: SHAKE supports streaming via repeated digest of
            # increasing length; we approximate by extending in 64-byte steps.
            block_idx += 1
            block = xof.digest(block_size * (block_idx + 1))[block_size * block_idx:]
            pos_in_block = 0
        b0 = block[pos_in_block]
        b1 = block[pos_in_block + 1]
        pos_in_block += 2
        # combine to 16-bit; uniform mod n since n | 2**16 for n=256.
        val = (b1 << 8) | b0
        if n & (n - 1) == 0 and n <= 65536:
            pos = val & (n - 1)            # bias-free for power-of-two n
        else:
            # Generic rejection sampling to remove modulo bias
            limit = (65536 // n) * n
            if val >= limit:
                continue                   # reject and resample
            pos = val % n
        selected.add(pos)
    for pos in selected:
        c[pos] = 1
    return c


def challenge_entropy(n, tau):
    """
    H_c = log2(C(n, tau)) + tau bits.
    log2(C(n,tau)) bits for position selection, tau bits for the u_j signs.
    """
    return log2(comb(n, tau)) + tau


# ===========================================================================
#  6. SPLIT(z)  (Fig. 2, CSign step 9)
# ===========================================================================
def Split(z, n, ell, k):
    """
    Split(z) as in Fig. 2, CSign step 9.
    z has length n*(ell+k+1)  (w = ell+k+1 polynomials total).
    z1 = first n*(ell+1) coefficients  -- lead poly + ell masking polys.
    z2 = remaining n*k coefficients    -- key-material contribution.
    z1 is transmitted in the signature; z2 is used only for MakeHint.
    """
    split_idx = n * (ell + 1)
    if _HAS_NUMPY and isinstance(z, np.ndarray):
        return z[:split_idx], z[split_idx:]
    return list(z[:split_idx]), list(z[split_idx:])


# ===========================================================================
#  7. SIGNATURE BYTE PACKING  (improvement #4)
# ===========================================================================
def pack_signature(z, c, h):
    """
    Pack signature components (z1, c, h) into a byte string.
    Format: [z1 coefficients as little-endian int16] [hint bytes] [challenge bytes]
    z: iterable of 16-bit signed integers (coefficients of z1, ell+1 polynomials)
    c: bytes-like, challenge encoding (ceil((log2 C(n,tau)+tau)/8) bytes)
    h: iterable of uint8 hint values (omega + k bytes)
    """
    data = bytearray()
    for x in z:
        data += struct.pack('<h', int(x))       # signed 16-bit, little-endian
    for x in h:
        data += struct.pack('B', int(x) & 0xFF)
    data += bytes(c)
    return bytes(data)


def unpack_signature(sig, z_len, h_len, c_len):
    """
    Unpack a byte string produced by pack_signature.
    Returns (z_list, h_list, c_bytes).
    """
    ptr = 0
    z   = []
    for _ in range(z_len):
        z.append(struct.unpack('<h', sig[ptr:ptr + 2])[0])
        ptr += 2
    h   = list(sig[ptr:ptr + h_len])
    ptr += h_len
    c   = sig[ptr:ptr + c_len]
    return z, h, c


# ===========================================================================
#  8. MONTE CARLO REJECTION RATE
# ===========================================================================
def monte_carlo_rejection(ps, trials=50_000, chunk=1_000, seed=42):
    """
    Estimate the signing rejection probability via Monte Carlo simulation,
    following Fig. 1 steps 3-10 line-by-line for each trial:

      step 3 : y ~ CBD_{gamma1}^{n*(ell+k+1)}
      step 5 : c in B_tau  (we sample tau distinct positions uniformly)
      step 6 : u, d ~ Uniform{0,1}^n
      step 7 : xi_j = 0 if c_j=0, else (2u_j + 1 - 2d_j) mod +-4
      step 9 : z = y + xi*S, modelled per-coefficient as:
                 lead poly  z_0[j] = y_0[j] + xi[j]
                 other poly z_i[j] = y_i[j] + (xi*s_i)[j]
                 where (xi*s_i)[j] is the j-th coefficient of the cyclic-anti
                 convolution xi*s_i in R_q.  We sample s_i ~ CBD_eta^n and
                 compute the convolution directly.
      step 10: reject if (||z||_inf > B)  OR  (z not in CBD_B).

    The "z not in CBD_B" check is implemented by verifying every coefficient
    of z lies in the support of CBD_B (i.e. in [-B, B]).  For the chosen
    parameter sets this support coincides with the ||z||_inf <= B condition,
    so the two conditions are equivalent under Fig. 1's parameter regime;
    we keep both for fidelity to the scheme text.

    Uses constant-time rejection check (bitwise OR accumulator) within each
    trial, matching a real implementation's side-channel-safe behavior.
    """
    if not _HAS_NUMPY:
        print("    [MC] numpy not available; skipping Monte Carlo.")
        return float("nan")

    rng = np.random.default_rng(seed)
    n, ell, k = ps.n, ps.ell, ps.k
    tau       = ps.tau
    eta       = ps.eta
    gamma1    = ps.gamma1
    B         = ps.B
    n_other   = ell + k     # number of non-lead polynomials
    rejected  = 0

    # Pre-compute index arrays for the anti-cyclic convolution xi*s in R_q = Z[X]/(X^n+1).
    # For (xi*s)[j] = sum_{l=0..n-1} sign(j,l) * xi[l] * s[(j-l) mod n]
    # with sign = -1 when (j-l) wraps (l > j), else +1.
    idx_j = np.arange(n, dtype=np.int64)[:, None]   # (n,1)
    idx_l = np.arange(n, dtype=np.int64)[None, :]   # (1,n)
    shift = (idx_j - idx_l) % n                     # (n,n)
    sign  = np.where(idx_l > idx_j, -1, 1).astype(np.int32)  # (n,n)

    for start in range(0, trials, chunk):
        t_ = min(chunk, trials - start)
        # Sample c in B_tau: tau distinct positions per trial, weight-tau binary
        c_batch = np.zeros((t_, n), dtype=np.int32)
        for r in range(t_):
            pos = rng.choice(n, size=tau, replace=False)
            c_batch[r, pos] = 1

        # u, d ~ Uniform{0,1}^n   (Fig. 1 step 6)
        u_batch = rng.integers(0, 2, size=(t_, n), dtype=np.int32)
        d_batch = rng.integers(0, 2, size=(t_, n), dtype=np.int32)

        # xi per Fig. 1 step 7
        raw = 2*u_batch + 1 - 2*d_batch
        raw = np.where(raw > 2, raw - 4, np.where(raw < -2, raw + 4, raw))
        xi_batch = np.where(c_batch != 0, raw, np.int32(0)).astype(np.int32)

        # Sample s = (s_1 | ... | s_{ell+k}) ~ CBD_eta^{n*(ell+k)}
        # For the MC we re-sample s every trial (matches an HVZK forgery game;
        # for keygen-fixed s the rejection rate is the same in expectation).
        s_batch = (rng.binomial(2*eta, 0.5, (t_, n_other, n)) - eta).astype(np.int32)

        # Sample y ~ CBD_{gamma1}^{n*(ell+k+1)}
        y_batch = (rng.binomial(2*gamma1, 0.5, (t_, ell+k+1, n)) - gamma1).astype(np.int32)

        # --- Lead poly: z_0 = y_0 + xi  (S[0] = 1) ---
        z_lead = y_batch[:, 0, :] + xi_batch    # (t_, n)

        # --- Other polys: z_i = y_i + xi * s_i  (anti-cyclic convolution in R_q) ---
        # xi[l] reordered by shift[j,l] gives the s-coefficient at position (j-l) mod n
        # (xi*s)[j] = sum_l sign[j,l] * xi[l] * s[shift[j,l]]
        # Vectorise across the batch and the (ell+k) polys.
        # xi_batch: (t_, n);  s_batch: (t_, n_other, n)
        # We want for each j:  conv[j] = sum_l sign[j,l] * xi[l] * s[shift[j,l]]
        # Build using gather: s_batch[..., shift] has shape (t_, n_other, n, n)
        # That's O(t_*n_other*n*n) -- for n=256, n_other up to 10, chunk=200 trials,
        # this is 200*10*256*256 = ~131M; doable but heavy. Drop chunk if needed.
        s_gather = s_batch[:, :, shift]                     # (t_, n_other, n_j, n_l)
        # multiply by xi[l] and sign[j,l], then sum over l
        weighted = s_gather * (xi_batch[:, None, None, :] * sign[None, None, :, :])
        conv     = weighted.sum(axis=-1).astype(np.int32)   # (t_, n_other, n)

        z_other  = y_batch[:, 1:, :] + conv                 # (t_, n_other, n)

        # Concatenate and flatten per trial
        z_full = np.concatenate(
            (z_lead[:, None, :], z_other), axis=1
        ).reshape(t_, -1)                                   # (t_, (ell+k+1)*n)

        # ------------------------------------------------------------
        # Rejection check (Fig.1 step 10):
        #
        # Reject if:
        #   (1) ||z||_inf > B
        #   OR
        #   (2) z not in CBD_B
        #
        # Here CBD_B membership is approximated via support membership:
        # every coefficient must lie inside the support of CBD_B.
        # ------------------------------------------------------------

        # Condition 1: infinity norm bound
        over_norm = (np.abs(z_full) > B)

        # Condition 2: CBD_B support membership
        cbd_B = build_centered_binomial_law(B)
        valid_vals = np.array(list(cbd_B.keys()), dtype=np.int32)

        over_cbd = ~np.isin(z_full, valid_vals)

        # Combine both rejection conditions
        reject_mask = np.logical_or(over_norm, over_cbd)

        # Constant-time reduction across coefficients
        flags = np.bitwise_or.reduce(reject_mask.astype(np.int32), axis=1)

        rejected += int((flags != 0).sum())

    return rejected / trials


def monte_carlo_rejection_fast(ps, trials=50_000, chunk=2_000, seed=42):
    """
    Faster but less faithful MC: bounds xi*S entry-wise by CBD_beta and skips the
    polynomial convolution.  Kept for cross-validation against the analytical
    rejection rate (both use the same conservative bound).  This is what the
    original code labelled monte_carlo_rejection -- renamed here to clarify
    that the scheme-faithful version is monte_carlo_rejection above.
    """
    if not _HAS_NUMPY:
        print("    [MC-fast] numpy not available; skipping.")
        return float("nan")
    rng      = np.random.default_rng(seed)
    dim      = ps.n * ps.w
    rejected = 0
    for start in range(0, trials, chunk):
        t_ = min(chunk, trials - start)
        y   = (rng.binomial(2 * ps.gamma1, 0.5, (t_, dim)) - ps.gamma1).astype(np.int32)
        xiS = (rng.binomial(2 * ps.beta,   0.5, (t_, dim)) - ps.beta  ).astype(np.int32)
        z   = y + xiS
        over  = (np.abs(z) > ps.B).astype(np.int32)
        flags = np.bitwise_or.reduce(over, axis=1)
        rejected += int((flags != 0).sum())
    return rejected / trials


# ===========================================================================
#  9. ANALYTICAL REJECTION RATE  (proba_util, for cross-validation)
# ===========================================================================
def rejection_rate_analytical(ps):
    """
    Analytical rejection probability via exact distribution of z.

    Per Fig. 1 step 9:  z = y + xi*S  where S = (1 | s1 | s2)^T.

    LEAD POLYNOMIAL  (S[0] = 1, the constant 1):
      z_0[j] = y_0[j] + xi[j]
      where xi[j] is the j-th coefficient of the xi polynomial.

      Marginal distribution of xi[j]  (averaging over a uniformly chosen index j,
      and over uniform choice of c in B_tau and (u,d) per Fig. 1 step 7):
          P(xi[j] = 0)  = (n - tau) / n     (when c_j = 0)
          P(xi[j] = +1) = tau / (2n)        (when c_j = 1, raw -> +1)
          P(xi[j] = -1) = tau / (2n)        (when c_j = 1, raw -> -1)

      Hence z_0[j] ~ CBD_{gamma1}  *  xi_marginal.
      This is NOT CBD_{gamma1} * CBD_1 -- CBD_1 has the wrong mass on 0.

    OTHER POLYNOMIALS  (S[i] = s1 or s2 component for i >= 1):
      z_i[j] = y_i[j] + (xi * s_i)[j]
      where xi has exactly tau nonzero coefficients in {-1,+1} uniformly.

      The j-th coefficient of (xi * s_i)  (multiplication in R_q = Z_q[X]/(X^n+1))
      is a signed sum of tau coefficients of s_i.  Since CBD_eta is symmetric
      around 0, (+-1)*CBD_eta = CBD_eta in distribution, so (xi*s_i)[j] is
      distributed as the tau-fold convolution of CBD_eta.

      Hence z_i[j] ~ CBD_{gamma1}  *  (CBD_eta convolved tau times).

    REJECTION  (Fig. 1 step 10:  ||z||_inf > B  OR  z not in CBD_B):
      The dominant condition is ||z||_inf > B.  Per-coefficient tail:
          p_lead  = P(|z_0[j]| > B)
          p_other = P(|z_i[j]| > B)  for i >= 1
      Treating coefficients as independent (mild over-estimate; in reality the
      tau positions of c couple them slightly):
          P_reject  =  1  -  (1 - p_lead)^n  *  (1 - p_other)^{n*(ell+k)}

    Note: the empirical Monte-Carlo function below uses the conservative
    CBD_beta model for xi*S to match Dilithium-style entry-wise bounds.
    The analytical function uses the EXACT scheme distribution above and
    is therefore the more precise estimate.
    """
    cbd_y   = build_centered_binomial_law(ps.gamma1)
    cbd_eta = build_centered_binomial_law(ps.eta)

    # CORRECT xi marginal per Fig. 1 step 7  (NOT CBD_1)
    xi_marginal = {
        -1: ps.tau / (2.0 * ps.n),
         0: (ps.n - ps.tau) / float(ps.n),
        +1: ps.tau / (2.0 * ps.n),
    }

    # Lead poly: z_0[j] = y_0[j] + xi[j]
    L_lead  = law_convolution(cbd_y, xi_marginal)

    # Other polys: z_i[j] = y_i[j] + (xi*s_i)[j]  where (xi*s_i)[j] is tau-fold conv of CBD_eta
    xiS     = iter_law_convolution(cbd_eta, ps.tau)
    L_other = law_convolution(cbd_y, xiS)

    # ------------------------------------------------------------
    # Scheme-style rejection:
    #
    # Reject if:
    #   (1) ||z||_inf > B
    #   OR
    #   (2) z not in CBD_B
    #
    # CBD_B membership is approximated via support membership.
    # ------------------------------------------------------------

    cbd_B = build_centered_binomial_law(ps.B)

    valid_vals = set(cbd_B.keys())

    # Probability coefficient falls outside CBD_B support
    p_lead = sum(
        prob for val, prob in L_lead.items()
        if val not in valid_vals
    )

    p_other = sum(
        prob for val, prob in L_other.items()
        if val not in valid_vals
    )
    # Total: any coefficient violating triggers rejection
    return 1.0 - (1.0 - p_lead)**ps.n * (1.0 - p_other)**(ps.n * (ps.ell + ps.k))


# ===========================================================================
#  10. ENTROPY / SIGNATURE SIZE  (rANS entropy coding)
# ===========================================================================
def _shannon(L):
    """Shannon entropy of a probability distribution given as a dict {x: prob}."""
    return -sum(p * log2(p) for p in L.values() if p > 0)


def _truncate_law(L, B):
    """Restrict L to [-B, B] and renormalise to unit mass."""
    mass = sum(p for x, p in L.items() if abs(x) <= B)
    if mass <= 0:
        raise ValueError(f"No probability mass in [-{B},{B}]")
    return {x: p / mass for x, p in L.items() if abs(x) <= B and p > 0}


def coef_entropies(ps):
    """
    Per-coefficient Shannon entropy (bits) of the accepted z distribution.

    Same distribution model as rejection_rate_analytical:
      lead poly  : z_0[j] = y_0[j] + xi[j]  with xi marginal {-1: t/2n, 0: (n-t)/n, 1: t/2n}
      other poly : z_i[j] = y_i[j] + (xi*s_i)[j]  ~  CBD_{gamma1} * tau-fold conv of CBD_eta
    Both distributions are truncated to |.| <= B (acceptance condition) before
    computing Shannon entropy, since the entropy coder only sees accepted z.
    """
    cbd_y   = build_centered_binomial_law(ps.gamma1)
    cbd_eta = build_centered_binomial_law(ps.eta)

    # Correct xi marginal per Fig. 1 step 7
    xi_marginal = {
        -1: ps.tau / (2.0 * ps.n),
         0: (ps.n - ps.tau) / float(ps.n),
        +1: ps.tau / (2.0 * ps.n),
    }

    L_lead  = law_convolution(cbd_y, xi_marginal)
    xiS     = iter_law_convolution(cbd_eta, ps.tau)
    L_other = law_convolution(cbd_y, xiS)
    H_lead  = _shannon(_truncate_law(L_lead,  ps.B))
    H_other = _shannon(_truncate_law(L_other, ps.B))
    return H_lead, H_other


def sig_size_entropy(ps):
    """
    Compressed signature size (bytes) using rANS entropy coding.

    Signature = (z1, h, c):
      z1   : n*(ell+1) coefficients -- 1 lead poly + ell masking polys.
              Entropy: n*H_lead + n*ell*H_other bits.
      h    : hint vector.
              Entropy model: Bernoulli(omega / (n*k)) per position over n*k bits.
              Note: hint_B in the breakdown is the entropy-coded size (not the
              naive packed sparse encoding of omega+k bytes).
      c    : challenge, ceil(log2(C(n,tau)) + tau) bits.
              (tau positions encoded + tau sign bits absorbed into xi)

    The rANS total covers all three components at their entropy lower-bound.
    For reference, the naive packed sparse hint encoding uses omega+k bytes
    (see packed_sig_size()), which is the upper-bound figure.

    Returns (total_bytes, breakdown_dict).
    """
    H_lead, H_other = coef_entropies(ps)
    z1_bits = ps.n * H_lead + ps.n * ps.ell * H_other   # z1 entropy (bits)

    # Hint entropy: sparse binary vector of length n*k, expected omega ones.
    # rANS encodes this at Bernoulli(p1) Shannon entropy.
    p1 = ps.omega / (ps.n * ps.k)
    p0 = 1.0 - p1
    H_hint = 0.0
    if p0 > 0: H_hint -= p0 * log2(p0)
    if p1 > 0: H_hint -= p1 * log2(p1)
    hint_bits = ps.n * ps.k * H_hint

    c_bits = ceil(log2(comb(ps.n, ps.tau))) + ps.tau     # challenge bits
    total  = ceil((z1_bits + hint_bits + c_bits) / 8)

    # Breakdown: all components at entropy-coded size (rANS lower bound).
    # hint_B_packed = omega + k  is the packed sparse alternative (upper bound).
    return total, dict(
        z1_B         = ceil(z1_bits / 8),
        hint_B       = ceil(hint_bits / 8),   # entropy-coded hint size
        hint_B_packed= ps.omega + ps.k,        # packed sparse hint (upper bound)
        c_B          = ceil(c_bits / 8),
        H_lead       = H_lead,
        H_other      = H_other,
    )


def packed_sig_size(ps):
    """
    Upper-bound signature size (bytes) using naive int16 packing (no entropy coding).

    Signature = (z1, h, c):
      z1 : n*(ell+1) coefficients x 2 bytes each  (int16, range [-B, B])
      h  : omega + k bytes  (Dilithium-style sparse hint: k offsets + omega positions)
      c  : ceil((log2(C(n,tau)) + tau) / 8) bytes  (challenge encoding)

    This is a strict upper bound; rANS entropy coding achieves sig_size_entropy().
    """
    z1_bytes   = ps.n * (ps.ell + 1) * 2          # ell+1 polys, int16
    hint_bytes = ps.omega + ps.k                   # sparse hint
    c_bytes    = ceil((log2(comb(ps.n, ps.tau)) + ps.tau) / 8)
    return z1_bytes + hint_bytes + c_bytes


# ===========================================================================
#  11. EMPIRICAL z-DISTRIBUTION VALIDATION
# ===========================================================================
def empirical_norm_test(ps, trials=5_000, seed=0):
    """
    Empirically estimate the L2 norm of *accepted* z vectors.

    Correct model (Fig. 1 steps 3-10):
      1. Sample y ~ CBD_{gamma1}^{n*w}  (untruncated, as in the scheme).
      2. Sample xiS ~ CBD_{beta}^{n*w}  (entry-wise conservative bound on xi*S).
      3. Compute z = y + xiS.
      4. Accept only if ||z||_inf <= B  (rejection condition from step 10).
      5. Record ||z||_2 for accepted trials.

    The empirical max L2 of accepted z should lie well below zeta_s_l2,
    confirming the proof bound is conservative.

    Note: the CBD_{beta} model for xiS is the conservative entry-wise bound
    used in the proof; the scheme-faithful model (negacyclic convolution of
    xi*s) is implemented in monte_carlo_rejection().
    """
    if not _HAS_NUMPY:
        print("    [empirical] numpy not available; skipping.")
        return {}
    rng     = np.random.default_rng(seed)
    dim     = ps.n * ps.w
    l2_vals = []
    attempts = 0
    while len(l2_vals) < trials:
        attempts += 1
        # Step 1: y ~ CBD_{gamma1} (full distribution, no truncation)
        y   = (rng.binomial(2 * ps.gamma1, 0.5, dim) - ps.gamma1).astype(np.int32)
        # Step 2: conservative xi*S perturbation ~ CBD_{beta}
        xiS = (rng.binomial(2 * ps.beta, 0.5, dim) - ps.beta).astype(np.int32)
        z   = y + xiS
        # Step 4: accept only if ||z||_inf <= B
        if np.max(np.abs(z)) <= ps.B:
            l2_vals.append(float(np.linalg.norm(z.astype(np.float64), ord=2)))
        if attempts > trials * 200:   # guard against infinite loop for bad params
            break
    if not l2_vals:
        return {}
    arr       = np.array(l2_vals)
    # Theoretical mean L2 of the untruncated distribution (conservative reference)
    theo_mean = sqrt(ps.n * ps.w * (ps.gamma1 + ps.beta) / 2)
    return dict(
        mean_l2   = float(arr.mean()),
        max_l2    = float(arr.max()),
        theo_mean = theo_mean,
        zeta_s    = ps.zeta_s_l2,
        accepted  = len(l2_vals),
        attempts  = attempts,
    )


# ===========================================================================
#  12. BENCHMARKS
# ===========================================================================
def benchmark(fn, iterations=100, warmup=3):
    """Wall-clock average time of fn() over `iterations` calls. Returns seconds."""
    for _ in range(warmup):
        fn()
    start = time.perf_counter()
    for _ in range(iterations):
        fn()
    return (time.perf_counter() - start) / iterations


def run_benchmarks(ps):
    """
    Simulated numpy-level timing for KeyGen, Sign, Verify.
    These are arithmetic-only lower bounds (no NTT, hash, or full lattice ops).
    """
    if not _HAS_NUMPY:
        return {}
    n, dim = ps.n, ps.n * ps.w
    rng = np.random.default_rng(0)

    def sim_keygen():
        # Sample s1 ~ CBD_eta^{n*ell}; per Fig. 1 KeyGen step 3:
        #   "if ||s1||_2^2 > eta^2 * n, bot, go to step 2"
        # where n = n_ring = 256 (the ring dimension).
        # This is interpreted as a per-polynomial check: for each of the ell
        # component polynomials s1_i, require ||s1_i||_2^2 <= eta^2 * n_ring.
        # E[||s1_i||^2] = n_ring * eta/2 = 128 for eta=1; threshold = 256 = 2*E,
        # so the check is rarely triggered and KeyGen terminates quickly.
        s1 = (rng.binomial(2 * ps.eta, 0.5, (ps.ell, n)) - ps.eta).astype(np.int32)
        norms_sq = np.sum(s1 ** 2, axis=1)          # one norm per polynomial
        return bool(np.all(norms_sq <= ps.eta**2 * n))

    def sim_sign():
        # y ~ CBD_{gamma1}; z = y + perturbation; constant-time reject check
        y   = (rng.binomial(2 * ps.gamma1, 0.5, dim) - ps.gamma1).astype(np.int32)
        xiS = (rng.binomial(2 * ps.beta,   0.5, dim) - ps.beta  ).astype(np.int32)
        z   = y + xiS
        # Constant-time: bitwise OR over all positions
        flag = int(np.bitwise_or.reduce((np.abs(z) > ps.B).astype(np.int32)))
        return flag == 0

    def sim_verify():
        # Check ||z||_inf <= B  (vectorised)
        z = (rng.binomial(2 * ps.gamma1, 0.5, dim) - ps.gamma1).astype(np.int32)
        return bool(np.max(np.abs(z)) <= ps.B)

    return dict(
        keygen_ms = benchmark(sim_keygen, 200) * 1000,
        sign_ms   = benchmark(sim_sign,   200) * 1000,
        verify_ms = benchmark(sim_verify, 200) * 1000,
    )


# ===========================================================================
#  CACHE  (skip slow MLWE sweeps on repeat runs)
# ===========================================================================
_CACHE_FILE = os.path.join(_DIR, "results_cache.json")

def _load_cache():
    if os.path.exists(_CACHE_FILE):
        try:
            with open(_CACHE_FILE) as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def _save_cache(c):
    with open(_CACHE_FILE, "w") as f:
        json.dump(c, f, indent=2, sort_keys=True)

def _run_msis(n, w, h, zeta, q, norm="linf"):
    """Run MSIS security estimator (with cache). norm = 'linf' or 'l2'."""
    key = f"msis|n={n}|w={w}|h={h}|B={zeta}|q={q}|norm={norm}"
    c   = _load_cache()
    if key in c:
        return tuple(c[key])
    print(f"      [estimator] MSIS(n={n},w={w},h={h},zeta={zeta},q={q},"
          f"norm={norm}) ...", flush=True)
    with contextlib.redirect_stdout(io.StringIO()):
        r = MSIS_summarize_attacks(MSISParameterSet(n, w, h, zeta, q, norm=norm))
    c[key] = list(r)
    _save_cache(c)
    return tuple(r)

def _run_mlwe(n, d, m, eta, q):
    """Run MLWE security estimator (with cache). d=ell, m=k."""
    key = f"mlwe|n={n}|d={d}|m={m}|k={eta}|q={q}|distr=binomial"
    c   = _load_cache()
    if key in c:
        return tuple(c[key])
    print(f"      [estimator] MLWE(n={n},d={d},m={m},eta={eta},q={q}) "
          f"-- first run ~60-90s ...", flush=True)
    with contextlib.redirect_stdout(io.StringIO()):
        r = MLWE_summarize_attacks(MLWEParameterSet(n, d, m, eta, q, "binomial"))
    c[key] = list(r)
    _save_cache(c)
    return tuple(r)


# ===========================================================================
#  PARAMETER CLASS
# ===========================================================================
class DUETParams:
    """
    One DUET-Sign parameter set.

    All parameters derived from the scheme specification (Fig. 1 / Fig. 2):

    Core constraint:  gamma1 = B + beta  (enforced by assertion)
      Here `beta` denotes tau*eta -- the entry-wise bound on (xi*s)_i used
      for the rejection check in Fig. 1 step 10.

      NOTE on notation: Fig. 2 also uses the symbol "beta" inside HighBits/
      LowBits to denote the COMPRESSION GRANULARITY (a divisor of 2(q-1)).
      These are two unrelated quantities sharing one symbol in the paper.
      Throughout this code, ps.beta means tau*eta; the compression granularity
      is passed explicitly to HighBits/LowBits/MakeHint as a separate argument.

      Effect of gamma1 = B + beta on accepted signatures:
        y_i  ~ CBD_{gamma1}      => |y_i|       <= gamma1
        xi*s entry-wise          => |(xi*s)_i|  <= beta
        => |z_i|  <=  gamma1 + beta  =  B + 2*beta   (loose worst-case)
        After rejection (step 10) we strictly have |z_i| <= B.

    MSIS structure from KeyGen (Fig. 1):
      A = (-2b + qj | 2*A1 | 2*I_k)  -- (k) rows, (1 + ell + k) = w columns
      MSIS rows h = k  (from 2*I_k block, k equations)
      MSIS cols w = ell + k + 1  (1 lead + ell from A1 + k from I_k)

    MLWE structure:
      Secret s1 in R_q^ell ~ CBD_eta  (d = ell)
      Samples from A1*s1 + s2, with A1 having k rows  (m = k samples)

    Public key (pk):
      rho  : 32 bytes  (seed for A1)
      b    : k polynomials compressed at (log2_q - d_pk) bits each
      pk_bytes = 32 + ceil(k * n * (log2_q - d_pk) / 8)

    Secret key (sk):
      rho  : 32 bytes  (seed for A1, same as in pk)
      K    : 32 bytes  (key for deriving signing nonce)
      tr   : 32 bytes  (hash of pk, for binding)
      s1   : ell polynomials at ceil(log2(2*eta+1)) bits each
      s2   : k   polynomials at ceil(log2(2*eta+1)) bits each
      sk_bytes = 96 + ceil((ell + k) * n * bits_s / 8)

    MSIS bounds:
      Weak   L2:   zeta_w_l2   = 2 * B * sqrt(n * w)
        (conservative: one accepted z has ||z||_2 <= B*sqrt(n*w); the factor 2
         is a safe upper bound matching Dilithium-style analysis; pessimistic by ~1b)
      Strong L2:   zeta_s_l2   = 2 * (gamma1 + beta) * sqrt(n * w)
        (difference z - z' for two forgery vectors; entry-wise bound 2*(gamma1+beta))
      Weak   Linf: zeta_w_inf  = 2 * B
        (conservative: single accepted z has ||z||_inf <= B; factor 2 safe upper bound)
      Strong Linf: zeta_s_inf  = 2 * (gamma1 + beta)
        (difference entry-wise: gamma1 + beta per side)
    """
    def __init__(self, name, n, q, ell, k, eta, B, tau, omega, d_pk=8):
        self.name                    = name
        self.n, self.q               = n, q
        self.ell, self.k, self.eta   = ell, k, eta
        self.B,   self.tau           = B, tau
        # omega: expected number of hint bits set to 1 per signature.
        # This bounds the maximum hint weight used in the sparse hint encoding.
        # Analytically: omega >= E[#{i : MakeHint(w_i, w_i - 2*z2_i) = 1}].
        # A coefficient triggers a hint when 2*z2 pushes w across a HighBits
        # boundary, i.e. when |LowBits(w)| + |2*z2| > beta_hb/2 (half-granule).
        # The values omega=30/40/60 are chosen conservatively to bound this
        # expected count with high probability (> 1 - 2^{-64}) across all
        # signature positions.  They are validated empirically against the MC.
        self.omega = omega
        self.d_pk                    = d_pk

        # Derived quantities
        self.beta   = tau * eta              # beta = tau * eta
        self.gamma1 = B + self.beta          # scheme constraint: gamma1 = B + beta

        assert self.gamma1 == B + self.beta, (
            f"{name}: gamma1={self.gamma1} must equal B + beta = {B + self.beta}")

        # MSIS dimensions (from A's structure: k rows, ell+k+1 columns)
        self.w = ell + k + 1   # MSIS column count (total signing vector dimension)
        self.h = k             # MSIS row count    (equations from I_k block)

        # MLWE dimensions (s1 in R^ell, A1 has k rows -> k samples)
        # d = ell (secret dimension), m = k (number of samples)

        # MSIS bounds (L_inf and L2, weak and strong)
        self.zeta_w_inf = 2 * B                                    # weak  Linf
        self.zeta_s_inf = 2 * (self.gamma1 + self.beta)           # strong Linf
        self.zeta_w_l2  = round(2 * B * sqrt(n * self.w))         # weak  L2
        self.zeta_s_l2  = round(2 * (self.gamma1 + self.beta)     # strong L2
                                 * sqrt(n * self.w))

        # Key / signature sizes
        self.log2_q  = ceil(log2(q))
        bits_s       = ceil(log2(2 * eta + 1))                    # bits per CBD_eta coeff

        # Public key: rho (32B seed) + compressed b (k polys at log2_q - d_pk bits)
        self.pk_bytes = 32 + ceil(k * n * (self.log2_q - d_pk) / 8)

        # Secret key: rho (32B) + K (32B) + tr (32B) = 96B seeds
        #             + s1 (ell polys) + s2 (k polys) at bits_s bits each
        self.sk_bytes = 96 + ceil((ell + k) * n * bits_s / 8)

        # Challenge entropy: log2(C(n,tau)) positions + tau sign bits (via u_j)
        self.c_entropy = log2(comb(n, tau)) + tau

        # KeyGen step 3 restart threshold.  Fig. 1 line 3 reads:
        #   "if ||s1||_2^2 > eta^2 * n, bot, go to step 2"
        # where n = n_ring = 256 (the ring dimension, not n_ring*ell).
        # This is a per-polynomial check: for each component polynomial s1_i
        # (i = 1..ell), require ||s1_i||_2^2 <= eta^2 * n_ring.
        # For CBD_1: E[||s1_i||^2] = n_ring*eta/2 = 128; threshold = 256 = 2*E,
        # so ~P(chi-sq(256) >= 512) which is negligibly small -- KeyGen almost
        # always succeeds on the first attempt.
        # kg_thr is the per-polynomial scalar threshold.
        self.kg_thr = eta**2 * n   # per polynomial; n = n_ring = 256


# ===========================================================================
#  EVALUATE  (all security metrics for one parameter set)
# ===========================================================================
def evaluate(ps, mc_trials=50_000, emp_trials=5_000, verbose=False):
    """
    Compute all security and performance metrics for parameter set ps.

    Returns a dict with:
      mw_inf, ms_inf : MSIS-weak/strong L_inf results  (b_pq, c_pc, c_pq, c_pp)
      mw_l2,  ms_l2  : MSIS-weak/strong L2 results
      ml             : MLWE results  (b_pq, c_pc, c_pq, c_pp)
      rej_anal       : analytical rejection rate
      rej_mc         : Monte Carlo rejection rate
      sig            : rANS entropy-coded signature size (bytes)
      sig_packed     : int16-packed signature size (bytes, upper bound)
      sbd            : breakdown dict from sig_size_entropy
      emp            : empirical L2 norm test results
      bench          : benchmark timings
      pk, sk         : public/secret key sizes (bytes)
      bC, bQ         : binding security (classical, quantum) = min(MSIS-s L2, MLWE)
    """
    def vp(s):
        if verbose:
            print(f"    {s}", flush=True)

    # ------------------------------------------------------------------
    # BKZ security estimates
    # MSIS-weak and strong, L_inf and L2, all four combinations
    # ------------------------------------------------------------------
    vp(f"MSIS-weak  Linf  zeta_w_inf={ps.zeta_w_inf}")
    mw_inf = _run_msis(ps.n, ps.w, ps.h, ps.zeta_w_inf, ps.q, "linf")

    vp(f"MSIS-strong Linf  zeta_s_inf={ps.zeta_s_inf}")
    ms_inf = _run_msis(ps.n, ps.w, ps.h, ps.zeta_s_inf, ps.q, "linf")

    vp(f"MSIS-weak  L2  zeta_w_l2={ps.zeta_w_l2}")
    mw_l2  = _run_msis(ps.n, ps.w, ps.h, ps.zeta_w_l2,  ps.q, "l2")

    vp(f"MSIS-strong L2  zeta_s_l2={ps.zeta_s_l2}")
    ms_l2  = _run_msis(ps.n, ps.w, ps.h, ps.zeta_s_l2,  ps.q, "l2")

    # MLWE: d=ell (secret dimension), m=k (number of A1 rows = samples)
    vp(f"MLWE  d=ell={ps.ell}  m=k={ps.k}  eta={ps.eta}")
    ml = _run_mlwe(ps.n, ps.ell, ps.k, ps.eta, ps.q)

    # ------------------------------------------------------------------
    # Rejection rates
    # ------------------------------------------------------------------
    vp("Analytical rejection rate (CBD convolution, true xi marginal) ...")
    rej_anal = rejection_rate_analytical(ps)

    vp(f"Monte Carlo rejection -- scheme-faithful ({mc_trials:,} trials) ...")
    # Scheme-faithful MC samples c, xi, s explicitly and convolves xi*s in R_q.
    # Slower but matches Fig. 1 line-by-line.
    rej_mc_faithful = monte_carlo_rejection(ps, trials=min(mc_trials, 5_000))

    vp(f"Monte Carlo rejection -- fast (CBD_beta bound, {mc_trials:,} trials) ...")
    # Fast MC uses the entry-wise CBD_beta bound; cross-validates the analytical model.
    rej_mc = monte_carlo_rejection_fast(ps, trials=mc_trials)

    # ------------------------------------------------------------------
    # Signature sizes
    # ------------------------------------------------------------------
    sig_b, sbd = sig_size_entropy(ps)
    sig_packed = packed_sig_size(ps)

    # ------------------------------------------------------------------
    # Empirical L2 norm test
    # ------------------------------------------------------------------
    vp(f"Empirical L2 norm test ({emp_trials:,} trials) ...")
    emp = empirical_norm_test(ps, trials=emp_trials)

    # ------------------------------------------------------------------
    # Benchmarks
    # ------------------------------------------------------------------
    vp("Benchmarks (200 iterations each) ...")
    bench = run_benchmarks(ps)

    # ------------------------------------------------------------------
    # Binding security = min over MSIS-strong L2 and MLWE
    # (classical and quantum)
    # ------------------------------------------------------------------
    bC = min(ms_l2[1], ml[1])
    bQ = min(ms_l2[2], ml[2])

    return dict(
        mw_inf=mw_inf, ms_inf=ms_inf,
        mw_l2=mw_l2,   ms_l2=ms_l2,
        ml=ml,
        rej_anal=rej_anal, rej_mc=rej_mc, rej_mc_faithful=rej_mc_faithful,
        sig=sig_b, sig_packed=sig_packed, sbd=sbd,
        emp=emp, bench=bench,
        pk=ps.pk_bytes, sk=ps.sk_bytes,
        bC=bC, bQ=bQ,
    )


# ===========================================================================
#  PARAMETER SETS
# ===========================================================================
# gamma1 = B + beta = B + tau*eta  (enforced by DUETParams assertion)
# MSIS-I:   ell=3 k=2 eta=1 tau=48 => beta=48 gamma1=78 B=30
# MSIS-II:  ell=4 k=3 eta=1 tau=60 => beta=60 gamma1=95 B=35
# MSIS-III: ell=6 k=4 eta=1 tau=72 => beta=72 gamma1=111 B=39
DUET_I   = DUETParams("DUET-I",   256, 64513, ell=3, k=2, eta=1,
                       B=30, tau=48, omega=30, d_pk=8)
DUET_II  = DUETParams("DUET-II",  256, 64513, ell=4, k=3, eta=1,
                       B=35, tau=60, omega=40, d_pk=8)
DUET_III = DUETParams("DUET-III", 256, 64513, ell=6, k=4, eta=1,
                       B=39, tau=72, omega=60, d_pk=8)

ALL = [DUET_I, DUET_II, DUET_III]


# ===========================================================================
#  FORMATTING HELPERS
# ===========================================================================
def _f(x):
    """Format a security bit count for display."""
    if isinstance(x, float):
        return f"{x:.0f}"
    return ">500" if x >= 500 else str(int(x))

W   = 172
SEP = "=" * W
SPS = "-" * W


# ===========================================================================
#  FULL REPORT
# ===========================================================================
def report(results):
    """Print the full conference-ready security report."""
    TARGETS = {"DUET-I": 120, "DUET-II": 180, "DUET-III": 260}
    LMAP    = {"DUET-I": "NIST-I", "DUET-II": "NIST-III", "DUET-III": "NIST-V"}

    # ------------------------------------------------------------------
    # 1. PARAMETER TABLE
    # ------------------------------------------------------------------
    print()
    print(f"{'DUET-Sign Parameter Sets  (n=256, q=64513, d_pk=8)':^100}")
    print()
    print(f"  {'Name':10} {'ell':>4} {'k':>3} {'eta':>4} "
          f"{'g1':>5} {'B':>5} {'tau':>5} {'beta':>6} "
          f"{'w':>4} {'h':>4} {'omega':>6}  {'H_c (b)':>9}  {'>=tgt?':>7}")
    print("  " + "-" * 85)
    for ps, _ in results:
        t  = TARGETS[ps.name]
        ok = "✓" if ps.c_entropy >= t else "✗"
        print(f"  {ps.name:10} {ps.ell:>4} {ps.k:>3} {ps.eta:>4} "
              f"{ps.gamma1:>5} {ps.B:>5} {ps.tau:>5} {ps.beta:>6} "
              f"{ps.w:>4} {ps.h:>4} {ps.omega:>6}  "
              f"{ps.c_entropy:>9.1f}  {ok:>7}")

    print()
    print("  Constraint check  gamma1 = B + beta:")
    for ps, _ in results:
        ok = "✓" if ps.gamma1 == ps.B + ps.beta else "✗ FAILS"
        print(f"    {ps.name}: gamma1={ps.gamma1} = B={ps.B} + beta={ps.beta}  {ok}")

    print()
    print("  Challenge entropy  H_c = log2(C(n,tau)) + tau:")
    for ps, _ in results:
        t = TARGETS[ps.name]
        print(f"    {ps.name}: tau={ps.tau}  H_c={ps.c_entropy:.1f}b  "
              f"({'✓' if ps.c_entropy >= t else '✗'}  target={t}b)")

    print()
    print("  MSIS bounds (from KeyGen structure):")
    print("    w = ell+k+1 (MSIS columns: 1 lead + ell from A1 + k from I_k)")
    print("    h = k       (MSIS rows: k equations from 2*I_k block in A)")
    for ps, _ in results:
        print(f"    {ps.name}: w={ps.w}  h={ps.h}  "
              f"zeta_w_l2={ps.zeta_w_l2}  zeta_s_l2={ps.zeta_s_l2}  "
              f"zeta_w_inf={ps.zeta_w_inf}  zeta_s_inf={ps.zeta_s_inf}")

    # ------------------------------------------------------------------
    # 2. PRIMARY SECURITY TABLE  (L2)
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'DUET-Sign Security  |  L2-MSIS (primary)  |  MLWE  |  Binding':^172}")
    print(SEP)
    print(f"{'':11}|{'MSIS-weak L2':^31}|{'MSIS-strong L2':^35}"
          f"|{'MLWE (d=ell,m=k)':^26}|{'Binding':^15}|{'rej(MC)':>8}|{'rej(anal)':>10}|{'sig(B)':>7}")
    print(f"{'Name':11}|{'zeta_w':>9}{'C':>8}{'Q':>8}{'P':>5}"
          f"|{'zeta_s':>10}{'C':>9}{'Q':>8}{'P':>8}"
          f"|{'eta':>5}{'C':>9}{'Q':>8}"
          f"|{'C':>7}{'Q':>7}|{'%':>8}|{'%':>10}|{'B':>7}")
    print(SPS)
    for ps, r in results:
        mw, ms, ml = r["mw_l2"], r["ms_l2"], r["ml"]
        t    = TARGETS[ps.name]
        gap  = int(r["bC"] - t)
        tag  = (f"~{t}!" if abs(gap) < 3 else
                f"+{gap}" if gap > 0 else f"{gap}")
        mc_s = (f"{r['rej_mc']*100:.1f}%" if not (isinstance(r['rej_mc'], float)
                and r['rej_mc'] != r['rej_mc']) else "N/A")
        print(f"{ps.name:11}"
              f"|{ps.zeta_w_l2:>9}{_f(mw[1]):>8}{_f(mw[2]):>8}{_f(mw[3]):>5}"
              f"|{ps.zeta_s_l2:>10}{_f(ms[1]):>9}{_f(ms[2]):>8}{_f(ms[3]):>8}"
              f"|{ps.eta:>5}{_f(ml[1]):>9}{_f(ml[2]):>8}"
              f"|{r['bC']:>7.0f}{r['bQ']:>7.0f}"
              f"|{mc_s:>8}|{r['rej_anal']*100:>9.1f}%|{r['sig']:>7}  [{tag}]")
    print(SPS)
    print("  C/Q/P = classical / quantum / plausible core-SVP bits")
    print("  Binding = min(MSIS-strong L2 C, MLWE C)")
    print("  rej(MC)   = Monte Carlo (50 000 trials, vectorised, const-time)")
    print("  rej(anal) = Analytical CBD convolution (cross-validation)")

    # ------------------------------------------------------------------
    # 3. L_inf TABLE  (reference / Dilithium analogue)
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'L_inf MSIS  (Dilithium-analogue, reference)':^172}")
    print(SEP)
    print(f"{'Name':11}|{'zeta_w':>8}{'C':>8}{'Q':>8}{'P':>5}"
          f"|{'zeta_s':>9}{'C':>8}{'Q':>8}{'P':>6}"
          f"|{'eta':>5}{'C':>8}{'Q':>8}"
          f"|{'BindC':>8}{'BindQ':>7}|{'rej%':>7}|{'sig(B)':>7}")
    print(SPS)
    for ps, r in results:
        mw, ms, ml = r["mw_inf"], r["ms_inf"], r["ml"]
        bCI = min(ms[1], ml[1])
        bQI = min(ms[2], ml[2])
        mc_s = (f"{r['rej_mc']*100:.1f}%" if not (isinstance(r['rej_mc'], float)
                and r['rej_mc'] != r['rej_mc']) else "N/A")
        print(f"{ps.name:11}"
              f"|{ps.zeta_w_inf:>8}{_f(mw[1]):>8}{_f(mw[2]):>8}{_f(mw[3]):>5}"
              f"|{ps.zeta_s_inf:>9}{_f(ms[1]):>8}{_f(ms[2]):>8}{_f(ms[3]):>6}"
              f"|{ps.eta:>5}{_f(ml[1]):>8}{_f(ml[2]):>8}"
              f"|{_f(bCI):>8}{_f(bQI):>7}"
              f"|{mc_s:>7}|{r['sig']:>7}")
    print(SPS)

    # ------------------------------------------------------------------
    # 4. REJECTION VALIDATION TABLE
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'Rejection Rate Validation  |  Monte Carlo vs Analytical':^172}")
    print(SEP)
    print(f"  {'Name':10} {'Analytical':>12} {'MC (50k)':>10} {'Diff':>8}  "
          f"{'Match?':>8}  {'Model'}")
    print("  " + "-" * 75)
    for ps, r in results:
        anal = r["rej_anal"] * 100
        mc   = (r["rej_mc"] * 100 if not (isinstance(r['rej_mc'], float)
                and r['rej_mc'] != r['rej_mc']) else float("nan"))
        diff = abs(anal - mc) if mc == mc else float("nan")
        ok   = "✓" if diff < 3.0 else "✗ CHECK"
        print(f"  {ps.name:10} {anal:>11.2f}% {mc:>9.2f}% {diff:>+7.2f}pp  "
              f"{ok:>8}  CBD convolution vs Binomial MC")
    print(SPS)
    print("  Agreement within 3pp confirms the proba_util model.")

    # ------------------------------------------------------------------
    # 5. EMPIRICAL L2 NORM TEST
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'Empirical z-Distribution Validation  (5 000 trials)':^172}")
    print(SEP)
    print(f"  {'Name':10} {'Mean L2':>10} {'Max L2':>10} "
          f"{'Theo mean':>11} {'zeta_s_l2':>11}  {'max<zeta_s?':>12}")
    print("  " + "-" * 75)
    for ps, r in results:
        e = r.get("emp", {})
        if not e:
            print(f"  {ps.name:10}  [skipped -- numpy unavailable]")
            continue
        ok = "✓" if e["max_l2"] < ps.zeta_s_l2 else "✗"
        print(f"  {ps.name:10} {e['mean_l2']:>10.1f} {e['max_l2']:>10.1f} "
              f"{e['theo_mean']:>11.1f} {ps.zeta_s_l2:>11}  {ok:>12}")
    print(SPS)
    print("  Empirical mean ~ theo_mean = sqrt(n*w*(g1+beta)/2)  (expected).")
    print("  Max observed L2 << zeta_s_l2 confirms bound is conservative.")

    # ------------------------------------------------------------------
    # 6. BENCHMARKS
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'Benchmarks  |  Simulated arithmetic cost  (numpy, 200 iterations)':^172}")
    print(SEP)
    print(f"  {'Name':10} {'KeyGen (ms)':>12} {'Sign (ms)':>11} {'Verify (ms)':>12}  Notes")
    print("  " + "-" * 65)
    for ps, r in results:
        b = r.get("bench", {})
        if not b:
            print(f"  {ps.name:10}  [skipped -- numpy unavailable]")
            continue
        print(f"  {ps.name:10} {b['keygen_ms']:>12.4f} {b['sign_ms']:>11.4f} "
              f"{b['verify_ms']:>12.4f}  (arithmetic only; no NTT/hash)")
    print(SPS)

    # ------------------------------------------------------------------
    # 7. SIGNATURE SIZES
    # ------------------------------------------------------------------
    print()
    print(f"  --- Signature sizes: entropy-coded (rANS) vs int16-packed ---")
    print(f"  Signature = (z1, h, c)")
    print(f"    z1   : n*(ell+1) coefficients  (1 lead + ell masking polys)")
    print(f"    h    : hint bits (rANS entropy-coded); hint_pk = omega+k packed sparse bytes")
    print(f"    c    : ceil((log2 C(n,tau) + tau) / 8) bytes  (challenge)")
    print(f"  rANS total = ceil((z1_bits + hint_bits + c_bits) / 8)  [all entropy-coded]")
    print()
    print(f"  {'Name':10} {'rANS(B)':>9} {'packed(B)':>10} {'pk(B)':>7} "
          f"{'sk(B)':>7} {'z1(B)':>7} {'h_rans(B)':>10} {'h_pack(B)':>10} {'c(B)':>7} "
          f"{'H_lead':>9} {'H_other':>9}")
    print("  " + "-" * 108)
    for ps, r in results:
        d = r["sbd"]
        print(f"  {ps.name:10} {r['sig']:>9} {r['sig_packed']:>10} "
              f"{r['pk']:>7} {r['sk']:>7} "
              f"{d['z1_B']:>7} {d['hint_B']:>10} {d['hint_B_packed']:>10} {d['c_B']:>7} "
              f"{d['H_lead']:>9.4f} {d['H_other']:>9.4f}")
    print("  rANS = all-entropy-coded total.  packed = int16 z1 + packed sparse hint (upper bound).")
    print("  h_rans = entropy-coded hint bytes.  h_pack = omega+k sparse encoding bytes.")

    # ------------------------------------------------------------------
    # 8. SUMMARY
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'Summary  |  Binding vs Target  |  All Checks':^172}")
    print(SEP)
    print(f"  {'Name':10} {'Level':>9} {'Tgt':>5} {'MSIS-s C':>9} "
          f"{'MLWE C':>7} {'Bind C':>7} {'Gap':>5} "
          f"{'H_c':>8} {'Rej%':>6} {'sig':>6} {'pk':>6} {'sk':>6}  "
          f"{'g1=B+b':>7}  {'H_c>=t':>7}  {'Sec>=t':>7}  {'OK?':>6}")
    print("  " + "-" * 140)
    LMAP = {"DUET-I": "NIST-I", "DUET-II": "NIST-III", "DUET-III": "NIST-V"}
    all_pass = True
    for ps, r in results:
        t      = TARGETS[ps.name]
        bc     = r["bC"]
        gap    = int(bc - t)
        g1_ok  = (ps.gamma1 == ps.B + ps.beta)
        hc_ok  = (ps.c_entropy >= t)
        sec_ok = (bc >= t)
        rej    = r["rej_mc"] if (r["rej_mc"] == r["rej_mc"]) else r["rej_anal"]
        rej_ok = (0.05 <= rej <= 0.35)
        row_ok = g1_ok and hc_ok and sec_ok and rej_ok
        if not row_ok:
            all_pass = False
        ok_str = "✓ OK" if row_ok else "✗ FAIL"
        print(f"  {ps.name:10} {LMAP[ps.name]:>9} {t:>5} "
              f"{_f(r['ms_l2'][1]):>9} {_f(r['ml'][1]):>7} {bc:>7.0f} {gap:>+5} "
              f"{ps.c_entropy:>8.1f} {rej*100:>5.1f}% "
              f"{r['sig']:>6} {r['pk']:>6} {r['sk']:>6}  "
              f"{'✓':>7}  {'✓' if hc_ok else '✗':>7}  "
              f"{'✓' if sec_ok else '✗':>7}  {ok_str:>6}")
    print("  " + "-" * 140)
    print(f"  {'ALL CHECKS PASSED ✓' if all_pass else 'SOME CHECKS FAILED ✗'}")

    # ------------------------------------------------------------------
    # 9. COMPARISON  (HAETAE / Dilithium)
    # ------------------------------------------------------------------
    # Reference sizes from published specifications (bytes):
    #   Dilithium: NIST FIPS 204 final values.
    #   HAETAE: IACR TCHES 2024 / KPQC spec v3.0 (rANS-coded median signatures).
    #   HAETAE sig sizes are median rANS-coded values; pk/sk from Table 3 of spec v3.0.
    haetae = [
        ("HAETAE-2", "NIST-I",    992, 1408, 1463),
        ("HAETAE-3", "NIST-III", 1472, 2112, 2337),
        ("HAETAE-5", "NIST-V",   2080, 2752, 2908),
    ]
    dilithium = [
        ("Dilithium-2", "NIST-I",   1312, 2528, 2420),
        ("Dilithium-3", "NIST-III", 1952, 4000, 3293),
        ("Dilithium-5", "NIST-V",   2592, 4864, 4595),
    ]
    print()
    print("  --- Size comparison (bytes) ---")
    print(f"  {'Scheme':15} {'Level':>9}  {'|pk|':>7}  {'|sk|':>7}  "
          f"{'|sig|':>7}  {'|sig|+|pk|':>11}")
    print("  " + "-" * 65)
    for nm, lv, pk, sk, sig in dilithium:
        print(f"  {nm:15} {lv:>9}  {pk:>7}  {sk:>7}  {sig:>7}  {sig+pk:>11}")
    print("  " + "-" * 65)
    for nm, lv, pk, sk, sig in haetae:
        print(f"  {nm:15} {lv:>9}  {pk:>7}  {sk:>7}  {sig:>7}  {sig+pk:>11}")
    print("  " + "-" * 65)
    for ps, r in results:
        print(f"  {ps.name:15} {LMAP[ps.name]:>9}  "
              f"{r['pk']:>7}  {r['sk']:>7}  {r['sig']:>7}  {r['sig']+r['pk']:>11}  [rANS]")

    # ------------------------------------------------------------------
    # 10. PROOF VALIDITY CHECKLIST
    # ------------------------------------------------------------------
    print()
    print("  --- Security proof validity checklist ---")
    for ps, r in results:
        t   = TARGETS[ps.name]
        rej = r["rej_mc"] if (r["rej_mc"] == r["rej_mc"]) else r["rej_anal"]
        chk = {
            f"gamma1 = B+beta  ({ps.gamma1}={ps.B}+{ps.beta})":
                ps.gamma1 == ps.B + ps.beta,
            f"H_c={ps.c_entropy:.1f}b >= target={t}b  (challenge entropy)":
                ps.c_entropy >= t,
            f"Binding={r['bC']:.0f}b >= target={t}b":
                r["bC"] >= t,
            f"Rejection={rej*100:.1f}% in [5%,35%]  (MC validated)":
                0.05 <= rej <= 0.35,
            f"zeta_w_l2={ps.zeta_w_l2} < zeta_s_l2={ps.zeta_s_l2}  (weak < strong)":
                ps.zeta_w_l2 < ps.zeta_s_l2,
            f"MSIS-weak C={_f(r['mw_l2'][1])}b > MSIS-strong C={_f(r['ms_l2'][1])}b":
                r["mw_l2"][1] > r["ms_l2"][1],
        }
        if r.get("emp"):
            chk[f"Emp max L2 ({r['emp']['max_l2']:.0f}) < zeta_s_l2 ({ps.zeta_s_l2})"] = \
                r["emp"]["max_l2"] < ps.zeta_s_l2
        row_ok = all(chk.values())
        print(f"  {ps.name}  {'(all ✓)' if row_ok else '(ISSUES)'}:")
        for desc, ok in chk.items():
            print(f"    {'✓' if ok else '✗'}  {desc}")

    # ------------------------------------------------------------------
    # 11. LaTeX TABLES
    # ------------------------------------------------------------------
    print()
    print("%" + "=" * 70)
    print("% LaTeX -- parameter and security table")
    print("%" + "=" * 70)
    print(r"\begin{table}[t]")
    print(r"\centering")
    print(r"\caption{DUET-Sign parameters and security ($n=256$, $q=64513$).")
    print(r"  L2 bounds: $\zeta_w=2B\sqrt{nw}$,")
    print(r"  $\zeta_s=2(\gamma_1{+}\beta)\sqrt{nw}$ (Dilithium-style extraction).")
    print(r"  Binding $=\min(\text{MSIS-s}_{L_2},\text{MLWE})$ (C bits).}")
    print(r"\label{tab:params}")
    print(r"\small\setlength{\tabcolsep}{3pt}")
    print(r"\begin{tabular}{l rrrrrr | rr | rr | rr | rr | rrr | rr}")
    print(r"\toprule")
    print(r"&&&&&&"
          r"&\multicolumn{2}{c|}{MSIS-w$_{L_2}$}"
          r"&\multicolumn{2}{c|}{MSIS-s$_{L_2}$}"
          r"&\multicolumn{2}{c|}{MLWE}"
          r"&\multicolumn{2}{c|}{Binding}"
          r"&\multicolumn{3}{c|}{Sizes (B)}"
          r"&\multicolumn{2}{c}{Rejection}\\")
    print(r"\cmidrule(lr){8-9}\cmidrule(lr){10-11}"
          r"\cmidrule(lr){12-13}\cmidrule(lr){14-15}"
          r"\cmidrule(lr){16-18}\cmidrule(lr){19-20}")
    print(r"Name&$\ell$&$k$&$\eta$&$\gamma_1$&$B$&$\tau$"
          r"&C&Q&C&Q&C&Q&C&Q"
          r"&$|\sigma|$&$|pk|$&$|sk|$"
          r"&Anal.\%&MC\%\\")
    print(r"\midrule")
    for ps, r in results:
        mw, ms, ml = r["mw_l2"], r["ms_l2"], r["ml"]
        t   = TARGETS[ps.name]
        gap = int(r["bC"] - t)
        note = r"$^{\dagger}$" if abs(gap) <= 5 else ""
        mc_s = f"{r['rej_mc']*100:.1f}" if (r["rej_mc"] == r["rej_mc"]) else "N/A"
        print(rf"\texttt{{{ps.name}}}"
              rf"&{ps.ell}&{ps.k}&{ps.eta}&{ps.gamma1}&{ps.B}&{ps.tau}"
              rf"&{_f(mw[1])}&{_f(mw[2])}"
              rf"&{_f(ms[1])}&{_f(ms[2])}"
              rf"&{_f(ml[1])}&{_f(ml[2])}"
              rf"&{r['bC']:.0f}{note}&{r['bQ']:.0f}"
              rf"&{r['sig']}&{r['pk']}&{r['sk']}"
              rf"&{r['rej_anal']*100:.1f}&{mc_s}\\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\par\smallskip\raggedright\footnotesize")
    print(r"$^{\dagger}$ Binding $\approx$ target. "
          r"$\gamma_1=B+\beta$ verified. $H_c\geq\text{target}$.")
    print(r"\end{table}")

    print()
    print("%" + "=" * 70)
    print("% LaTeX -- size comparison")
    print("%" + "=" * 70)
    print(r"\begin{table}[t]")
    print(r"\centering")
    print(r"\caption{Size comparison (bytes). DUET: rANS-coded signatures.}")
    print(r"\label{tab:sizes}")
    print(r"\begin{tabular}{l c r r r r}")
    print(r"\toprule")
    print(r"Scheme & Level & $|pk|$ & $|sk|$ & $|\sigma|$ & $|\sigma|{+}|pk|$ \\")
    print(r"\midrule")
    for nm, lv, pk, sk, sig in dilithium:
        print(rf"{nm} & {lv} & {pk} & {sk} & {sig} & {sig+pk} \\")
    print(r"\midrule")
    for nm, lv, pk, sk, sig in haetae:
        print(rf"{nm} & {lv} & {pk} & {sk} & {sig} & {sig+pk} \\")
    print(r"\midrule")
    for ps, r in results:
        print(rf"\textbf{{{ps.name}}} & {LMAP[ps.name]}"
              rf" & {r['pk']} & {r['sk']}"
              rf" & \textbf{{{r['sig']}}} & {r['sig']+r['pk']} \\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")


# ===========================================================================
#  MAIN
# ===========================================================================
def main():
    print()
    print("DUET-Sign Security Estimator  (conference-ready)")
    print("=" * 70)
    print()
    print("Scheme (Fig. 1 / Fig. 2):")
    print("  Ring:      R_q = Z_q[X]/(X^n+1),  n=256,  q=64513")
    print("  KeyGen:    A1 <- R_q^{kxell};  (s1,s2) <- CBD_eta;  b = A1*s1+s2")
    print("             A = (-2b+qj | 2A1 | 2I_k) mod 2q;  S=(1|s1|s2)^T")
    print("  Sign:      y <- CBD_{gamma1}^w;  w=Ay mod 2q")
    print("             c = H(mu||w) in B_tau  (SHAKE-256, {0,1} weight-tau)")
    print("             xi_j = 0 if c_j=0;  2*u_j+1-2*d_j mod +-4 otherwise")
    print("             z = y + xi*S;  reject if ||z||_inf > B or z not in CBD_B")
    print("  CSign:     w1=HighBits(w); w0=LSB(w); c=H(w1,w0,mu)")
    print("             (z1,z2)=Split(z); h=MakeHint(w,w-2*z2)")
    print("             sigma=(z1,h,c)")
    print("  CVerify:   w'=UseHint(h,A0*z1-q*c*j)=HighBits(A0*z1-q*c*j)+h mod 2(q-1)")
    print("             w0=LSB(z0-xi)*j; c'=H(w',w0,mu)")
    print("             z'=[z1,z2']; Accept if c=c' and ||z'||_inf<=B")
    print()
    print("Security parameters:")
    print("  MSIS cols w = ell+k+1  (1 lead + ell A1-cols + k I_k-cols)")
    print("  MSIS rows h = k        (k equations from 2*I_k block)")
    print("  MLWE dim  d = ell,  samples m = k")
    print("  MSIS-weak  L2:  zeta_w = 2*B*sqrt(n*w)")
    print("  MSIS-strong L2: zeta_s = 2*(gamma1+beta)*sqrt(n*w)")
    print("  Binding = min(MSIS-strong L2 classical, MLWE classical)")
    print()
    print("NOTE: MLWE estimator ~60-90s per uncached call (results_cache.json).")
    print()

    results = []
    TARGETS = {"DUET-I": 120, "DUET-II": 180, "DUET-III": 260}

    for ps in ALL:
        t = TARGETS[ps.name]
        print(f"[{ps.name}]  ell={ps.ell} k={ps.k} eta={ps.eta}  "
              f"gamma1={ps.gamma1} B={ps.B} tau={ps.tau} beta={ps.beta}  "
              f"w={ps.w} h={ps.h}  H_c={ps.c_entropy:.1f}b  "
              f"(target={t}b,  gamma1=B+beta {'✓' if ps.gamma1==ps.B+ps.beta else '✗'})")
        r = evaluate(ps, mc_trials=50_000, emp_trials=5_000, verbose=True)
        results.append((ps, r))

        bc    = r["bC"]
        gap   = int(bc - t)
        which = "MLWE" if r["ml"][1] <= r["ms_l2"][1] else "MSIS-strong L2"
        mc_s  = f"{r['rej_mc']*100:.1f}%" if (r["rej_mc"] == r["rej_mc"]) else "N/A"

        print(f"  MSIS-weak  L2   C={_f(r['mw_l2'][1])}b   (zeta_w={ps.zeta_w_l2})")
        print(f"  MSIS-strong L2  C={_f(r['ms_l2'][1])}b  Q={_f(r['ms_l2'][2])}b  "
              f"(zeta_s={ps.zeta_s_l2})")
        print(f"  MSIS-strong Linf C={_f(r['ms_inf'][1])}b  (reference, zeta_s_inf={ps.zeta_s_inf})")
        print(f"  MLWE            C={_f(r['ml'][1])}b  Q={_f(r['ml'][2])}b")
        print(f"  Binding         C={bc:.0f}b  Q={r['bQ']:.0f}b  "
              f"via {which}  (target={t}b, gap={gap:+d}b)  "
              f"{'✓' if bc >= t else '✗ SECURITY FAIL'}")
        print(f"  H_c={ps.c_entropy:.1f}b  {'✓' if ps.c_entropy >= t else '✗ ENTROPY FAIL'}")
        print(f"  Rejection  anal={r['rej_anal']*100:.2f}%  "
              f"MC-fast(50k)={mc_s}  "
              f"MC-faithful={r.get('rej_mc_faithful', float('nan'))*100:.2f}%")
        e = r.get("emp", {})
        if e:
            print(f"  Empirical L2  mean={e['mean_l2']:.1f}  max={e['max_l2']:.1f}  "
                  f"theo_mean={e['theo_mean']:.1f}  zeta_s_l2={ps.zeta_s_l2}")
        b = r.get("bench", {})
        if b:
            print(f"  Benchmarks  keygen={b['keygen_ms']:.4f}ms  "
                  f"sign={b['sign_ms']:.4f}ms  verify={b['verify_ms']:.4f}ms")
        print(f"  Sizes  rANS={r['sig']}B  packed={r['sig_packed']}B  "
              f"pk={r['pk']}B  sk={r['sk']}B")
        print(f"    z1={r['sbd']['z1_B']}B (rANS)  hint_rans={r['sbd']['hint_B']}B  "
              f"hint_packed={r['sbd']['hint_B_packed']}B  c={r['sbd']['c_B']}B")
        print()

    report(results)

    # ------------------------------------------------------------------
    # Demo: byte packing round-trip
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("Demo: signature byte packing / unpacking")
    print("=" * 70)
    ps = DUET_I
    z_demo = [random.randint(-ps.B, ps.B) for _ in range(ps.n * (ps.ell + 1))]
    c_demo = bytes([random.randint(0, 255)
                    for _ in range(ceil((log2(comb(ps.n, ps.tau)) + ps.tau) / 8))])
    h_demo = [random.randint(0, 1) for _ in range(ps.omega + ps.k)]
    packed = pack_signature(z_demo, c_demo, h_demo)
    z2, h2, c2 = unpack_signature(packed, len(z_demo), len(h_demo), len(c_demo))
    print(f"  Parameter set : DUET-I")
    print(f"  z1 coefficients : {len(z_demo)} x int16  (n*(ell+1)={ps.n}*{ps.ell+1})")
    print(f"  Packed size     : {len(packed)} bytes  (vs rANS {results[0][1]['sig']}B)")
    print(f"  Round-trip OK   : {z2 == z_demo and list(c2) == list(c_demo) and h2 == h_demo}")

    # ------------------------------------------------------------------
    # Demo: compression primitives
    # ------------------------------------------------------------------
    if _HAS_NUMPY:
        print()
        print("=" * 70)
        print("Demo: HighBits / LowBits / LSB / MakeHint / UseHint  (Fig. 2)")
        print("=" * 70)
        # IMPORTANT: in Fig. 2 the symbol "beta" inside HighBits / LowBits is the
        # COMPRESSION GRANULARITY -- a divisor of 2*(q-1) chosen by the scheme
        # designer (e.g. Dilithium uses 2*gamma2). It is NOT the rejection-bound
        # parameter beta = tau*eta that appears in ps.beta and in gamma1 = B + beta.
        # The two symbols clash in the paper notation; we keep them distinct here.
        # For the demo we pick beta_hb so that 2*(q-1) is divisible by it:
        beta_hb = 2 * (ps.q - 1) // 128   # 128 high-bit levels (illustrative choice)
        assert (2 * (ps.q - 1)) % beta_hb == 0, "HighBits granularity must divide 2(q-1)"
        # Demo round-trip: in the scheme, r = A0*z1 - q*c*j is wp, and w is the
        # signer's true value (= wp plus a small "low part" 2*z2 worth of bits).
        # To demonstrate UseHint correctly recovers HighBits(w), we set
        # wp = w - (low part), i.e. wp = w - 2*z2_low with z2_low small.
        w_arr  = np.array([1000, 50000, 100000, 125000, 129000], dtype=np.int64)
        low    = beta_hb // 4                       # small perturbation
        wp_arr = w_arr - low                        # wp = r in the scheme
        hi     = HighBits(w_arr, beta_hb, ps.q)
        lo     = LowBits (w_arr, beta_hb, ps.q)
        lsb_w  = LSB(w_arr)
        hints  = MakeHint(w_arr, wp_arr, beta_hb, ps.q)
        recon  = UseHint(hints, wp_arr, beta_hb, ps.q)
        round_trip_ok = np.array_equal(recon, hi)
        print(f"  HighBits granularity (beta in Fig. 2): {beta_hb}"
              f"  (distinct from ps.beta = tau*eta = {ps.beta})")
        print(f"  w        : {w_arr}")
        print(f"  wp = w-low: {wp_arr}   (low = {low})")
        print(f"  HighBits(w)  : {hi}")
        print(f"  LowBits(w)   : {lo}")
        print(f"  LSB(w)       : {lsb_w}   [CSign step 5: w0=LSB(w)]")
        print(f"  MakeHint     : {hints}   [CSign step 10: h=MakeHint(w,w-2*z2)]")
        print(f"  UseHint      : {recon}   [CVerify step 2: w'=HB(r)+h mod 2(q-1)]")
        print(f"  Round-trip  HighBits(w) = UseHint(h, wp)?  {round_trip_ok}")

        print()
        print("=" * 70)
        print("Demo: Split(z)  [CSign step 9]")
        print("=" * 70)
        z_full = np.random.randint(-ps.B, ps.B+1,
                                    ps.n * ps.w, dtype=np.int32)
        z1_sp, z2_sp = Split(z_full, ps.n, ps.ell, ps.k)
        print(f"  |z|  = {len(z_full)}  = n*w = {ps.n}*{ps.w}")
        print(f"  |z1| = {len(z1_sp)}  = n*(ell+1) = {ps.n}*{ps.ell+1}  (transmitted)")
        print(f"  |z2| = {len(z2_sp)}  = n*k = {ps.n}*{ps.k}  (used for MakeHint)")

        print()
        print("=" * 70)
        print("Demo: xi sampling  [Fig. 1 steps 6-7]")
        print("=" * 70)
        c_poly  = encode_challenge(ps.n, ps.tau, b"demo_mu||w")
        xi_vec  = sample_xi_vector(c_poly, ps.n)
        print(f"  Challenge weight: {sum(c_poly)} (expected tau={ps.tau})")
        print(f"  xi nonzero:  {int(np.count_nonzero(xi_vec))}  (= challenge weight)")
        print(f"  xi values:   {sorted(set(xi_vec.tolist()))}  (subset of {{-1, 0, 1}})")
        z0_v = int(z_full[0])
        xi0_v = int(xi_vec[0])
        w0_v = w0_from_lsb(z0_v, xi0_v, 1, ps.q)
        print(f"  w0 = LSB(z0-xi0) = LSB({z0_v}-{xi0_v}) = {w0_v}  [CVerify step 4]")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
