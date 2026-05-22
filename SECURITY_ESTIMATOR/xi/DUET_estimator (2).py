
# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import sys, os, json, contextlib, io, struct, time, random, hashlib
from math import ceil, log2, comb, sqrt, exp

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
# pq-crystals security-estimates
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
#  3. COMPRESSION PRIMITIVES 
# ===========================================================================
def HighBits(w, beta, q):
    

    if _HAS_NUMPY and isinstance(w, np.ndarray):
        return (np.floor_divide(w.astype(np.int64), beta) * beta) % (2 * (q - 1))

    return ((int(w) // beta) * beta) % (2 * (q - 1))


def LowBits(w, beta, q):
    

    if _HAS_NUMPY and isinstance(w, np.ndarray):
        w64 = w.astype(np.int64)
        hi  = HighBits(w64, beta, q)
        return (w64 - hi) % (2 * q)

    hi = HighBits(w, beta, q)
    return (int(w) - hi) % (2 * q)


def LSB(w):
    
    if _HAS_NUMPY and isinstance(w, np.ndarray):
        return (w.astype(np.int64) % 2).astype(np.int32)
    return int(w) % 2

def MakeHint(w, wp, beta, q):
    

    modulus = 2 * (q - 1)

    hw  = HighBits(w,  beta, q)
    hwp = HighBits(wp, beta, q)

    if _HAS_NUMPY and isinstance(w, np.ndarray):
        return (hw.astype(np.int64) - hwp.astype(np.int64)) % modulus

    return (int(hw) - int(hwp)) % modulus


def UseHint(h, r, beta, q):
    

    modulus = 2 * (q - 1)

    if _HAS_NUMPY and isinstance(r, np.ndarray):

        hi = HighBits(r.astype(np.int64), beta, q)

        return (hi.astype(np.int64) + h.astype(np.int64)) % modulus

    hi = HighBits(int(r), beta, q)

    return (int(hi) + int(h)) % modulus

# ===========================================================================
#  4. XI SAMPLING  
# ===========================================================================
def sample_xi(c_coeff, u_j, d_j):
    
    if c_coeff == 0:
        return 0
    
    
    val = 2 * u_j + 1 - 2 * d_j
    
    if val == 3: return -1
    if val == -3: return 1
    return val

def sample_xi_vector(c, n, rng=None):
    
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
    
    return int(LSB(int(z0) - int(xi0)) * int(j))


# ===========================================================================
#  5. CHALLENGE ENCODING  
# ===========================================================================
def encode_challenge(n, tau, seed_bytes):
    
    c = [0] * n
    xof = hashlib.shake_256(seed_bytes)
    selected = set()
    
    block_size = 64
    block_idx  = 0
    pos_in_block = 0
    block = xof.digest(block_size)
    while len(selected) < tau:
        if pos_in_block + 2 > len(block):
            
            block_idx += 1
            block = xof.digest(block_size * (block_idx + 1))[block_size * block_idx:]
            pos_in_block = 0
        b0 = block[pos_in_block]
        b1 = block[pos_in_block + 1]
        pos_in_block += 2
       
        val = (b1 << 8) | b0
        if n & (n - 1) == 0 and n <= 65536:
            pos = val & (n - 1)            
        else:
            
            limit = (65536 // n) * n
            if val >= limit:
                continue                  
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
#  6. SPLIT(z)  
# ===========================================================================
def Split(z, n, ell, k):
    
    split_idx = n * (ell + 1)
    if _HAS_NUMPY and isinstance(z, np.ndarray):
        return z[:split_idx], z[split_idx:]
    return list(z[:split_idx]), list(z[split_idx:])


# ===========================================================================
#  7. SIGNATURE BYTE PACKING  
# ===========================================================================
def pack_signature(z1, c, h):
    
    data = bytearray()
    for x in z1:
        data += struct.pack('<h', int(x))       
    for x in h:
        data += struct.pack('B', int(x) & 0xFF)
    data += bytes(c)
    return bytes(data)


def unpack_signature(sig, z1_len, h_len, c_len):
    
    ptr = 0
    z1 = []

    for _ in range(z1_len):
        z1.append(struct.unpack('<h', sig[ptr:ptr + 2])[0])
        ptr += 2

    h = list(sig[ptr:ptr + h_len])
    ptr += h_len

    c = sig[ptr:ptr + c_len]

    return z1, h, c

# ===========================================================================
#  8. MONTE CARLO REJECTION RATE
# ===========================================================================
def monte_carlo_rejection(ps, trials=10_000, chunk=1_000, seed=42):
    
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

   
    idx_j = np.arange(n, dtype=np.int64)[:, None]   # (n,1)
    idx_l = np.arange(n, dtype=np.int64)[None, :]   # (1,n)
    shift = (idx_j - idx_l) % n                     # (n,n)
    sign  = np.where(idx_l > idx_j, -1, 1).astype(np.int32)  # (n,n)

    for start in range(0, trials, chunk):
        t_ = min(chunk, trials - start)

        c_batch = np.zeros((t_, n), dtype=np.int32)
        for r in range(t_):
            pos = rng.choice(n, size=tau, replace=False)
            c_batch[r, pos] = 1

    
        u_batch = (rng.binomial(2, 0.5, size=(t_, n)) - 1).astype(np.int32)
        d_batch = rng.integers(0, 2, size=(t_, n), dtype=np.int32)

    
        raw = 2*u_batch + 1 - 2*d_batch
        raw = np.where(raw > 2, raw - 4, np.where(raw < -2, raw + 4, raw))
        xi_batch = np.where(c_batch != 0, raw, np.int32(0)).astype(np.int32)

        
        s_batch = (rng.binomial(2*eta, 0.5, (t_, n_other, n)) - eta).astype(np.int32)

        
        y_batch = (rng.binomial(2*gamma1, 0.5, (t_, ell+k+1, n)) - gamma1).astype(np.int32)

       
        z_lead = y_batch[:, 0, :] + xi_batch    # (t_, n)

        
        s_gather = s_batch[:, :, shift]                     # (t_, n_other, n_j, n_l)
       
        weighted = s_gather * (xi_batch[:, None, None, :] * sign[None, None, :, :])
        conv     = weighted.sum(axis=-1).astype(np.int32)   # (t_, n_other, n)

        z_other  = y_batch[:, 1:, :] + conv                 # (t_, n_other, n)

       
        z_full = np.concatenate(
            (z_lead[:, None, :], z_other), axis=1
        ).reshape(t_, -1)                                   # (t_, (ell+k+1)*n)




def monte_carlo_rejection_fast(ps, trials=10_000, chunk=2_000, seed=42):
    
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
        beta = ps.tau * ps.eta

    r = beta / float(ps.B)

    M = np.exp(r + 0.5 * r * r)

    p_accept = 1.0 / M

    randvals = rng.random(t_)

    reject_mask = randvals > p_accept

    rejected += int(reject_mask.sum())
    return rejected / trials


# ===========================================================================
#  9. ANALYTICAL REJECTION RATE  
# ===========================================================================
def rejection_rate_analytical(ps):
    

    beta = ps.tau * ps.eta

    r = beta / float(ps.B)

    M = exp(r + 0.5 * r * r)

    return 1.0 - (1.0 / M)


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
    
    cbd_y   = build_centered_binomial_law(ps.gamma1)
    cbd_eta = build_centered_binomial_law(ps.eta)

    
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
    
    H_lead, H_other = coef_entropies(ps)
    z1_bits = ps.n * H_lead + ps.n * ps.ell * H_other   # z1 entropy (bits)

    
    p1 = ps.omega / (ps.n * ps.k)
    p0 = 1.0 - p1
    H_hint = 0.0
    if p0 > 0: H_hint -= p0 * log2(p0)
    if p1 > 0: H_hint -= p1 * log2(p1)
    hint_bits = ps.n * ps.k * H_hint

    c_bits = ceil(log2(comb(ps.n, ps.tau))) + ps.tau     
    total  = ceil((z1_bits + hint_bits + c_bits) / 8)

   
    return total, dict(
        z1_B         = ceil(z1_bits / 8),
        hint_B       = ceil(hint_bits / 8),   # entropy-coded hint size
        hint_B_packed= ps.omega + ps.k,        # packed sparse hint (upper bound)
        c_B          = ceil(c_bits / 8),
        H_lead       = H_lead,
        H_other      = H_other,
    )


def packed_sig_size(ps):
    
    z1_bytes   = ps.n * (ps.ell + 1) * 2          # ell+1 polys, int16
    hint_bytes = ps.omega + ps.k                   # sparse hint
    c_bytes    = ceil((log2(comb(ps.n, ps.tau)) + ps.tau) / 8)
    return z1_bytes + hint_bytes + c_bytes


# ===========================================================================
#  11. EMPIRICAL z-DISTRIBUTION VALIDATION
# ===========================================================================
def empirical_norm_test(ps, trials=3_000, seed=0):
    
    if not _HAS_NUMPY:
        print("    [empirical] numpy not available; skipping.")
        return {}
    rng     = np.random.default_rng(seed)
    dim     = ps.n * ps.w
    l2_vals = []
    attempts = 0
    while len(l2_vals) < trials:
        attempts += 1
       
        y   = (rng.binomial(2 * ps.gamma1, 0.5, dim) - ps.gamma1).astype(np.int32)
       
        xiS = (rng.binomial(2 * ps.beta, 0.5, dim) - ps.beta).astype(np.int32)
        z   = y + xiS
        
        if np.max(np.abs(z)) <= ps.B:
            l2_vals.append(float(np.linalg.norm(z.astype(np.float64), ord=2)))
        if attempts > trials * 200:   
            break
    if not l2_vals:
        return {}
    arr       = np.array(l2_vals)
    
    theo_mean = sqrt(ps.n * ps.w * (ps.gamma1 + ps.beta) / 2)
    return dict(
        mean_l2   = float(arr.mean()),
        max_l2    = float(arr.max()),
        theo_mean = theo_mean,
        zeta_s    = ps.zeta_s_inf,
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
    
    if not _HAS_NUMPY:
        return {}
    n, dim = ps.n, ps.n * ps.w
    rng = np.random.default_rng(0)

    def sim_keygen():
        
        s1 = (rng.binomial(2 * ps.eta, 0.5, (ps.ell, n)) - ps.eta).astype(np.int32)
        norms_sq = np.sum(s1 ** 2, axis=1)          # one norm per polynomial
        return bool(np.all(norms_sq <= ps.eta**2 * n))

    def sim_sign():
    
        y   = (rng.binomial(2 * ps.gamma1, 0.5, dim) - ps.gamma1).astype(np.int32)
        xiS = (rng.binomial(2 * ps.beta,   0.5, dim) - ps.beta  ).astype(np.int32)
        z   = y + xiS
        
        flag = int(np.bitwise_or.reduce((np.abs(z) > ps.B).astype(np.int32)))
        return flag == 0

    def sim_verify():
        
        z = (rng.binomial(2 * ps.gamma1, 0.5, dim) - ps.gamma1).astype(np.int32)
        return bool(np.max(np.abs(z)) <= ps.B)

    return dict(
        keygen_ms = benchmark(sim_keygen, 200) * 1000,
        sign_ms   = benchmark(sim_sign,   200) * 1000,
        verify_ms = benchmark(sim_verify, 200) * 1000,
    )


# ===========================================================================
#  13. CACHE  
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
    
    print(f"      [estimator] MSIS(n={n},w={w},h={h},zeta={zeta},q={q},"
          f"norm={norm}) ...", flush=True)
    
    r = MSIS_summarize_attacks(MSISParameterSet(n, w, h, zeta, q, norm=norm))
    c[key] = list(r)
    _save_cache(c)
    return tuple(r)

def _run_mlwe(n, d, m, eta, q):
    """Run MLWE security estimator (with cache). d=ell, m=k."""
    key = f"mlwe|n={n}|d={d}|m={m}|k={eta}|q={q}|distr=binomial"
    c   = _load_cache()
   
    print(f"      [estimator] MLWE(n={n},d={d},m={m},eta={eta},q={q}) "
          f"-- first run ~60-90s ...", flush=True)
    
    r = MLWE_summarize_attacks(MLWEParameterSet(n, d, m, eta, q, "binomial"))
    c[key] = list(r)
    _save_cache(c)
    return tuple(r)


# ===========================================================================
# 14.  PARAMETER CLASS
# ===========================================================================
class DUETParams:
    
    def __init__(self, name, n, q, ell, k, eta, B, tau, omega, d_pk=8):
        self.name                    = name
        self.n, self.q               = n, q
        self.ell, self.k, self.eta   = ell, k, eta
        self.B,   self.tau           = B, tau
        
        
        self.omega = omega
        self.d_pk  = d_pk

        # Derived quantities
        self.beta   = tau * eta              # beta = tau * eta
        self.gamma1 = B + self.beta          # scheme constraint: gamma1 = B + beta

        assert self.gamma1 == B + self.beta, (
            f"{name}: gamma1={self.gamma1} must equal B + beta = {B + self.beta}")

        # MSIS dimensions (from A's structure: k rows, ell+k+1 columns)
        self.w = ell + k + 1   # MSIS column count (total signing vector dimension)
        self.h = k             # MSIS row count    (equations from I_k block)

        

       
        self.zeta_w_inf = 2 * B                                    # weak  Linf
        self.zeta_s_inf = 2 * (self.gamma1 + self.beta)           # strong Linf
        
       
        self.log2_q  = ceil(log2(q))
        bits_s       = ceil(log2(2 * eta + 1))                    # bits per CBD_eta coeff

        
        self.pk_bytes = 32 + ceil(k * n * (self.log2_q - d_pk) / 8)

        
        self.sk_bytes = 96 + ceil((ell + k) * n * bits_s / 8)

        # Challenge entropy: log2(C(n,tau)) positions + tau sign bits (via u_j)
        self.c_entropy = log2(comb(n, tau)) + tau

        
        self.kg_thr = eta**2 * n   # per polynomial; n = n_ring = 256


# ==========================================================================
#  15. EVALUATE  (all security metrics for one parameter set)
# ===========================================================================
def evaluate(ps, mc_trials=10_000, emp_trials=3_000, verbose=False):
    
    def vp(s):
        if verbose:
            print(f"    {s}", flush=True)

   
    vp(f"MSIS-weak  Linf  zeta_w_inf={ps.zeta_w_inf}")
    mw_inf = _run_msis(ps.n, ps.w, ps.h, ps.zeta_w_inf, ps.q, "linf")

    vp(f"MSIS-strong Linf  zeta_s_inf={ps.zeta_s_inf}")
    ms_inf = _run_msis(ps.n, ps.w, ps.h, ps.zeta_s_inf, ps.q, "linf")

    vp(f"MLWE  d=ell={ps.ell}  m=k={ps.k}  eta={ps.eta}")
    ml = _run_mlwe(ps.n, ps.ell, ps.k, ps.eta, ps.q)

    # ------------------------------------------------------------------
    # Rejection rates
    # ------------------------------------------------------------------
    vp("Analytical rejection rate (CBD convolution, true xi marginal) ...")
    rej_anal = rejection_rate_analytical(ps)

    vp(f"Monte Carlo rejection -- scheme-faithful ({mc_trials:,} trials) ...")
   
    rej_mc_faithful = monte_carlo_rejection(ps, trials=min(mc_trials, 3_000))

    vp(f"Monte Carlo rejection -- fast (CBD_beta bound, {mc_trials:,} trials) ...")
    # Fast MC uses the entry-wise CBD_beta bound; cross-validates the analytical model.
    rej_mc = monte_carlo_rejection_fast(ps, trials=mc_trials)

    # ------------------------------------------------------------------
    # Signature sizes
    # ------------------------------------------------------------------
    sig_b, sbd = sig_size_entropy(ps)
    sig_packed = packed_sig_size(ps)

    # ------------------------------------------------------------------
    # Empirical  norm test
    # ------------------------------------------------------------------
    vp(f"Empirical L2 norm test ({emp_trials:,} trials) ...")
    emp = empirical_norm_test(ps, trials=emp_trials)

    # ------------------------------------------------------------------
    # Benchmarks
    # ------------------------------------------------------------------
    vp("Benchmarks (200 iterations each) ...")
    bench = run_benchmarks(ps)

    # ------------------------------------------------------------------
    # Binding security
    # (classical and quantum)
    # ------------------------------------------------------------------
    bC = min(ms_inf[1], ml[1])
    bQ = min(ms_inf[2], ml[2])

    return dict(
        mw_inf=mw_inf, ms_inf=ms_inf,
        #mw_l2=mw_l2,   ms_l2=ms_l2,
        ml=ml,
        rej_anal=rej_anal, rej_mc=rej_mc, rej_mc_faithful=rej_mc_faithful,
        sig=sig_b, sig_packed=sig_packed, sbd=sbd,
        emp=emp, bench=bench,
        pk=ps.pk_bytes, sk=ps.sk_bytes,
        bC=bC, bQ=bQ,
    )


# ===========================================================================
# 16.  PARAMETER SETS
# ===========================================================================
#
DUET_I   = DUETParams("DUET-I",   256, 64513, ell=3, k=2, eta=1,
                       B=300, tau=39, omega=5, d_pk=8)
DUET_II  = DUETParams("DUET-II",  256, 64513, ell=4, k=3, eta=1,
                       B=420, tau=49, omega=5, d_pk=8)
DUET_III = DUETParams("DUET-III", 256, 64513, ell=6, k=4, eta=1,
                       B=450, tau=60, omega=5, d_pk=8)

ALL = [DUET_I, DUET_II, DUET_III]


# ===========================================================================
# 17. FORMATTING HELPERS
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
# 18. FULL REPORT
# ===========================================================================
def report(results):
    """Print the full conference-ready security report."""
    TARGETS = {"DUET-I": 120, "DUET-II": 180, "DUET-III": 260}
    LMAP    = {"DUET-I": "NIST-I", "DUET-II": "NIST-III", "DUET-III": "NIST-V"}

    # ------------------------------------------------------------------
    #  PARAMETER TABLE
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
              f"zeta_w_inf={ps.zeta_w_inf}  zeta_s_inf={ps.zeta_s_inf}  "
              f"zeta_w_inf={ps.zeta_w_inf}  zeta_s_inf={ps.zeta_s_inf}")

    # ------------------------------------------------------------------
    #  PRIMARY SECURITY TABLE  (Linf)
    # ------------------------------------------------------------------
    print()
    print(SEP)
    print(f"{'DUET-Sign Security  |  Linf-MSIS (primary)  |  MLWE  |  Binding':^172}")
    print(SEP)
    print(f"{'':11}|{'MSIS-weak Linf':^31}|{'MSIS-strong Linf':^35}"
          f"|{'MLWE (d=ell,m=k)':^26}|{'Binding':^15}|{'rej(MC)':>8}|{'rej(anal)':>10}|{'sig(B)':>7}")
    print(f"{'Name':11}|{'zeta_w':>9}{'BKZ':>7}{'C':>8}{'Q':>8}{'P':>5}"
      f"|{'zeta_s':>10}{'BKZ':>7}{'C':>9}{'Q':>8}{'P':>8}"
      f"|{'eta':>5}{'BKZ':>7}{'C':>9}{'Q':>8}"
          f"|{'C':>7}{'Q':>7}|{'%':>8}|{'%':>10}|{'B':>7}")
    print(SPS)
    for ps, r in results:
        mw, ms, ml = r["mw_inf"], r["ms_inf"], r["ml"]
        t    = TARGETS[ps.name]
        gap  = int(r["bC"] - t)
        tag  = (f"~{t}!" if abs(gap) < 3 else
                f"+{gap}" if gap > 0 else f"{gap}")
        mc_s = f"{r['rej_mc']*100:.2f}%"
        print(f"{ps.name:11}"
      f"|{ps.zeta_w_inf:>9}{mw[0]:>7.0f}{_f(mw[1]):>8}{_f(mw[2]):>8}{_f(mw[3]):>5}"
      f"|{ps.zeta_s_inf:>10}{ms[0]:>7.0f}{_f(ms[1]):>9}{_f(ms[2]):>8}{_f(ms[3]):>8}"
      f"|{ps.eta:>5}{ml[0]:>7.0f}{_f(ml[1]):>9}{_f(ml[2]):>8}"
              f"|{r['bC']:>7.0f}{r['bQ']:>7.0f}"
              f"|{mc_s:>8}|{r['rej_anal']*100:>9.1f}%|{r['sig']:>7}  [{tag}]")
    print(SPS)
    print("  C/Q/P = classical / quantum / plausible core-SVP bits")
  
    print("  rej(MC)   = Monte Carlo (30 000 trials, vectorised, const-time)")
    print("  rej(anal) = Analytical CBD convolution ")

    # ------------------------------------------------------------------
    #  -------------------L_inf TABLE -------------------------------
    print()
    print(SEP)
    print(f"{'L_inf MSIS  ':^172}")
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
        mc_s = f"{r['rej_mc']*100:.2f}%"
        print(f"{ps.name:11}"
              f"|{ps.zeta_w_inf:>8}{_f(mw[1]):>8}{_f(mw[2]):>8}{_f(mw[3]):>5}"
              f"|{ps.zeta_s_inf:>9}{_f(ms[1]):>8}{_f(ms[2]):>8}{_f(ms[3]):>6}"
              f"|{ps.eta:>5}{_f(ml[1]):>8}{_f(ml[2]):>8}"
              f"|{_f(bCI):>8}{_f(bQI):>7}"
              f"|{mc_s:>7}|{r['sig']:>7}")
    print(SPS)

    
    

    # ------------------------------------------------------------------
    #  SIGNATURE SIZES
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
    print("  rANS = all-entropy-coded total.  packed = int16 z1 + packed sparse hint .")
    print("  h_rans = entropy-coded hint bytes.  h_pack = omega+k sparse encoding bytes.")

    

# ===========================================================================
#  MAIN
# ===========================================================================
def main():
    print()
    print("DUET-Sign Security Estimator  ")
    print("=" * 70)
    print()
    print("Scheme :")
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
   
    print("  MSIS-weak  Linf:  zeta_w = 2*B")
    print("  MSIS-strong Linf: zeta_s = 2*(gamma1+beta)")
    
    print()
   
    print()

    results = []
    TARGETS = {"DUET-I": 120, "DUET-II": 180, "DUET-III": 260}

    for ps in ALL:
        t = TARGETS[ps.name]
        print(f"[{ps.name}]  ell={ps.ell} k={ps.k} eta={ps.eta}  "
              f"gamma1={ps.gamma1} B={ps.B} tau={ps.tau} beta={ps.beta}  "
              f"w={ps.w} h={ps.h}  H_c={ps.c_entropy:.1f}b  "
              f"(target={t}b,  gamma1=B+beta {'✓' if ps.gamma1==ps.B+ps.beta else '✗'})")
        r = evaluate(ps, mc_trials=10_000, emp_trials=3_000, verbose=True)
        results.append((ps, r))

        bc    = r["bC"]
        gap   = int(bc - t)
        which = "MLWE" if r["ml"][1] <= r["ms_inf"][1] else "MSIS-strong L2"
        mc_s  = f"{r['rej_mc']*100:.1f}%" if (r["rej_mc"] == r["rej_mc"]) else "N/A"

        print(f"  MSIS-strong Linf C={_f(r['ms_inf'][1])}b  (reference, zeta_s_inf={ps.zeta_s_inf})")
        print(f"  MLWE            BKZ={r['ml'][0]:.0f}  "
              f"C={_f(r['ml'][1])}b  Q={_f(r['ml'][2])}b")
        
        print(f"  H_c={ps.c_entropy:.1f}b  {'✓' if ps.c_entropy >= t else '✗ ENTROPY FAIL'}")
        print(f"  Rejection  anal={r['rej_anal']*100:.2f}%  "
              f"MC-fast(10k)={mc_s}  "
              f"MC-faithful={r.get('rej_mc_faithful', float('nan'))*100:.2f}%")
        e = r.get("emp", {})
        if e:
            print(f"  Empirical L2  mean={e['mean_l2']:.1f}  max={e['max_l2']:.1f}  "
                  f"theo_mean={e['theo_mean']:.1f} zeta_s_inf={ps.zeta_s_inf}")
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

    
    
    print("Done.")


if __name__ == "__main__":
    main()
