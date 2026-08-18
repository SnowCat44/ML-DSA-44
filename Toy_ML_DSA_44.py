import hashlib
import os

Q      = 8380417
D      = 13
GAMMA1 = (1 << 17)
GAMMA2 = ((Q-1)//88)

N   = 8 # 다항식 차수 (표준: 256)
K   = 4
L   = 4
TAU = 5
ETA = 2

# NTT용 zeta 값, ZETAS[0] = 0은 미사용 (NTT에서 m = 1부터 접근)
ZETAS = [0, 4808194, 4614810, 4618904, 2883726, 5178923, 5178987, 3145678]


# ================================================================================
# ====================== Hash Function: SHAKE-256/SHAKE-128 ======================
# ================================================================================

# ========================= 해시 함수: SHAKE-256/SHAKE-128 =========================
def H(data, l):
    ctx = hashlib.shake_256()
    ctx.update(data)

    return ctx.digest(8 * l)

def G(data: bytes, length: int) -> bytes:
    ctx = hashlib.shake_128()
    ctx.update(data)

    return ctx.digest(length)

def H_Init():
    return hashlib.shake_256(), 0

def G_Init():
    return hashlib.shake_128(), 0

def H_Absorb(ctx_state, data: bytes):
    ctx, offset = ctx_state
    ctx.update(data)

    return ctx, offset

def G_Absorb(ctx_state, data: bytes):
    ctx, offset = ctx_state
    ctx.update(data)

    return ctx, offset

def H_Squeeze(ctx_state, length: int) -> tuple:
    ctx, offset = ctx_state
    out = ctx.digest(offset + length)[offset:]

    return (ctx, offset + length), out

def G_Squeeze(ctx_state, length: int) -> tuple:
    ctx, offset = ctx_state
    out = ctx.digest(offset + length)[offset:]

    return (ctx, offset + length), out
# ================================================================================



# ================================================================================
# ============================== Utility Functions ===============================
# ================================================================================

# =================================== 덧셈 연산 ====================================
def Add(a, b): # 계수
    return (a + b) % Q

def Add_Poly(f, g): # 다항식
    return [Add(f[i], g[i]) for i in range(N)]

def Add_PolyVec(v, w): # 다항식 벡터
    return [Add_Poly(v[i], w[i]) for i in range(len(v))]
# ================================================================================

# =================================== 뺄셈 연산 ====================================
def Sub(a, b): # 계수
    return (a - b) % Q

def Sub_Poly(f, g): # 다항식
    return [Sub(f[i], g[i]) for i in range(N)]

def Sub_PolyVec(v, w): # 다항식 벡터
    return [Sub_Poly(v[i], w[i]) for i in range(len(v))]
# ================================================================================

# =================================== 곱셈 연산 ====================================
def Mul(a, b): # 계수
    return (a * b) % Q

def Mul_Point(f, g): # point-wise
    return [Mul(f[i], g[i]) for i in range(N)]
# ================================================================================



# ================================================================================
# =============================== External Functions =============================
# ================================================================================

def ML_DSA_KeyGen():
    try:
        nu = os.urandom(4)  # 4 bytes 랜덤 시드 생성
    except Exception:
        return None  # 시드 생성 실패 → Fail

    if not nu:
        return None  # 시드가 NULL → Fail

    return ML_DSA_KeyGen_internal(nu)

def ML_DSA_Sign(sk, M):
    sigma = ML_DSA_Sign_internal(sk, M)
    if sigma is None:
        return None  # 서명 없음 → Fail

    return sigma

def ML_DSA_Verify(pk, M, sigma):
    if sigma is None:
        return False  # 서명 없음 → Fail

    return ML_DSA_Verify_internal(pk, M, sigma)

def HashML_DSA_Sign(sk, M):
    sigma = ML_DSA_Sign_internal(sk, M)
    if sigma is None:
        return None  # 서명 실패 → Fail

    c_tilde, z, h = sigma
    h_bytes = bytes(bit for poly_h in h for bit in poly_h)
    return c_tilde + b''.join(coef.to_bytes(4, 'little') for poly in z for coef in poly) + h_bytes

def HashML_DSA_Verify(pk, M, signature_bytes):

    if signature_bytes is None or len(signature_bytes) < 8 + L * N * 4 + K * N:
        return False

    c_tilde = signature_bytes[:8]
    z_bytes = signature_bytes[8:8 + L * N * 4]
    h_bytes = signature_bytes[8 + L * N * 4:8 + L * N * 4 + K * N]

    z = [
        [int.from_bytes(z_bytes[(i * N + j) * 4:(i * N + j + 1) * 4], 'little') for j in range(N)]
        for i in range(L)
    ]
    h = [list(h_bytes[i * N:(i + 1) * N]) for i in range(K)]

    return ML_DSA_Verify_internal(pk, M, (c_tilde, z, h))
# ================================================================================



# ================================================================================
# =============================== Internal Functions =============================
# ================================================================================

def ML_DSA_KeyGen_internal(nu):
    seed = H(nu + K.to_bytes(1, 'little') + L.to_bytes(1, 'little'), 2)  # 8*2 = 16 bytes
    rho   = seed[0:4]   # ρ  (4 bytes) → ExpandA 입력
    rho_p = seed[4:12]  # ρ' (8 bytes) → ExpandS 입력
    k_key = seed[12:16] # k  (4 bytes) → 서명 시 사용

    A_hat  = ExpandA(rho)
    s1, s2 = ExpandS(rho_p)

    t_hat  = MatrixMulNTT_PolyVec(A_hat, NTT_PolyVec(s1))
    t      = Add_PolyVec(InvNTT_PolyVec(t_hat), s2)
    t1, t0 = Power2Round_PolyVec(t)

    pk = (rho, t1)
    pk_bytes = rho + b''.join(c.to_bytes(4, 'little') for poly in t1 for c in poly)
    tr = H(pk_bytes, 1)
    sk = (rho, k_key, tr, s1, s2, t0)

    return pk, sk

def ML_DSA_Sign_internal(sk, M_):
    rho, k_key, tr, s1, s2, t0 = sk
    s1_hat = NTT_PolyVec(s1)
    s2_hat = NTT_PolyVec(s2)
    t0_hat = NTT_PolyVec(t0)
    A_hat  = ExpandA(rho)
    mu = H(tr + M_, 1)
    rho__ = H(k_key + mu, 1)
    kappa = 0
    z, h = None, None

    while z is None or h is None:
        y = ExpandMask(rho__, kappa)
        w = InvNTT_PolyVec(MatrixMulNTT_PolyVec(A_hat, NTT_PolyVec(y)))
        w1 = HighBits_PolyVec(w)

        c_tilde = H(mu + b''.join(c.to_bytes(4, 'little') for poly in w1 for c in poly), 1)
        c = SampleInBall(c_tilde)
        c_hat = NTT_Poly(c)

        cs1 = InvNTT_PolyVec(MulNTT_PolyVec([c_hat] * L, s1_hat))
        cs2 = InvNTT_PolyVec(MulNTT_PolyVec([c_hat] * K, s2_hat))
        z = Add_PolyVec(cs1, y)
        r0 = LowBits_PolyVec(Sub_PolyVec(w, cs2))

        if any(abs(v if v <= Q // 2 else v - Q) >= GAMMA1 - TAU * ETA for poly in z for v in poly) or \
           any(abs(v if v <= Q // 2 else v - Q) >= GAMMA2 - TAU * ETA for poly in r0 for v in poly):
            z, h = None, None
        else:
            ct0 = InvNTT_PolyVec(MulNTT_PolyVec([c_hat] * K, t0_hat))
            neg_ct0 = Sub_PolyVec([[0] * N for _ in range(K)], ct0)

            h = MakeHint_PolyVec(neg_ct0, Add_PolyVec(Sub_PolyVec(w, cs2), ct0))
        kappa += L

    sigma = (c_tilde, z, h)

    return sigma

def ML_DSA_Verify_internal(pk, M_, sigma):
    rho, t1 = pk
    c_tilde, z, h = sigma

    if h is None:
        return False

    A_hat  = ExpandA(rho)
    tr = H(rho + b''.join(c.to_bytes(4, 'little') for poly in t1 for c in poly), 1)
    mu = H(tr + M_, 1)

    c = SampleInBall(c_tilde)
    t1_2d = [[coef * (1 << D) % Q for coef in poly] for poly in t1]
    w_Approx = InvNTT_PolyVec(Sub_PolyVec(MatrixMulNTT_PolyVec(A_hat, NTT_PolyVec(z)), MulNTT_PolyVec([NTT_Poly(c)] * K, NTT_PolyVec(t1_2d))))
    w_ = UseHint_PolyVec(h, w_Approx)

    c_tilde_ = H(mu + b''.join(c.to_bytes(4, 'little') for poly in w_ for c in poly), 1)

    if any(abs(v if v <= Q // 2 else v - Q) >= GAMMA1 - TAU * ETA for poly in z for v in poly):
        return False
    return c_tilde == c_tilde_

# ================================================================================



# ================================================================================
# ============================== Auxiliary Functions =============================
# ================================================================================

def CoeffFromThreeBytes(b0, b1, b2):
    b2_ = b2

    if b2_ > 127:
        b2_ -= 128
    z = (2 ** 16) * b2_ + (2 ** 8) * b1 + b0
    if z < Q:
        return z
    else:
        return None

def CoeffFromHalfByte(b):
    if ETA == 2 and b < 15:
        return 2 - (b % 5)
    elif ETA == 4 and b < 9:
        return 4 - b
    else:
        return None

def BitUnpack(v, a, b):
    c = (a + b).bit_length()
    bits = []
    for byte in v:
        for i in range(8):
            bits.append((byte >> i) & 1)
    w = [0] * N
    for i in range(N):
        z = 0
        for j in range(c):
            z |= bits[i * c + j] << j
        w[i] = b - z

    return w

# ============================== Sampling Functions ==============================
def SampleInBall(rho):
    c = [0] * N
    ctx = H_Init()
    ctx = H_Absorb(ctx, rho)
    ctx, s = H_Squeeze(ctx, 8)

    h = []
    for byte in s:
        for i in range(8):
            h.append((byte >> i) & 1)

    for i in range(N - TAU, N):
        ctx, j_bytes = H_Squeeze(ctx, 1)
        j = j_bytes[0]
        while j > i:
            ctx, j_bytes = H_Squeeze(ctx, 1)
            j = j_bytes[0]
        c[i] = c[j]
        c[j] = (-1) ** h[i + TAU - N]

    return c

def RejNTTPoly(rho):
    j = 0
    a_hat = [0] * N
    ctx = G_Init()
    ctx = G_Absorb(ctx, rho)

    while j < N:
        (ctx, s) = G_Squeeze(ctx, 3)
        a_hat[j] = CoeffFromThreeBytes(s[0], s[1], s[2])
        if a_hat[j] is not None:
            j += 1

    return a_hat

def RejBoundedPoly(rho):
    j = 0
    a = [0] * N
    ctx = H_Init()
    ctx = H_Absorb(ctx, rho)

    while j < N:
        ctx, z_bytes = H_Squeeze(ctx, 1)
        z = z_bytes[0]
        z0 = CoeffFromHalfByte(z % 16)
        z1 = CoeffFromHalfByte(z >> 4)
        if z0 is not None:
            a[j] = z0
            j += 1
        if z1 is not None and j < N:
            a[j] = z1
            j += 1

    return a
# ================================================================================

# ============================== Expand Functions ================================
def ExpandA(rho):
    A_hat = [[None] * L for _ in range(K)]

    for r in range(K):
        for s in range(L):
            rho_ = rho + s.to_bytes(1, 'little') + r.to_bytes(1, 'little')
            A_hat[r][s] = RejNTTPoly(rho_)

    return A_hat

def ExpandS(rho):
    s1 = [[None] * N for _ in range(L)]
    s2 = [[None] * N for _ in range(K)]

    for r in range(L):
        s1[r] = RejBoundedPoly(rho + r.to_bytes(2, 'little'))
    for r in range(K):
        s2[r] = RejBoundedPoly(rho + (r + L).to_bytes(2, 'little'))

    return s1, s2

def ExpandMask(rho, kappa):
    c = 1 + (GAMMA1 - 1).bit_length()
    y = [[None] * N for _ in range(L)]

    for r in range(L):
        rho_ = rho + (kappa + r).to_bytes(2, 'little')
        v = H(rho_, c << 5)
        y[r] = BitUnpack(v, GAMMA1 - 1, GAMMA1)

    return y
# ================================================================================

# ================================= Power2Round ==================================
def Power2Round(r): # 계수
    r_plus = r % Q
    r0 = r_plus % (1 << D)

    if r0 > (1 << (D - 1)):
        r0 -= (1 << D)
    r1 = (r_plus - r0) >> D

    return r1, r0

def Power2Round_Poly(r): # 다항식
    r1 = [0] * N
    r0 = [0] * N

    for i in range(N):
        r1[i], r0[i] = Power2Round(r[i])

    return r1, r0

def Power2Round_PolyVec(r): # 다항식 벡터
    r1 = []
    r0 = []

    for i in range(len(r)):
        p1, p0 = Power2Round_Poly(r[i])
        r1.append(p1)
        r0.append(p0)

    return r1, r0
# ================================================================================

# =============================== Decompose / Bits ===============================
def Decompose(r): # 계수
    r_plus = r % Q
    r0 = r_plus % (2 * GAMMA2)

    if r0 > GAMMA2:
        r0 -= (2 * GAMMA2)
    if (r_plus - r0) == Q - 1:
        return 0, r0 - 1
    r1 = (r_plus - r0) // (2 * GAMMA2)

    return r1, r0

def Decompose_Poly(r): # 다항식
    r1 = [0] * N
    r0 = [0] * N

    for i in range(N):
        r1[i], r0[i] = Decompose(r[i])

    return r1, r0

def HighBits(r): # 계수
    return Decompose(r)[0]

def LowBits(r): # 계수
    return Decompose(r)[1]

def HighBits_Poly(r): # 다항식
    return Decompose_Poly(r)[0]

def LowBits_Poly(r): # 다항식
    return Decompose_Poly(r)[1]

def HighBits_PolyVec(r): # 다항식 벡터
    return [HighBits_Poly(r[i]) for i in range(len(r))]

def LowBits_PolyVec(r): # 다항식 벡터
    return [LowBits_Poly(r[i]) for i in range(len(r))]
# ================================================================================

# ===================================== Hint =====================================
def MakeHint(z, r):
    return HighBits(r) != HighBits(r + z)

def MakeHint_Poly(z, r):
    return [MakeHint(z[i], r[i]) for i in range(N)]

def MakeHint_PolyVec(z, r):
    return [MakeHint_Poly(z[i], r[i]) for i in range(len(z))]

def UseHint(h, r):
    m = (Q - 1) // (2 * GAMMA2)
    r1, r0 = Decompose(r)

    if h != 1:
        return r1
    if r0 > 0:
        return (r1 + 1) % m
    return (r1 - 1) % m

def UseHint_Poly(h, r):
    return [UseHint(h[i], r[i]) for i in range(N)]

def UseHint_PolyVec(h, r):
    return [UseHint_Poly(h[i], r[i]) for i in range(len(h))]
# ================================================================================

# ===================================  NTT 변환 ===================================
def NTT_Poly(a): # 다항식 변환
    a_hat = list(a)

    m = 0
    length = N >> 1

    while length >= 1:
        start = 0
        while start < N:
            m += 1
            z = ZETAS[m]
            for j in range(start, start + length):
                t = Mul(z, a_hat[j + length])
                a_hat[j + length] = Sub(a_hat[j], t)
                a_hat[j] = Add(a_hat[j], t)
            start += 2 * length
        length >>= 1

    return a_hat

def InvNTT_Poly(a_hat): # 다항식 변환
    a = list(a_hat)

    m = N
    length = 1

    while length < N:
        start = 0
        while start < N:
            m -= 1
            zeta = -ZETAS[m]
            for j in range(start, start + length):
                t = a[j]
                a[j] = Add(t, a[j + length])
                a[j + length] = Sub(t, a[j + length])
                a[j + length] = Mul(zeta, a[j + length])
            start += 2 * length
        length <<= 1

    f = pow(N, -1, Q)
    for j in range(N):
        a[j] = Mul(f, a[j])

    return a

def NTT_PolyVec(a): # 다항식 벡터 변환
    return [NTT_Poly(a[i]) for i in range(len(a))]

def InvNTT_PolyVec(a_hat): # 다항식 벡터 변환
    return [InvNTT_Poly(a_hat[i]) for i in range(len(a_hat))]
# ================================================================================

# =================================== NTT 연산 ====================================
def AddNTT_Poly(f_hat, g_hat): # 다항식
    return [Add(f_hat[i], g_hat[i]) for i in range(N)]

def AddNTT_PolyVec(v_hat, w_hat): # 다항식 벡터
    return [AddNTT_Poly(v_hat[i], w_hat[i]) for i in range(len(v_hat))]

def MulNTT_Poly(f_hat, g_hat): # Point-wise
    return Mul_Point(f_hat, g_hat)

def MulNTT_PolyVec(v_hat, w_hat): # 다항식 벡터
    return [MulNTT_Poly(v_hat[i], w_hat[i]) for i in range(len(v_hat))]

def ScalarMulNTT_Poly(f_hat, scalar): # 다항식
    return [MulNTT_Poly(scalar, f_hat[i]) for i in range(N)]

def MatrixMulNTT_PolyVec(M_hat, v_hat): # 행렬(K x L) x 벡터(L)
    w = []

    for i in range(K):
        acc = [0] * N
        for j in range(L):
            acc = AddNTT_Poly(acc, MulNTT_Poly(M_hat[i][j], v_hat[j]))
        w.append(acc)

    return w
# ================================================================================
