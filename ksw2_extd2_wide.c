/* ksw2_extd2_wide.c — AVX2 / AVX-512 widening of minibwa's ksw_extd2_sse
 * (two-piece / dual-gap affine extension, Suzuki-Kasahara). This is minibwa's
 * PRIMARY base-alignment kernel for gap filling between seeds (align.c), and the
 * dominant cost for long reads (~52-66% of HiFi/ONT compute).
 *
 * Only the path minibwa's aligner actually uses for gap filling is widened:
 *   cigar (not SCORE_ONLY) + APPROX_MAX + left-alignment (not RIGHT) +
 *   exact match/mismatch scoring (not GENERIC_SC) + not EXTZ_ONLY.
 * Everything else dispatches to the unmodified original (extd2_ref_impl).
 *
 * Correctness is bit-identical to the 128-bit kernel by construction: the band
 * is rounded to 16 bytes exactly as the original, and that 16-aligned byte range
 * is processed in W-byte chunks + a 16-byte SSE tail, with the three
 * anti-diagonal carries (x1, x21, v1) threaded as scalars. The direction byte
 * `d` is computed identically and stored W-at-a-time at the same byte offsets in
 * the `p` array, so ksw_backtrack (and n_col_, off[]) are untouched. See
 * ksw2_wide_simd.h for the width abstraction and the cross-lane byte shift.
 */
#include <string.h>
#include <stdint.h>
#include "ksw2.h"
#include "ksw2_wide_simd.h"

/* width-specific compare/select so the d-flag chain is written once */
#if KSW_VBYTES == 64
typedef __mmask64 ksw_mask;
#  define KSW_CMPGT(a, b)      _mm512_cmpgt_epi8_mask((a), (b))   /* a > b ? 1:0 per lane */
#  define KSW_SEL(g, val, d)   _mm512_mask_mov_epi8((d), (g), (val)) /* g ? val : d */
#  define KSW_MVAL(g, v)       _mm512_maskz_set1_epi8((g), (v))      /* g ? v   : 0 */
#  define KSW_MVEC(g, a)       _mm512_maskz_mov_epi8((g), (a))       /* g ? a   : 0 */
#else  /* AVX2 */
typedef __m256i ksw_mask;
#  define KSW_CMPGT(a, b)      _mm256_cmpgt_epi8((a), (b))
#  define KSW_SEL(g, val, d)   _mm256_blendv_epi8((d), (val), (g))
#  define KSW_MVAL(g, v)       _mm256_and_si256((g), KSW_SET1(v))
#  define KSW_MVEC(g, a)       _mm256_and_si256((g), (a))
#endif

void extd2_ref_impl(void *km, int qlen, const uint8_t *query, int tlen,
                    const uint8_t *target, int8_t m, const int8_t *mat,
                    int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop,
                    int end_bonus, int flag, ksw_extz_t *ez);

typedef struct { ksw_vec q, q2, qe, qe2, zero, schi; } dconsts;  /* schi = mat[0] cap */

/* one wide chunk of the left-alignment cigar recurrence at byte offset tb */
static inline void d_step_wide(const dconsts *c, uint8_t *u, uint8_t *v,
                               uint8_t *x, uint8_t *y, uint8_t *x2, uint8_t *y2,
                               const uint8_t *s, uint8_t *pr, int tb,
                               int *cx, int *cv, int *cx2) {
    ksw_vec z   = KSW_LDU(s + tb);
    ksw_vec xt1 = KSW_LDU(x + tb);  int nx = x[tb + KSW_VBYTES - 1];
    xt1 = KSW_OR(ksw_vslli1(xt1), ksw_setlane0(*cx));  *cx = nx;
    ksw_vec vt1 = KSW_LDU(v + tb);  int nv = v[tb + KSW_VBYTES - 1];
    vt1 = KSW_OR(ksw_vslli1(vt1), ksw_setlane0(*cv));  *cv = nv;
    ksw_vec a  = KSW_ADD8(xt1, vt1);
    ksw_vec ut = KSW_LDU(u + tb);
    ksw_vec b  = KSW_ADD8(KSW_LDU(y + tb), ut);
    ksw_vec x2t1 = KSW_LDU(x2 + tb);  int nx2 = x2[tb + KSW_VBYTES - 1];
    x2t1 = KSW_OR(ksw_vslli1(x2t1), ksw_setlane0(*cx2));  *cx2 = nx2;
    ksw_vec a2 = KSW_ADD8(x2t1, vt1);
    ksw_vec b2 = KSW_ADD8(KSW_LDU(y2 + tb), ut);
    /* z = max(z,a,b,a2,b2) capped at mat[0]; d records the argmax (1..4) */
    ksw_vec d = KSW_MVAL(KSW_CMPGT(a, z), 1);  z = KSW_MAXS8(z, a);
    d = KSW_SEL(KSW_CMPGT(b,  z), KSW_SET1(2), d);  z = KSW_MAXS8(z, b);
    d = KSW_SEL(KSW_CMPGT(a2, z), KSW_SET1(3), d);  z = KSW_MAXS8(z, a2);
    d = KSW_SEL(KSW_CMPGT(b2, z), KSW_SET1(4), d);  z = KSW_MAXS8(z, b2);
    z = KSW_MINS8(z, c->schi);
    /* block2: store u,v; update a,b,a2,b2 */
    KSW_STU(u + tb, KSW_SUB8(z, vt1));
    KSW_STU(v + tb, KSW_SUB8(z, ut));
    ksw_vec tmp = KSW_SUB8(z, c->q);  a = KSW_SUB8(a, tmp);  b = KSW_SUB8(b, tmp);
    tmp = KSW_SUB8(z, c->q2);         a2 = KSW_SUB8(a2, tmp); b2 = KSW_SUB8(b2, tmp);
    /* store x,y,x2,y2 (positive part minus open+extend) and the open-gap d bits */
    ksw_mask g;
    g = KSW_CMPGT(a,  c->zero); KSW_STU(x  + tb, KSW_SUB8(KSW_MVEC(g, a),  c->qe));  d = KSW_OR(d, KSW_MVAL(g, 0x08));
    g = KSW_CMPGT(b,  c->zero); KSW_STU(y  + tb, KSW_SUB8(KSW_MVEC(g, b),  c->qe));  d = KSW_OR(d, KSW_MVAL(g, 0x10));
    g = KSW_CMPGT(a2, c->zero); KSW_STU(x2 + tb, KSW_SUB8(KSW_MVEC(g, a2), c->qe2)); d = KSW_OR(d, KSW_MVAL(g, 0x20));
    g = KSW_CMPGT(b2, c->zero); KSW_STU(y2 + tb, KSW_SUB8(KSW_MVEC(g, b2), c->qe2)); d = KSW_OR(d, KSW_MVAL(g, 0x40));
    KSW_STU(pr + tb, d);
}

/* 16-byte SSE tail (SSE4.1): the wide chunk is 32/64 bytes, so the 16-aligned
 * band remainder (one or two 16-byte blocks) is finished here. */
static inline void d_step_128(int q, int q2, int qe, int qe2, int schi,
                              uint8_t *u, uint8_t *v, uint8_t *x, uint8_t *y,
                              uint8_t *x2, uint8_t *y2, const uint8_t *s,
                              uint8_t *pr, int tb, int *cx, int *cv, int *cx2) {
    __m128i q_ = _mm_set1_epi8(q), q2_ = _mm_set1_epi8(q2);
    __m128i qe_ = _mm_set1_epi8(qe), qe2_ = _mm_set1_epi8(qe2);
    __m128i zero_ = _mm_setzero_si128(), schi_ = _mm_set1_epi8(schi);
    __m128i z = _mm_loadu_si128((const __m128i*)(s + tb));
    __m128i xt1 = _mm_loadu_si128((const __m128i*)(x + tb)); int nx = x[tb+15];
    xt1 = _mm_or_si128(_mm_slli_si128(xt1,1), _mm_cvtsi32_si128((unsigned char)*cx)); *cx = nx;
    __m128i vt1 = _mm_loadu_si128((const __m128i*)(v + tb)); int nv = v[tb+15];
    vt1 = _mm_or_si128(_mm_slli_si128(vt1,1), _mm_cvtsi32_si128((unsigned char)*cv)); *cv = nv;
    __m128i a = _mm_add_epi8(xt1, vt1);
    __m128i ut = _mm_loadu_si128((const __m128i*)(u + tb));
    __m128i b = _mm_add_epi8(_mm_loadu_si128((const __m128i*)(y + tb)), ut);
    __m128i x2t1 = _mm_loadu_si128((const __m128i*)(x2 + tb)); int nx2 = x2[tb+15];
    x2t1 = _mm_or_si128(_mm_slli_si128(x2t1,1), _mm_cvtsi32_si128((unsigned char)*cx2)); *cx2 = nx2;
    __m128i a2 = _mm_add_epi8(x2t1, vt1);
    __m128i b2 = _mm_add_epi8(_mm_loadu_si128((const __m128i*)(y2 + tb)), ut);
    __m128i d = _mm_and_si128(_mm_cmpgt_epi8(a, z), _mm_set1_epi8(1));      z = _mm_max_epi8(z, a);
    d = _mm_blendv_epi8(d, _mm_set1_epi8(2), _mm_cmpgt_epi8(b,  z));        z = _mm_max_epi8(z, b);
    d = _mm_blendv_epi8(d, _mm_set1_epi8(3), _mm_cmpgt_epi8(a2, z));        z = _mm_max_epi8(z, a2);
    d = _mm_blendv_epi8(d, _mm_set1_epi8(4), _mm_cmpgt_epi8(b2, z));        z = _mm_max_epi8(z, b2);
    z = _mm_min_epi8(z, schi_);
    _mm_storeu_si128((__m128i*)(u + tb), _mm_sub_epi8(z, vt1));
    _mm_storeu_si128((__m128i*)(v + tb), _mm_sub_epi8(z, ut));
    __m128i tmp = _mm_sub_epi8(z, q_);  a = _mm_sub_epi8(a, tmp);  b = _mm_sub_epi8(b, tmp);
    tmp = _mm_sub_epi8(z, q2_);         a2 = _mm_sub_epi8(a2, tmp); b2 = _mm_sub_epi8(b2, tmp);
    __m128i g;
    g = _mm_cmpgt_epi8(a,  zero_); _mm_storeu_si128((__m128i*)(x  + tb), _mm_sub_epi8(_mm_and_si128(g,a),  qe_));  d = _mm_or_si128(d, _mm_and_si128(g, _mm_set1_epi8(0x08)));
    g = _mm_cmpgt_epi8(b,  zero_); _mm_storeu_si128((__m128i*)(y  + tb), _mm_sub_epi8(_mm_and_si128(g,b),  qe_));  d = _mm_or_si128(d, _mm_and_si128(g, _mm_set1_epi8(0x10)));
    g = _mm_cmpgt_epi8(a2, zero_); _mm_storeu_si128((__m128i*)(x2 + tb), _mm_sub_epi8(_mm_and_si128(g,a2), qe2_)); d = _mm_or_si128(d, _mm_and_si128(g, _mm_set1_epi8(0x20)));
    g = _mm_cmpgt_epi8(b2, zero_); _mm_storeu_si128((__m128i*)(y2 + tb), _mm_sub_epi8(_mm_and_si128(g,b2), qe2_)); d = _mm_or_si128(d, _mm_and_si128(g, _mm_set1_epi8(0x40)));
    _mm_storeu_si128((__m128i*)(pr + tb), d);
}

/* wide match/mismatch score-setup chunk (may overshoot into the 128-byte pad) */
static inline void score_step_wide(uint8_t *s, const uint8_t *sf, const uint8_t *qrr,
                                   ksw_vec m1, ksw_vec mch, ksw_vec mis, ksw_vec scN) {
    ksw_vec sq = KSW_LDU(sf), sti = KSW_LDU(qrr);
#if KSW_VBYTES == 64
    __mmask64 meq = _mm512_cmpeq_epi8_mask(sq, sti);
    __mmask64 mN  = _kor_mask64(_mm512_cmpeq_epi8_mask(sq, m1), _mm512_cmpeq_epi8_mask(sti, m1));
    ksw_vec r = _mm512_mask_blend_epi8(meq, mis, mch);
    r = _mm512_mask_blend_epi8(mN, r, scN);
#else
    __m256i mask = _mm256_or_si256(_mm256_cmpeq_epi8(sq, m1), _mm256_cmpeq_epi8(sti, m1));
    __m256i r = _mm256_blendv_epi8(mis, mch, _mm256_cmpeq_epi8(sq, sti));
    r = _mm256_blendv_epi8(r, scN, mask);
#endif
    KSW_STU(s, r);
}

void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target,
                   int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
                   int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
    /* widen only minibwa's gap-fill path: cigar + approx + left + exact-sc + !extz */
    if ((flag & KSW_EZ_SCORE_ONLY) || !(flag & KSW_EZ_APPROX_MAX) ||
        (flag & KSW_EZ_GENERIC_SC) || (flag & KSW_EZ_RIGHT) || (flag & KSW_EZ_EXTZ_ONLY)) {
        extd2_ref_impl(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w,
                       zdrop, end_bonus, flag, ez);
        return;
    }

    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr;
    int max_sc, min_sc, long_thres, long_diff;
    int32_t H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0, *u8, *v8, *x8, *y8, *x28, *y28, *s8, *pp;
    int mat0 = mat[0], matmis = mat[1];

    ksw_reset_extz(ez);
    if (m <= 1 || qlen <= 0 || tlen <= 0) return;
    if (q2 + e2 < q + e) { int8_t tt; tt=q,q=q2,q2=tt; tt=e,e=e2,e2=tt; }
    int scN = mat[m*m-1] == 0 ? -e2 : mat[m*m-1];

    if (w < 0) w = tlen > qlen ? tlen : qlen;
    wl = wr = w;
    tlen_ = (tlen + 15) / 16;
    n_col_ = qlen < tlen ? qlen : tlen;
    n_col_ = ((n_col_ < w + 1 ? n_col_ : w + 1) + 15) / 16 + 1;
    qlen_ = (qlen + 15) / 16;
    for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
        max_sc = max_sc > mat[t] ? max_sc : mat[t];
        min_sc = min_sc < mat[t] ? min_sc : mat[t];
    }
    if (-min_sc > 2 * (q + e)) return;
    (void)max_sc;

    long_thres = e != e2 ? (q2 - q) / (e - e2) - 1 : 0;
    if (q2 + e2 + long_thres * e2 > q + e + long_thres * e) ++long_thres;
    long_diff = long_thres * (e - e2) - (q2 - q) - e2;

    dconsts dc;
    dc.q = KSW_SET1((int8_t)q); dc.q2 = KSW_SET1((int8_t)q2);
    dc.qe = KSW_SET1((int8_t)(q + e)); dc.qe2 = KSW_SET1((int8_t)(q2 + e2));
    dc.zero = KSW_SET1(0); dc.schi = KSW_SET1((int8_t)mat0);
    ksw_vec Wm1 = KSW_SET1((int8_t)(m - 1)), Wmch = KSW_SET1((int8_t)mat0);
    ksw_vec Wmis = KSW_SET1((int8_t)matmis), WN = KSW_SET1((int8_t)scN);

    /* padded arrays (128B pad per array so wide chunks can overshoot safely) */
    int tcap = tlen_ * 16 + 128, qcap = qlen_ * 16 + 128;
    mem = (uint8_t*)kmalloc(km, (size_t)7 * tcap);
    u8 = mem; v8 = u8+tcap; x8 = v8+tcap; y8 = x8+tcap; x28 = y8+tcap; y28 = x28+tcap;
    s8 = y28+tcap;
    sf = (uint8_t*)kmalloc(km, (size_t)tcap); qr = (uint8_t*)kmalloc(km, (size_t)qcap);
    memset(u8,  -q  - e,  tcap); memset(v8,  -q  - e,  tcap);
    memset(x8,  -q  - e,  tcap); memset(y8,  -q  - e,  tcap);
    memset(x28, -q2 - e2, tcap); memset(y28, -q2 - e2, tcap);
    memset(s8, 0, tcap); memset(sf, 0, tcap); memset(qr, 0, qcap);
    /* cigar direction matrix + per-row offsets */
    mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
    pp = mem2;
    off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
    off_end = off + qlen + tlen - 1;

    for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
    memcpy(sf, target, tlen);

    for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
        int st = 0, en = tlen - 1, st0, en0;
        int8_t x1, x21, v1;
        uint8_t *qrr = qr + (qlen - 1 - r);
        if (st < r - qlen + 1) st = r - qlen + 1;
        if (en > r) en = r;
        if (st < (r-wr+1)>>1) st = (r-wr+1)>>1;
        if (en > (r+wl)>>1)   en = (r+wl)>>1;
        if (st > en) { ez->zdropped = 1; break; }
        st0 = st, en0 = en;
        st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
        if (st > 0) {
            if (st - 1 >= last_st && st - 1 <= last_en)
                x1 = x8[st-1], x21 = x28[st-1], v1 = v8[st-1];
            else { x1 = -q - e, x21 = -q2 - e2; v1 = -q - e; }
        } else {
            x1 = -q - e, x21 = -q2 - e2;
            v1 = r == 0 ? -q - e : r < long_thres ? -e : r == long_thres ? long_diff : -e2;
        }
        if (en >= r) {
            y8[r] = -q - e, y28[r] = -q2 - e2;
            u8[r] = r == 0 ? -q - e : r < long_thres ? -e : r == long_thres ? long_diff : -e2;
        }
        /* set scores (wide; may overshoot into pad) */
        for (t = st0; t <= en0; t += KSW_VBYTES)
            score_step_wide(s8 + t, sf + t, qrr + t, Wm1, Wmch, Wmis, WN);
        /* cigar core: wide chunks over 16-aligned [st,en], then 16-byte tail */
        {
            uint8_t *pr = pp + ((size_t)r * n_col_ - (st/16)) * 16;  /* byte base */
            int cx = x1, cv = v1, cx2 = x21, tb = st;
            off[r] = st, off_end[r] = en;
            for (; tb + KSW_VBYTES - 1 <= en; tb += KSW_VBYTES)
                d_step_wide(&dc, u8, v8, x8, y8, x28, y28, s8, pr, tb, &cx, &cv, &cx2);
            for (; tb + 15 <= en; tb += 16)
                d_step_128(q, q2, q+e, q2+e2, mat0, u8, v8, x8, y8, x28, y28, s8, pr, tb, &cx, &cv, &cx2);
        }
        /* approximate max (minibwa's gap-fill uses APPROX_MAX). NB: in extd2 the
         * u/v deltas are SIGNED (init -q-e); read them as int8_t, not uint8_t. */
        if (r > 0) {
            if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
                int32_t d0 = (int8_t)v8[last_H0_t], d1 = (int8_t)u8[last_H0_t + 1];
                if (d0 > d1) H0 += d0; else H0 += d1, ++last_H0_t;
            } else if (last_H0_t >= st0 && last_H0_t <= en0) {
                H0 += (int8_t)v8[last_H0_t];
            } else { ++last_H0_t, H0 += (int8_t)u8[last_H0_t]; }
        } else H0 = (int8_t)v8[0] - qe, last_H0_t = 0;
        if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
        if (r == qlen + tlen - 2 && en0 == tlen - 1) ez->score = H0;
        last_st = st, last_en = en;
    }
    kfree(km, sf); kfree(km, qr); kfree(km, mem);
    /* backtrack (gap-fill: !zdropped && !EXTZ_ONLY) */
    if (!ez->zdropped) {
        int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
        ksw_backtrack(km, 1, rev_cigar, 0, pp, off, off_end, n_col_*16, tlen-1, qlen-1,
                      &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    }
    kfree(km, mem2); kfree(km, off);
    (void)end_bonus;
}
