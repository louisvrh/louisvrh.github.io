// Physics core: elastic X-ray scattering from an N x N square lattice (a = 1).
// Units: lengths in a; wavevectors in 1/a; K = |k_in| is expressed in units of pi/a.
// Reciprocal lattice of the square lattice: G = (2pi/a)(h1,h2).
const Core = (function () {
  const PI = Math.PI, DEG = PI / 180;
  const KMIN = 0.3, KMAX = 6;

  function kin(K, phiDeg) { return [K * PI * Math.cos(phiDeg * DEG), K * PI * Math.sin(phiDeg * DEG)]; }

  // |sum_{j=0}^{N-1} e^{-i x j}|^2 / N^2  (normalised Dirichlet / Laue function), exact closed form
  function dirichlet2(N, x) {
    const s = Math.sin(x / 2);
    if (Math.abs(s) < 1e-12) return 1;
    const r = Math.sin(N * x / 2) / (N * s);
    return r * r;
  }

  // Normalised intensity |sum_R e^{-i q.R}|^2 / N^4 into direction theta (rad), |k_out| = |k_in|.
  // q = k_out - k_in. The 2D sum over R=(i,j) factorises into two 1D Laue functions.
  function intensity(N, K, phiDeg, theta) {
    const k = K * PI, kv = kin(K, phiDeg);
    const qx = k * Math.cos(theta) - kv[0], qy = k * Math.sin(theta) - kv[1];
    return dirichlet2(N, qx) * dirichlet2(N, qy);
  }

  // Direct double-sum over all N^2 atoms (used only to verify the closed form)
  function bruteForce(N, K, phiDeg, theta) {
    const k = K * PI, kv = kin(K, phiDeg);
    const qx = k * Math.cos(theta) - kv[0], qy = k * Math.sin(theta) - kv[1];
    let re = 0, im = 0;
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      const ph = -(qx * i + qy * j);
      re += Math.cos(ph); im += Math.sin(ph);
    }
    return (re * re + im * im) / Math.pow(N, 4);
  }

  function gcd(a, b) { a = Math.abs(a); b = Math.abs(b); while (b) { const t = a % b; a = b; b = t; } return a; }

  function describe(h1, h2, K, phiDeg) {
    const kv = kin(K, phiDeg), kmag = K * PI;
    const ox = kv[0] + 2 * PI * h1, oy = kv[1] + 2 * PI * h2, omag = Math.hypot(ox, oy);
    const defl = Math.abs(Math.atan2(kv[0] * oy - kv[1] * ox, kv[0] * ox + kv[1] * oy)) / DEG;
    const n = gcd(h1, h2);
    const p1 = n ? h1 / n : 0, p2 = n ? h2 / n : 0;
    return {
      h1: h1, h2: h2, n: n, p1: p1, p2: p2,
      d: n ? 1 / Math.hypot(p1, p2) : Infinity,           // plane spacing in units of a
      eps: (omag - kmag) / PI,                            // |k_in+G| - |k_in|, units pi/a (0 = on the Ewald circle)
      dirRad: Math.atan2(oy, ox), dirDeg: ((Math.atan2(oy, ox) / DEG) % 360 + 360) % 360,
      deflDeg: defl, thetaB: defl / 2, lambda: 2 * PI / kmag
    };
  }

  // Reciprocal lattice points exactly on the Ewald circle (2 k.G + G^2 = 0). Includes G=0 (forward beam).
  function exactBeams(K, phiDeg, tol) {
    tol = tol || 1e-9;
    const m = Math.ceil(K) + 1, kv = kin(K, phiDeg), kmag = K * PI, out = [];
    for (let h1 = -m; h1 <= m; h1++) for (let h2 = -m; h2 <= m; h2++) {
      const ox = kv[0] + 2 * PI * h1, oy = kv[1] + 2 * PI * h2;
      if (Math.abs(Math.hypot(ox, oy) - kmag) < tol) out.push(describe(h1, h2, K, phiDeg));
    }
    return out;
  }

  // Main-lobe envelope of the Laue function: (sin(u/2)/(u/2))^2, zero beyond the first zero |u| = 2pi
  function lobe(u) { if (Math.abs(u) >= 2 * PI) return 0; if (Math.abs(u) < 1e-12) return 1; const r = Math.sin(u / 2) / (u / 2); return r * r; }

  // Finite-N beams: G != 0 whose lattice point lies within the finite-size peak (|delta q| < 2pi/(N a) per axis)
  // of the Ewald circle. delta = q_on_circle - G is the radial mismatch; I = peak intensity relative to N^4.
  function visibleBeams(N, K, phiDeg, thr) {
    thr = thr === undefined ? 0.05 : thr;
    const m = Math.ceil(K) + 1, kv = kin(K, phiDeg), kmag = K * PI, out = [];
    for (let h1 = -m; h1 <= m; h1++) for (let h2 = -m; h2 <= m; h2++) {
      if (h1 === 0 && h2 === 0) continue;
      const ox = kv[0] + 2 * PI * h1, oy = kv[1] + 2 * PI * h2, om = Math.hypot(ox, oy);
      const dx = (kmag - om) * ox / om, dy = (kmag - om) * oy / om;
      const I = lobe(N * dx) * lobe(N * dy);
      if (I >= thr) { const b = describe(h1, h2, K, phiDeg); b.I = I; out.push(b); }
    }
    out.sort(function (a, b) { return a.dirDeg - b.dirDeg; });
    return out;
  }

  // Polar profile of intensity, M uniform samples of the outgoing direction
  function polar(N, K, phiDeg, M) {
    const th = new Float64Array(M), f = new Float64Array(M);
    for (let i = 0; i < M; i++) { th[i] = 2 * PI * i / M; f[i] = intensity(N, K, phiDeg, th[i]); }
    return { theta: th, f: f };
  }

  // Nearest |k_in| (keeping the angle) for which some G != 0 lies exactly on the Ewald circle: K (h.u) = -h^2
  function snapK(K, phiDeg) {
    const c = Math.cos(phiDeg * DEG), s = Math.sin(phiDeg * DEG);
    let best = null;
    for (let h1 = -8; h1 <= 8; h1++) for (let h2 = -8; h2 <= 8; h2++) {
      if (!h1 && !h2) continue;
      const d = h1 * c + h2 * s;
      if (d > -1e-9) continue;
      const Kc = -(h1 * h1 + h2 * h2) / d;
      if (Kc < KMIN || Kc > KMAX) continue;
      if (!best || Math.abs(Kc - K) < Math.abs(best.K - K)) best = { K: Kc, h1: h1, h2: h2 };
    }
    return best;
  }

  // Nearest angle (keeping |k_in|) for which some G != 0 lies exactly on the circle: cos(phi - atan2(h)) = -|h|/K
  function snapPhi(K, phiDeg) {
    let best = null;
    for (let h1 = -8; h1 <= 8; h1++) for (let h2 = -8; h2 <= 8; h2++) {
      if (!h1 && !h2) continue;
      const r = Math.hypot(h1, h2);
      if (r > K) continue;
      const base = Math.atan2(h2, h1) / DEG, a = Math.acos(-r / K) / DEG;
      [base + a, base - a].forEach(function (p) {
        p = ((p % 360) + 360) % 360;
        let dd = Math.abs(p - phiDeg) % 360; if (dd > 180) dd = 360 - dd;
        if (!best || dd < best.dd) best = { phi: p, h1: h1, h2: h2, dd: dd };
      });
    }
    return best;
  }

  // Half-width (rad) of the forward peak at half maximum
  function forwardHalfWidth(N, K, phiDeg) {
    const ph = phiDeg * DEG, step = 2e-6;
    let t = 0, prev = 1;
    while (t < 1) {
      t += step;
      const v = intensity(N, K, phiDeg, ph + t);
      if (v < 0.5) return t - step + step * (prev - 0.5) / (prev - v);
      prev = v;
    }
    return NaN;
  }


  // ---- G-triangle / Bragg-plane overlay geometry (screen pixels: x right, y down, origin = triangle tail) ----
  // pxU = pixels per unit of 2pi/a (shared by k_in, G, k_out). ad = pixels per lattice constant a.
  function pickBeam(beams, key) {
    if (!beams || !beams.length) return null;
    if (key) { for (let i = 0; i < beams.length; i++) if (beams[i].h1 + ',' + beams[i].h2 === key) return beams[i]; }
    let best = beams[0];
    for (let i = 1; i < beams.length; i++) {
      const b = beams[i], gb = Math.hypot(b.h1, b.h2), gbest = Math.hypot(best.h1, best.h2);
      const wb = b.I * (b.W === undefined ? 1 : b.W), wbest = best.I * (best.W === undefined ? 1 : best.W);
      if (wb > wbest + 1e-9 || (Math.abs(wb - wbest) <= 1e-9 && gb < gbest - 1e-12)) best = b;
    }
    return best;
  }
  function overlay(b, K, phiDeg, pxU, ad, N) {
    const pp = pxU / (2 * PI), kv = kin(K, phiDeg), kmag = K * PI;
    const kinS = [kv[0] * pp, -kv[1] * pp];
    const G = [2 * PI * b.h1 * pp, -2 * PI * b.h2 * pp];
    const koutS = [Math.cos(b.dirRad) * kmag * pp, -Math.sin(b.dirRad) * kmag * pp];
    const tipG = [kinS[0] + G[0], kinS[1] + G[1]];
    const pm = Math.hypot(b.p1, b.p2);
    const nS = [b.p1 / pm, -b.p2 / pm], tS = [b.p2 / pm, b.p1 / pm];
    let c = -(b.p1 + b.p2) * (N - 1) / 2; c = c - Math.floor(c);
    const Gmag = 2 * PI * Math.hypot(b.h1, b.h2);                 // |G| in 1/a
    const gx = 2 * PI * b.h1 + kv[0], gy = 2 * PI * b.h2 + kv[1];
    return { kin: kinS, G: G, kout: koutS, tipKin: kinS, tipG: tipG, tipKout: koutS,
             gap: [koutS[0] - tipG[0], koutS[1] - tipG[1]],
             n: nS, t: tS, spacingPx: ad / pm, offset: c, pm: pm,
             Gmag: Gmag, d: b.d, thetaDeg: b.thetaB, twoThetaDeg: b.deflDeg,
             lambda: 2 * PI / kmag,
             check: (Math.hypot(gx, gy) - kmag) / PI,              // |k_in+G| - |k_out| in pi/a
             sinThetaG: Gmag / (2 * kmag) };
  }
  // Bragg-plane line segments (relative to centre) clipped to a disc of radius Rc
  function planeSegs(ov, Rc) {
    const sp = ov.spacingPx, out = [];
    if (!(sp > 0)) return out;
    const mm = Math.ceil(Rc / sp) + 2;
    for (let m = -mm; m <= mm; m++) {
      const s = (m + ov.offset) * sp; if (Math.abs(s) >= Rc) continue;
      const h = Math.sqrt(Rc * Rc - s * s);
      out.push({ s: s, seg: [ov.n[0] * s - ov.t[0] * h, ov.n[1] * s - ov.t[1] * h, ov.n[0] * s + ov.t[0] * h, ov.n[1] * s + ov.t[1] * h] });
    }
    return out;
  }


  // ---- Bragg-plane family: atoms R = (i,j) a (i,j = 0..N-1, j up) with g.R = 2pi v, v = p1 i + p2 j (integer) ----
  // m = v - v0 is the plane index; v0 = plane closest to the crystal centre (the overlay origin). R is then measured
  // from an atom on plane m = 0, so g.R = 2pi m.
  function planeFamily(p1, p2, N) {
    const c0 = (N - 1) / 2, v0 = Math.round((p1 + p2) * c0), atoms = [], counts = {};
    let vmin = Infinity, vmax = -Infinity;
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
      const v = p1 * i + p2 * j; atoms.push({ i: i, j: j, v: v, m: v - v0 });
      counts[v - v0] = (counts[v - v0] || 0) + 1; if (v < vmin) vmin = v; if (v > vmax) vmax = v;
    }
    return { v0: v0, mmin: vmin - v0, mmax: vmax - v0, atoms: atoms, counts: counts, nPlanes: vmax - vmin + 1 };
  }
  function atomsOnPlane(fam, m) { return fam.atoms.filter(function (a) { return a.m === m; }); }
  // signed distance (px) of plane m from the centre along n (screen normal), for pixel pitch ad
  function planeOffsetPx(fam, p1, p2, N, m, ad) { return (fam.v0 + m - (p1 + p2) * (N - 1) / 2) / Math.hypot(p1, p2) * ad; }
  // nearest plane index to a signed distance s (px) along n, clamped to existing planes
  function planeFromOffset(fam, p1, p2, N, s, ad) {
    const m = Math.round(s * Math.hypot(p1, p2) / ad + (p1 + p2) * (N - 1) / 2) - fam.v0;
    return Math.max(fam.mmin, Math.min(fam.mmax, m));
  }


  // ================= Non-lattice structures: direct sum I(q) = |sum_j f_j e^{-i q.r_j}|^2 =================
  // Lengths in units of a; the N x N lattice sites are R = ((i-(N-1)/2), (j-(N-1)/2)) a, site index s = i*N + j (j up).
  // Intensity is returned normalised by (sum f)^2 (so the forward beam is 1); Itilde = I / N_sites^2 = Ihat (sumF/N_sites)^2.
  function rng(seed) { let a = (seed * 2654435761 + 1013904223) >>> 0; return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
  function gauss(r) { let u = 0; while (u < 1e-12) u = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * PI * r()); }

  function amorphousPositions(N, r, minSep) {
    const n = N * N, H = N / 2, x = new Float64Array(n), y = new Float64Array(n);
    for (let k = 0; k < n; k++) { x[k] = (r() * 2 - 1) * H; y[k] = (r() * 2 - 1) * H; }
    if (n < 2) return { x: x, y: y, minSep: Infinity, iters: 0 };
    const cs = minSep, nc = Math.max(1, Math.ceil(N / cs)), cells = new Array(nc * nc);
    let iters = 0, ms = 0;
    for (; iters < 2000; iters++) {
      for (let c = 0; c < cells.length; c++) cells[c] = [];
      for (let k = 0; k < n; k++) { const cx = Math.min(nc - 1, Math.floor((x[k] + H) / cs)), cy = Math.min(nc - 1, Math.floor((y[k] + H) / cs)); cells[cy * nc + cx].push(k); }
      let moved = false;
      for (let k = 0; k < n; k++) {
        const cx = Math.min(nc - 1, Math.floor((x[k] + H) / cs)), cy = Math.min(nc - 1, Math.floor((y[k] + H) / cs));
        for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
          const ax = cx + dx, ay = cy + dy; if (ax < 0 || ay < 0 || ax >= nc || ay >= nc) continue;
          const lst = cells[ay * nc + ax];
          for (let q = 0; q < lst.length; q++) { const m = lst[q]; if (m <= k) continue;
            let ex = x[m] - x[k], ey = y[m] - y[k], d = Math.hypot(ex, ey);
            if (d < minSep) { if (d < 1e-9) { ex = r() - 0.5; ey = r() - 0.5; d = Math.hypot(ex, ey); }
              const push = (minSep - d) / 2 * 1.02 + 1e-9, ux = ex / d, uy = ey / d;
              x[k] -= ux * push; y[k] -= uy * push; x[m] += ux * push; y[m] += uy * push; moved = true; } } }
        x[k] = Math.max(-H, Math.min(H, x[k])); y[k] = Math.max(-H, Math.min(H, y[k]));
      }
      if (!moved) break;
    }
    // exact minimum separation (brute force)
    ms = Infinity; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = Math.hypot(x[i] - x[j], y[i] - y[j]); if (d < ms) ms = d; }
    return { x: x, y: y, minSep: ms, iters: iters };
  }

  // par = {p (0..1), sigma (a), tx, ty (a), rho = f2/f1}
  function makeStructure(mode, N, par, seed) {
    par = par || {}; const ns = N * N, c0 = (N - 1) / 2, r = rng(seed || 1);
    const sx = new Float64Array(ns), sy = new Float64Array(ns), ord = new Uint8Array(ns).fill(1);
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { sx[i * N + j] = i - c0; sy[i * N + j] = j - c0; }
    // fixed random streams (independent of the sliders) so that sliders move atoms smoothly
    const rx = new Float64Array(ns), ry = new Float64Array(ns), gx = new Float64Array(ns), gy = new Float64Array(ns), rank = new Float64Array(ns), perm = new Int32Array(ns);
    for (let s = 0; s < ns; s++) { rx[s] = (r() * 2 - 1) * N / 2; ry[s] = (r() * 2 - 1) * N / 2; gx[s] = gauss(r); gy[s] = gauss(r); perm[s] = s; }
    for (let s = ns - 1; s > 0; s--) { const t = Math.floor(r() * (s + 1)), tmp = perm[s]; perm[s] = perm[t]; perm[t] = tmp; }
    for (let s = 0; s < ns; s++) rank[perm[s]] = s;
    let x = [], y = [], f = [], kind = [], site = [];
    const push = function (xx, yy, ff, kk, ss) { x.push(xx); y.push(yy); f.push(ff); kind.push(kk); site.push(ss); };
    let nOrd = ns;
    if (mode === 'semi') {
      const p = Math.max(0, Math.min(1, par.p === undefined ? 1 : par.p)); nOrd = Math.round(p * ns);
      for (let s = 0; s < ns; s++) { if (rank[s] < nOrd) push(sx[s], sy[s], 1, 0, s); else { ord[s] = 0; push(rx[s], ry[s], 1, 1, -1); } }
    } else if (mode === 'amorph') {
      const am = amorphousPositions(N, rng((seed || 1) + 7919), 0.8); nOrd = 0; ord.fill(0);
      for (let s = 0; s < ns; s++) push(am.x[s], am.y[s], 1, 1, -1);
    } else if (mode === 'thermal') {
      const sg = par.sigma || 0;
      for (let s = 0; s < ns; s++) { sx[s] += sg * gx[s]; sy[s] += sg * gy[s]; push(sx[s], sy[s], 1, 3, s); }
    } else if (mode === 'basis') {
      const tx = par.tx === undefined ? 0.5 : par.tx, ty = par.ty === undefined ? 0.5 : par.ty, rho = par.rho === undefined ? 1 : par.rho;
      for (let s = 0; s < ns; s++) push(sx[s], sy[s], 1, 0, s);
      for (let s = 0; s < ns; s++) push(sx[s] + tx, sy[s] + ty, rho, 2, s);
    } else { for (let s = 0; s < ns; s++) push(sx[s], sy[s], 1, 0, s); }
    let sumF = 0; for (let k = 0; k < f.length; k++) sumF += f[k];
    return { mode: mode, N: N, nSites: ns, nOrd: nOrd, n: x.length, x: Float64Array.from(x), y: Float64Array.from(y), f: Float64Array.from(f),
             kind: Uint8Array.from(kind), site: Int32Array.from(site), sx: sx, sy: sy, ord: ord, sumF: sumF, par: par };
  }
  function amplitudeDirect(st, K, phiDeg, theta) {
    const k = K * PI, kv = kin(K, phiDeg), qx = k * Math.cos(theta) - kv[0], qy = k * Math.sin(theta) - kv[1];
    let re = 0, im = 0; const x = st.x, y = st.y, f = st.f;
    for (let a = 0; a < st.n; a++) { const ph = -(qx * x[a] + qy * y[a]); re += f[a] * Math.cos(ph); im += f[a] * Math.sin(ph); }
    return [re, im];
  }
  // |sum f e^{-iq.r}|^2 / (sum f)^2
  function intensityDirect(st, K, phiDeg, theta) { const A = amplitudeDirect(st, K, phiDeg, theta); return (A[0] * A[0] + A[1] * A[1]) / (st.sumF * st.sumF); }
  function polarDirect(st, K, phiDeg, M) {
    const th = new Float64Array(M), fo = new Float64Array(M);
    for (let i = 0; i < M; i++) { th[i] = 2 * PI * i / M; fo[i] = intensityDirect(st, K, phiDeg, th[i]); }
    return { theta: th, f: fo };
  }
  // structure factor modulus |1 + rho e^{i G.tau}| for G = (2pi/a)(h1,h2), tau = (tx,ty) a
  function structFactor(h1, h2, tx, ty, rho) { const ph = 2 * PI * (h1 * tx + h2 * ty); return Math.hypot(1 + rho * Math.cos(ph), rho * Math.sin(ph)); }

  // Off-peak statistics of a polar profile. Excludes +-win unit widths (unit = 2/(N K) rad) around the forward beam and all listed beams.
  // Returns means in Itilde = I/N_sites^2 units. median/ln2 is the speckle-robust estimate of the mean diffuse level.
  function offPeak(st, P, K, phiDeg, beams, win) {
    win = win || 4; const N = st.N, w = Math.max(win * 2 / (N * K), 1.2 * DEG), M = P.f.length, norm = Math.pow(st.sumF / st.nSites, 2);
    const dirs = [phiDeg * DEG].concat((beams || []).map(function (b) { return b.dirRad; }));
    const vals = [], qs = [];
    for (let i = 0; i < M; i++) { let ok = true; for (let d = 0; d < dirs.length; d++) { let dd = Math.abs(P.theta[i] - dirs[d]) % (2 * PI); if (dd > PI) dd = 2 * PI - dd; if (dd < w) { ok = false; break; } }
      if (ok) { vals.push(P.f[i] * norm); const k = K * PI, kv = kin(K, phiDeg); qs.push(Math.hypot(k * Math.cos(P.theta[i]) - kv[0], k * Math.sin(P.theta[i]) - kv[1])); } }
    if (!vals.length) return { n: 0, mean: NaN, medianEst: NaN, meanQ2: NaN, qs: qs };
    let s = 0; for (let i = 0; i < vals.length; i++) s += vals[i];
    const sorted = vals.slice().sort(function (a, b) { return a - b; }), med = sorted[Math.floor(sorted.length / 2)];
    return { n: vals.length, mean: s / vals.length, medianEst: med / Math.LN2, qs: qs };
  }
  // Predicted diffuse level  Itilde_bg * N_sites  (= I_bg / N_sites) for the model, averaged over the same off-peak directions
  function predictedDiffuse(st, qs) {
    const m = st.mode;
    if (m === 'semi') { const pe = st.nOrd / st.nSites; return 1 - pe * pe; }                       // (1-p) random atoms + p(1-p) vacancies
    if (m === 'thermal') { const s2 = Math.pow(st.par.sigma || 0, 2); let a = 0; for (let i = 0; i < qs.length; i++) a += 1 - Math.exp(-qs[i] * qs[i] * s2); return qs.length ? a / qs.length : 0; }
    if (m === 'amorph') return NaN;
    return 0;
  }
  function predictedPeak(st, h1, h2) {   // expected Itilde = I/N_sites^2 at an exact Bragg direction (large-N average)
    const m = st.mode, ns = st.nSites, G2 = Math.pow(2 * PI, 2) * (h1 * h1 + h2 * h2);
    if (m === 'semi') { const pe = st.nOrd / ns; return pe * pe + (1 - pe * pe) / ns; }
    if (m === 'thermal') { const e = Math.exp(-G2 * Math.pow(st.par.sigma || 0, 2)); return e + (1 - e) / ns; }
    if (m === 'basis') { const F = structFactor(h1, h2, st.par.tx, st.par.ty, st.par.rho); return F * F; }
    return 1;
  }

  return { rng: rng, amorphousPositions: amorphousPositions, makeStructure: makeStructure, amplitudeDirect: amplitudeDirect, intensityDirect: intensityDirect, polarDirect: polarDirect, structFactor: structFactor, offPeak: offPeak, predictedDiffuse: predictedDiffuse, predictedPeak: predictedPeak, planeFamily: planeFamily, atomsOnPlane: atomsOnPlane, planeOffsetPx: planeOffsetPx, planeFromOffset: planeFromOffset, pickBeam: pickBeam, overlay: overlay, planeSegs: planeSegs, KMIN: KMIN, KMAX: KMAX, kin: kin, dirichlet2: dirichlet2, intensity: intensity, bruteForce: bruteForce,
           describe: describe, exactBeams: exactBeams, visibleBeams: visibleBeams, polar: polar,
           snapK: snapK, snapPhi: snapPhi, forwardHalfWidth: forwardHalfWidth };
})();
if (typeof module !== 'undefined') module.exports = Core;
