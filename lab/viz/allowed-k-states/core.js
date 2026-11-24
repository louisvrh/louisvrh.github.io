// Physics core: allowed k-states of a monatomic chain with Born-von Karman (ring) or box BCs.
// Units: a = 1, omega_D = sqrt(kappa/m) = 1, so L = N.
var Core = (function () {
  var TWO_PI = 2 * Math.PI;

  // Eq. 4.10: omega_k = 2 omega_D |sin(ka/2)|
  function omega(k) { return 2 * Math.abs(Math.sin(k / 2)); }

  // Eq. 5.2 range of p for -pi/a < k <= pi/a
  function pRange(N) {
    if (N % 2 === 1) { var m = (N - 1) / 2; return [-m, m]; }
    return [-N / 2 + 1, N / 2];
  }

  // Allowed states. periodic: k = p 2pi/(Na). box (fixed ends u_0 = u_N = 0): k = p pi/(Na), p = 1..N-1
  // (N-1 mobile atoms -> N-1 modes; "N" in the large-N limit).
  function allowedK(N, bc) {
    var out = [], p;
    if (bc === 'box') {
      for (p = 1; p <= N - 1; p++) out.push({ p: p, k: p * Math.PI / N });
    } else {
      var r = pRange(N);
      for (p = r[0]; p <= r[1]; p++) out.push({ p: p, k: p * TWO_PI / N });
    }
    return out;
  }
  function spacing(N, bc) { return (bc === 'box' ? Math.PI : TWO_PI) / N; }      // 2pi/L or pi/L
  function density(N, bc) { return (bc === 'box' ? N / Math.PI : N / TWO_PI); }  // L/2pi or L/pi

  // Displacement of atom n at time t in mode k, amplitude A.
  function modeU(n, t, k, A, bc) {
    var w = omega(k);
    if (bc === 'box') return A * Math.sin(k * n) * Math.cos(w * t);
    return A * Math.cos(k * n - w * t);
  }
  // Residual of  u_ddot_n = u_{n+1} - 2u_n + u_{n-1}  (m = kappa = 1), neighbours taken with
  // the ring identification n+N -> n (periodic) or the fixed-wall values (box). Max over atoms.
  function eomResidual(N, k, bc, t) {
    var w = omega(k), worst = 0, n, lo, hi, upl, umi, u, acc, rhs;
    lo = (bc === 'box') ? 1 : 0; hi = (bc === 'box') ? N - 1 : N - 1;
    for (n = lo; n <= hi; n++) {
      var np = (bc === 'box') ? n + 1 : (n + 1) % N;
      var nm = (bc === 'box') ? n - 1 : (n - 1 + N) % N;
      u = modeU(n, t, k, 1, bc); upl = modeU(np, t, k, 1, bc); umi = modeU(nm, t, k, 1, bc);
      acc = -w * w * u; rhs = upl - 2 * u + umi;
      worst = Math.max(worst, Math.abs(acc - rhs));
    }
    return worst;
  }
  // Eq. 5.1 check: |u_{n+N} - u_n| for the *unwrapped* analytic wave (A=1). Zero iff k N = 2 pi p.
  function seamMismatch(N, k, t) {
    return Math.abs(Math.cos(k * N - omega(k) * t) - Math.cos(-omega(k) * t));
  }
  function seamAmplitude(N, k) { return 2 * Math.abs(Math.sin(k * N / 2)); } // max over t
  function wavesAround(N, k) { return k * N / TWO_PI; } // p' = kL/2pi

  // Number of allowed k in [k0, k0+dk)
  function windowCount(N, bc, k0, dk) {
    var ks = allowedK(N, bc), c = 0, i, e = 1e-12;
    for (i = 0; i < ks.length; i++) if (ks[i].k >= k0 - e && ks[i].k < k0 + dk - e) c++;
    return c;
  }
  function nearest(N, bc, k) {
    var ks = allowedK(N, bc), best = 0, d = 1e9, i;
    var kk = (bc === 'box') ? Math.abs(k) : k;
    for (i = 0; i < ks.length; i++) { var dd = Math.abs(ks[i].k - kk); if (dd < d) { d = dd; best = i; } }
    return ks[best];
  }
  // (1/N) sum_k omega_k  (box: mean over the N-1 modes); integral limit is 4/pi
  function meanOmega(N, bc) {
    var ks = allowedK(N, bc), s = 0, i;
    for (i = 0; i < ks.length; i++) s += omega(ks[i].k);
    return s / ks.length;
  }
  var MEAN_OMEGA_INF = 4 / Math.PI;

  // Velocity-Verlet evolution of the ring, initial data = mode k at t=0. Independent check of the
  // analytic mode and of energy conservation.
  function verletRing(N, k, A, dt, steps) {
    var u = [], v = [], a = [], n, s, w = omega(k);
    for (n = 0; n < N; n++) { u.push(A * Math.cos(k * n)); v.push(A * w * Math.sin(k * n)); a.push(0); }
    function acc(u, a) { for (var i = 0; i < N; i++) a[i] = u[(i + 1) % N] - 2 * u[i] + u[(i - 1 + N) % N]; }
    function energy() { var e = 0, i; for (i = 0; i < N; i++) { var d = u[(i + 1) % N] - u[i]; e += 0.5 * v[i] * v[i] + 0.5 * d * d; } return e; }
    acc(u, a);
    var E0 = energy(), maxDrift = 0;
    for (s = 0; s < steps; s++) {
      for (n = 0; n < N; n++) { v[n] += 0.5 * dt * a[n]; u[n] += dt * v[n]; }
      acc(u, a);
      for (n = 0; n < N; n++) v[n] += 0.5 * dt * a[n];
      if (E0 > 0) maxDrift = Math.max(maxDrift, Math.abs(energy() - E0) / E0);
    }
    var T = steps * dt, err = 0;
    for (n = 0; n < N; n++) err = Math.max(err, Math.abs(u[n] - A * Math.cos(k * n - w * T)));
    return { posErr: err, drift: maxDrift };
  }

  // Thermal-like mixture: every allowed mode present with random phase, |c_k| ~ 1/omega_k.
  function rng(seed) { var a = seed >>> 0; return function () { a |= 0; a = a + 0x6D2B79F5 | 0; var t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function gauss(r) { return Math.sqrt(-2 * Math.log(1 - r())) * Math.cos(TWO_PI * r()); }
  function thermal(N, bc, seed, rms) {
    var r = rng(seed), ks = allowedK(N, bc), modes = [], i;
    for (i = 0; i < ks.length; i++) {
      var w = omega(ks[i].k);
      if (w < 1e-9) continue;               // k=0 is a rigid translation, not a vibration
      modes.push({ p: ks[i].p, k: ks[i].k, omega: w, cr: gauss(r) / w, ci: gauss(r) / w });
    }
    var m = { modes: modes, bc: bc, N: N, scale: 1 };
    var s = 0, n; for (n = 0; n < N; n++) { var u = thermalU(m, n, 0); s += u * u; }
    var cur = Math.sqrt(s / N) || 1; m.scale = (rms || 0.25) / cur;
    return m;
  }
  function thermalU(m, n, t) {
    var s = 0, i, md, ph;
    for (i = 0; i < m.modes.length; i++) {
      md = m.modes[i];
      if (m.bc === 'box') s += Math.sin(md.k * n) * (md.cr * Math.cos(md.omega * t) - md.ci * Math.sin(md.omega * t));
      else { ph = md.k * n - md.omega * t; s += md.cr * Math.cos(ph) - md.ci * Math.sin(ph); }
    }
    return s * m.scale;
  }

  return { omega: omega, pRange: pRange, allowedK: allowedK, spacing: spacing, density: density, modeU: modeU,
    eomResidual: eomResidual, seamMismatch: seamMismatch, seamAmplitude: seamAmplitude, wavesAround: wavesAround,
    windowCount: windowCount, nearest: nearest, meanOmega: meanOmega, MEAN_OMEGA_INF: MEAN_OMEGA_INF,
    verletRing: verletRing, thermal: thermal, thermalU: thermalU, TWO_PI: TWO_PI };
})();
if (typeof module !== 'undefined') module.exports = Core;
