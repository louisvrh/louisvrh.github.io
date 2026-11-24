/* Monatomic chain core: units a = m = kappa = 1  (omega_D = 1, time in 1/omega_D).
   u'' _n = u_{n+1} - 2u_n + u_{n-1}.  BCs: 'fixed' (walls at both ends), 'free', 'periodic'. */
var Core = (function () {
  'use strict';

  function uAt(u, i, bc) {
    var N = u.length;
    if (bc === 'periodic') return u[((i % N) + N) % N];
    if (bc === 'fixed') return (i < 0 || i >= N) ? 0 : u[i];
    return u[i < 0 ? 0 : (i >= N ? N - 1 : i)]; // free: zero-force ghost (u_{-1}=u_0)
  }

  /* F_n = kappa (u_{n+1} - 2 u_n + u_{n-1}) */
  function force(u, bc, F) {
    var N = u.length;
    for (var n = 0; n < N; n++) F[n] = uAt(u, n + 1, bc) - 2 * u[n] + uAt(u, n - 1, bc);
    return F;
  }

  /* Springs as bonds (a,b); a=-1 is the left wall, b=N the right wall (fixed BC). Extension = u_b - u_a. */
  function bonds(N, bc) {
    var B = [], i;
    if (bc === 'fixed') { for (i = 0; i <= N; i++) B.push([i - 1, i]); }
    else { for (i = 0; i < N - 1; i++) B.push([i, i + 1]); if (bc === 'periodic') B.push([N - 1, 0]); }
    return B;
  }
  function extensions(u, bc) {
    var B = bonds(u.length, bc), e = new Array(B.length);
    for (var i = 0; i < B.length; i++) e[i] = uAt(u, B[i][1], bc) - uAt(u, B[i][0], bc);
    return e;
  }
  function energy(u, v, bc) {
    var K = 0, P = 0, i, e = extensions(u, bc);
    for (i = 0; i < v.length; i++) K += 0.5 * v[i] * v[i];
    for (i = 0; i < e.length; i++) P += 0.5 * e[i] * e[i];
    return { K: K, P: P, E: K + P };
  }

  /* velocity Verlet (symplectic); optional non-physical damping via exact exponential half-kicks */
  function Sim(N, bc) {
    this.N = N; this.bc = bc; this.t = 0; this.gamma = 0;
    this.u = new Float64Array(N); this.v = new Float64Array(N);
    this.F = new Float64Array(N); this.u0 = null; this.v0 = null;
  }
  Sim.prototype.setState = function (u, v) {
    this.u.set(u); this.v.set(v); this.t = 0;
    this.u0 = Float64Array.from(u); this.v0 = Float64Array.from(v);
    force(this.u, this.bc, this.F);
  };
  Sim.prototype.step = function (dt) {
    var N = this.N, u = this.u, v = this.v, F = this.F, n;
    var d = this.gamma > 0 ? Math.exp(-0.5 * this.gamma * dt) : 1;
    if (d !== 1) for (n = 0; n < N; n++) v[n] *= d;
    for (n = 0; n < N; n++) v[n] += 0.5 * dt * F[n];
    for (n = 0; n < N; n++) u[n] += dt * v[n];
    force(u, this.bc, F);
    for (n = 0; n < N; n++) v[n] += 0.5 * dt * F[n];
    if (d !== 1) for (n = 0; n < N; n++) v[n] *= d;
    this.t += dt;
  };
  Sim.prototype.energy = function () { return energy(this.u, this.v, this.bc); };

  /* ---- exact normal modes: orthonormal real basis, omega_j = 2 sin(k/2) ---- */
  function modes(N, bc) {
    var M = [], j, n, vec, c;
    if (bc === 'fixed') {
      c = Math.sqrt(2 / (N + 1));
      for (j = 1; j <= N; j++) {
        vec = new Float64Array(N);
        for (n = 0; n < N; n++) vec[n] = c * Math.sin(Math.PI * j * (n + 1) / (N + 1));
        var k = Math.PI * j / (N + 1);
        M.push({ j: j, k: k, omega: 2 * Math.sin(k / 2), vec: vec });
      }
    } else if (bc === 'free') {
      for (j = 0; j < N; j++) {
        vec = new Float64Array(N);
        c = j === 0 ? Math.sqrt(1 / N) : Math.sqrt(2 / N);
        for (n = 0; n < N; n++) vec[n] = c * Math.cos(Math.PI * j * (n + 0.5) / N);
        var kf = Math.PI * j / N;
        M.push({ j: j, k: kf, omega: 2 * Math.sin(kf / 2), vec: vec });
      }
    } else {
      for (j = 0; j <= Math.floor(N / 2); j++) {
        var kp = 2 * Math.PI * j / N, om = 2 * Math.sin(kp / 2);
        var vc = new Float64Array(N), vs = null;
        if (j === 0) { for (n = 0; n < N; n++) vc[n] = Math.sqrt(1 / N); }
        else if (2 * j === N) { for (n = 0; n < N; n++) vc[n] = (n % 2 ? -1 : 1) * Math.sqrt(1 / N); }
        else {
          vs = new Float64Array(N);
          for (n = 0; n < N; n++) { vc[n] = Math.sqrt(2 / N) * Math.cos(kp * n); vs[n] = Math.sqrt(2 / N) * Math.sin(kp * n); }
        }
        M.push({ j: j, k: kp, omega: om, vec: vc });
        if (vs) M.push({ j: j, k: kp, omega: om, vec: vs });
      }
    }
    return M;
  }
  function dot(a, b) { var s = 0; for (var i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

  /* Exact solution u(t) from (u0, v0): each mode c cos(wt) + (d/w) sin(wt); w = 0: c + d t */
  function exact(N, bc, u0, v0, t) {
    var M = modes(N, bc), out = new Float64Array(N), m, i, c, d, a;
    for (m = 0; m < M.length; m++) {
      c = dot(M[m].vec, u0); d = dot(M[m].vec, v0); var w = M[m].omega;
      a = w > 1e-12 ? c * Math.cos(w * t) + d / w * Math.sin(w * t) : c + d * t;
      for (i = 0; i < N; i++) out[i] += a * M[m].vec[i];
    }
    return out;
  }

  /* Energy content per wavevector k (merging cos/sin pairs): [{k, omega, E}] with k in [0, pi] (a = 1) */
  function spectrum(u, v, bc) {
    var N = u.length, M = modes(N, bc), res = [], m, last = null, c, d, E;
    for (m = 0; m < M.length; m++) {
      c = dot(M[m].vec, u); d = dot(M[m].vec, v);
      E = 0.5 * (d * d + M[m].omega * M[m].omega * c * c);
      if (last && bc === 'periodic' && last.j === M[m].j) { last.E += E; }
      else { last = { j: M[m].j, k: M[m].k, omega: M[m].omega, E: E }; res.push(last); }
    }
    return res;
  }

  /* ---- presets: return {u, v} ---- */
  function preset(name, N, bc, amp, opt) {
    opt = opt || {};
    var u = new Float64Array(N), v = new Float64Array(N), n, j, M, k;
    if (name === 'pluck') { u[Math.floor(N / 2)] = amp; }
    else if (name === 'uniform') { for (n = 0; n < N; n++) u[n] = amp; }
    else if (name === 'mode' || name === 'longwave' || name === 'zoneedge') {
      var jmax = bc === 'fixed' ? N : (bc === 'free' ? N - 1 : Math.floor(N / 2));
      j = name === 'longwave' ? 1 : (name === 'zoneedge' ? jmax : Math.max(1, Math.min(jmax, opt.j || 1)));
      M = modes(N, bc);
      var sel = M.filter(function (x) { return x.j === j; });
      var vec = sel[0].vec, mx = 0;
      for (n = 0; n < N; n++) mx = Math.max(mx, Math.abs(vec[n]));
      if (bc === 'periodic' && sel.length === 2 && opt.travel) {
        // right-moving wave u = A cos(k n - w t): u = A cos(kn), v = A w sin(kn)
        k = sel[0].k;
        for (n = 0; n < N; n++) { u[n] = amp * Math.cos(k * n); v[n] = amp * sel[0].omega * Math.sin(k * n); }
      } else {
        for (n = 0; n < N; n++) u[n] = amp * vec[n] / mx;
      }
    }
    return { u: u, v: v };
  }

  /* Bessel J_m(x) by quadrature (benchmark: delta pluck on an infinite chain gives u_n = J_{2n}(2t)) */
  function besselJ(m, x) {
    var S = 4000, s = 0, i, tau;
    for (i = 0; i < S; i++) { tau = (i + 0.5) * Math.PI / S; s += Math.cos(m * tau - x * Math.sin(tau)); }
    return s / S;
  }

  return { uAt: uAt, force: force, bonds: bonds, extensions: extensions, energy: energy, Sim: Sim,
           modes: modes, exact: exact, spectrum: spectrum, preset: preset, besselJ: besselJ, DT: 0.01 };
})();
if (typeof module !== 'undefined') module.exports = Core;
