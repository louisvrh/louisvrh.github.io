// Physics core: monatomic chain with periodic BCs, phonon (Bose-Einstein) occupations.
// Units: hbar = k_B = a = m = kappa = omega_D = 1.  Displacements u_n in units of a.
const Core = (function () {
  'use strict';

  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Eq. 4.10:  omega_k = 2 omega_D |sin(ka/2)|
  function omega(k) { return 2 * Math.abs(Math.sin(k / 2)); }

  // Eq. 5.2: k = p 2pi/(N a), N values in the first Brillouin zone (-pi, pi]
  function allowedK(N) {
    const pmin = Math.floor(-N / 2) + 1, pmax = Math.floor(N / 2);
    const p = [], k = [], w = [];
    for (let q = pmin; q <= pmax; q++) {
      const kk = 2 * Math.PI * q / N;
      p.push(q); k.push(kk); w.push(omega(kk));
    }
    return { N: N, p: p, k: k, w: w, pmin: pmin, pmax: pmax };
  }

  // Eq. 5.12: Bose-Einstein, no chemical potential
  function nB(w, T) {
    if (w <= 0) return Infinity;
    const x = w / T;
    if (x > 700) return 0;
    return 1 / Math.expm1(x);
  }

  // Sample n from P(n) = (1-x) x^n, x = exp(-w/T) (geometric = BE for one mode). U in (0,1].
  function sampleBE(w, T, U) {
    if (w <= 0) return 0;
    const lx = -w / T;
    if (lx < -700) return 0;
    return Math.floor(Math.log(U) / lx);
  }

  // Classical amplitude of a travelling wave u_n = A cos(kna - wt + phi) carrying energy hbar w (n + 1/2)
  // (E = (1/2) m w^2 A^2 N)  ->  A = sqrt(2 hbar (n+1/2) / (N m w))
  function amp(n, w, N, zp) { return Math.sqrt(2 * (n + (zp ? 0.5 : 0)) / (N * w)); }

  function makeModes(al, nArr, phi, zp, only) {
    const m = [];
    for (let i = 0; i < al.k.length; i++) {
      if (al.w[i] <= 1e-12) continue;           // k = 0: rigid translation, omega = 0, not a vibration
      if (only !== null && only !== undefined && i !== only) continue;
      const A = amp(nArr[i], al.w[i], al.N, zp);
      if (A > 0) m.push({ i: i, k: al.k[i], w: al.w[i], A: A, phi: phi[i] });
    }
    return m;
  }

  // Eqs. 4.5-4.6: u(x,t) = sum_k A_k cos(k x - w_k t + phi_k)   (x in units of a, integer for atoms)
  function field(modes, x, t) {
    let u = 0;
    for (let j = 0; j < modes.length; j++) { const m = modes[j]; u += m.A * Math.cos(m.k * x - m.w * t + m.phi); }
    return u;
  }
  function velocity(modes, x, t) {
    let v = 0;
    for (let j = 0; j < modes.length; j++) { const m = modes[j]; v += m.A * m.w * Math.sin(m.k * x - m.w * t + m.phi); }
    return v;
  }
  function accel(modes, x, t) {
    let a = 0;
    for (let j = 0; j < modes.length; j++) { const m = modes[j]; a -= m.w * m.w * m.A * Math.cos(m.k * x - m.w * t + m.phi); }
    return a;
  }

  function atoms(modes, N, t) {
    const u = new Float64Array(N), v = new Float64Array(N), a = new Float64Array(N);
    for (let n = 0; n < N; n++) { u[n] = field(modes, n, t); v[n] = velocity(modes, n, t); a[n] = accel(modes, n, t); }
    return { u: u, v: v, a: a };
  }

  // max_n | m u_n'' - kappa (u_{n+1} - 2 u_n + u_{n-1}) |, ring indices
  function eomResidual(modes, N, t) {
    const s = atoms(modes, N, t);
    let r = 0;
    for (let n = 0; n < N; n++) {
      const lap = s.u[(n + 1) % N] - 2 * s.u[n] + s.u[(n - 1 + N) % N];
      r = Math.max(r, Math.abs(s.a[n] - lap));
    }
    return r;
  }

  // Classical energy of the chain: kinetic + elastic
  function energy(modes, N, t) {
    const s = atoms(modes, N, t);
    let ke = 0, pe = 0;
    for (let n = 0; n < N; n++) {
      ke += 0.5 * s.v[n] * s.v[n];
      const d = s.u[(n + 1) % N] - s.u[n];
      pe += 0.5 * d * d;
    }
    return ke + pe;
  }

  // sum_k hbar w_k (n_k + 1/2 zp) over nonzero modes
  function quantumEnergy(al, nArr, zp) {
    let E = 0;
    for (let i = 0; i < al.k.length; i++) { if (al.w[i] <= 1e-12) continue; E += al.w[i] * (nArr[i] + (zp ? 0.5 : 0)); }
    return E;
  }

  // Eq. 5.13: Delta U = sum_k n_B hbar w_k   (k = 0 translation mode omitted)
  function deltaU(N, T) {
    const al = allowedK(N); let U = 0;
    for (let i = 0; i < N; i++) if (al.w[i] > 1e-12) U += nB(al.w[i], T) * al.w[i];
    return U;
  }
  function meanCount(N, T) {
    const al = allowedK(N); let c = 0;
    for (let i = 0; i < N; i++) if (al.w[i] > 1e-12) c += nB(al.w[i], T);
    return c;
  }

  // Velocity-Verlet (kick-drift-kick) for the ring, for cross-checking the analytic mode solution
  function verlet(u0, v0, dt, steps) {
    const N = u0.length;
    let u = Float64Array.from(u0), v = Float64Array.from(v0);
    const acc = (uu) => { const a = new Float64Array(N); for (let n = 0; n < N; n++) a[n] = uu[(n + 1) % N] - 2 * uu[n] + uu[(n - 1 + N) % N]; return a; };
    let a = acc(u);
    for (let s = 0; s < steps; s++) {
      for (let n = 0; n < N; n++) v[n] += 0.5 * dt * a[n];
      for (let n = 0; n < N; n++) u[n] += dt * v[n];
      a = acc(u);
      for (let n = 0; n < N; n++) v[n] += 0.5 * dt * a[n];
    }
    return { u: u, v: v };
  }

  function regime(T) {
    if (T < 0.3) return 'low T: k_BT << ħω_max, only small-k (low-ω) modes are occupied';
    if (T < 2) return 'crossover: modes with ħω ≲ k_BT are occupied, the rest are frozen out';
    return 'high T: k_BT ≳ ħω_max = 2ħω_D, every mode carries ≈ k_BT (equipartition)';
  }

  return { mulberry32, omega, allowedK, nB, sampleBE, amp, makeModes, field, velocity, accel, atoms,
           eomResidual, energy, quantumEnergy, deltaU, meanCount, verlet, regime };
})();
if (typeof module !== 'undefined') module.exports = Core;
