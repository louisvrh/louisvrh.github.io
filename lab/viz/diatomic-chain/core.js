/* Diatomic chain core. No DOM. Units: a = 1 (cell length), atom 1 at x = n, atom 2 at x = n + 1/2.
   Equations of motion (nearest-neighbour springs):
     m1 u1_n'' = C1 (u2_n - u1_n) + C2 (u2_{n-1} - u1_n)
     m2 u2_n'' = C2 (u1_{n+1} - u2_n) + C1 (u1_n - u2_n)
   Trial u_n^(j) = A_j exp(i(k n - w t)):
     (S - m1 w^2) A1 = g  A2,   (S - m2 w^2) A2 = g* A1,   S = C1 + C2,  g = C1 + C2 exp(-ik)
   => m1 m2 w^4 - S (m1+m2) w^2 + 4 C1 C2 sin^2(k/2) = 0.
   Complex numbers are [re, im] pairs. */
const Core = (function () {
  'use strict';
  const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const cabs = a => Math.hypot(a[0], a[1]);
  const cconj = a => [a[0], -a[1]];
  const cdiv = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
  const cexp = th => [Math.cos(th), Math.sin(th)];

  // p = {m1, m2, C1, C2}; k in units of 1/a. Returns [w^2_acoustic, w^2_optical].
  function omega2(k, p) {
    const S = p.C1 + p.C2, M = p.m1 + p.m2, s2 = Math.sin(k / 2) * Math.sin(k / 2);
    const disc = Math.max(0, S * S * M * M - 16 * p.m1 * p.m2 * p.C1 * p.C2 * s2);
    const rt = Math.sqrt(disc);
    const wOpt = (S * M + rt) / (2 * p.m1 * p.m2);
    // stable form of the small root (no cancellation as k -> 0)
    const wAc = 8 * p.C1 * p.C2 * s2 / (S * M + rt);
    return [wAc, wOpt];
  }
  const omega = (k, p, branch) => Math.sqrt(omega2(k, p)[branch]); // branch 0 acoustic, 1 optical

  // Eigenvector, normalised so max(|A1|,|A2|) = 1. Gauge: dominant-at-edge atom real (A1 for acoustic, A2 for optical).
  function eigvec(k, p, branch) {
    const w2 = omega2(k, p)[branch], S = p.C1 + p.C2;
    const g = [p.C1 + p.C2 * Math.cos(k), -p.C2 * Math.sin(k)];  // C1 + C2 e^{-ik}
    const d1 = S - p.m1 * w2, d2 = S - p.m2 * w2;
    let A1, A2;
    if (Math.abs(d1) >= Math.abs(d2)) { A1 = g; A2 = [d1, 0]; }      // from eq. 1: A1 d1 = g A2
    else { A1 = [d2, 0]; A2 = cconj(g); }                             // from eq. 2: A2 d2 = g* A1
    if (cabs(A1) < 1e-14 && cabs(A2) < 1e-14) { A1 = [1, 0]; A2 = [1, 0]; }
    const ref = branch === 0 ? A1 : A2;
    const rot = cconj(ref); const rm = cabs(ref);
    if (rm > 1e-14) { const u = [rot[0] / rm, rot[1] / rm]; A1 = cmul(A1, u); A2 = cmul(A2, u); }
    const mx = Math.max(cabs(A1), cabs(A2));
    return { A1: [A1[0] / mx, A1[1] / mx], A2: [A2[0] / mx, A2[1] / mx], w: Math.sqrt(w2), w2 };
  }

  // Intrinsic ratio B1/B2 = (A1/A2) e^{ik/2}: removes the e^{ik/2} lag expected for atoms a/2 apart.
  function ratioB(ev, k) {
    if (cabs(ev.A2) < 1e-12) return null;
    return cmul(cdiv(ev.A1, ev.A2), cexp(k / 2));
  }

  // Displacement of atom j (1|2) in cell n at time t for eigenvector ev, wavenumber k (physical amplitude 1).
  function dispPhase(ev, k, j, n, theta) { // theta = omega t (or any display phase)
    const A = j === 1 ? ev.A1 : ev.A2, th = k * n - theta;
    return A[0] * Math.cos(th) - A[1] * Math.sin(th);   // Re[A e^{i th}]
  }
  function disp(ev, k, j, n, t) { return dispPhase(ev, k, j, n, ev.w * t); }

  // Observables for a mode (all in units of max atom amplitude)
  function observables(ev, k, p) {
    const intra = cabs([ev.A2[0] - ev.A1[0], ev.A2[1] - ev.A1[1]]);        // |u2_n - u1_n| amplitude
    const e = cexp(k), A1e = cmul(ev.A1, e);
    const inter = cabs([A1e[0] - ev.A2[0], A1e[1] - ev.A2[1]]);            // |u1_{n+1} - u2_n| amplitude
    const M = p.m1 + p.m2;
    const com = cabs([(p.m1 * ev.A1[0] + p.m2 * ev.A2[0]) / M, (p.m1 * ev.A1[1] + p.m2 * ev.A2[1]) / M]);
    return { intra, inter, com };
  }

  function soundSpeed(p) { return Math.sqrt(p.C1 * p.C2 / ((p.C1 + p.C2) * (p.m1 + p.m2))); }

  // Residual of the equations of motion for a mode, evaluated on an explicit chain (cells n-1, n, n+1) at time t.
  function residual(k, p, branch, t, n) {
    const ev = eigvec(k, p, branch), w2 = ev.w2, d = (j, m) => disp(ev, k, j, m, t);
    const acc = (j, m) => -w2 * d(j, m); // u'' of the trial solution
    const r1 = p.m1 * acc(1, n) - (p.C1 * (d(2, n) - d(1, n)) + p.C2 * (d(2, n - 1) - d(1, n)));
    const r2 = p.m2 * acc(2, n) - (p.C2 * (d(1, n + 1) - d(2, n)) + p.C1 * (d(1, n) - d(2, n)));
    return Math.max(Math.abs(r1), Math.abs(r2));
  }

  // ---- Independent check tools: brute-force periodic chain of N cells ----
  function accel(u, N, p, out) { // u: Float64Array length 2N, index 2n = atom1_n, 2n+1 = atom2_n
    for (let n = 0; n < N; n++) {
      const u1 = u[2 * n], u2 = u[2 * n + 1], u2m = u[2 * ((n - 1 + N) % N) + 1], u1p = u[2 * ((n + 1) % N)];
      out[2 * n] = (p.C1 * (u2 - u1) + p.C2 * (u2m - u1)) / p.m1;
      out[2 * n + 1] = (p.C2 * (u1p - u2) + p.C1 * (u1 - u2)) / p.m2;
    }
  }
  function energy(u, v, N, p) {
    let E = 0;
    for (let n = 0; n < N; n++) {
      E += 0.5 * p.m1 * v[2 * n] ** 2 + 0.5 * p.m2 * v[2 * n + 1] ** 2;
      E += 0.5 * p.C1 * (u[2 * n + 1] - u[2 * n]) ** 2;
      E += 0.5 * p.C2 * (u[2 * ((n + 1) % N)] - u[2 * n + 1]) ** 2;
    }
    return E;
  }
  // Velocity-Verlet integration of the full chain (symplectic).
  function verlet(u0, v0, N, p, dt, steps) {
    const u = Float64Array.from(u0), v = Float64Array.from(v0), a = new Float64Array(2 * N), a2 = new Float64Array(2 * N);
    accel(u, N, p, a);
    for (let s = 0; s < steps; s++) {
      for (let i = 0; i < 2 * N; i++) u[i] += v[i] * dt + 0.5 * a[i] * dt * dt;
      accel(u, N, p, a2);
      for (let i = 0; i < 2 * N; i++) { v[i] += 0.5 * (a[i] + a2[i]) * dt; a[i] = a2[i]; }
    }
    return { u, v };
  }
  // Jacobi eigenvalues of a real symmetric matrix (array of arrays)
  function jacobiEig(A) {
    const n = A.length, M = A.map(r => r.slice());
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += M[i][j] * M[i][j];
      if (off < 1e-26) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(M[p][q]) < 1e-300) continue;
        const th = (M[q][q] - M[p][p]) / (2 * M[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const kp = M[k][p], kq = M[k][q]; M[k][p] = c * kp - s * kq; M[k][q] = s * kp + c * kq; }
        for (let k = 0; k < n; k++) { const pk = M[p][k], qk = M[q][k]; M[p][k] = c * pk - s * qk; M[q][k] = s * pk + c * qk; }
      }
    }
    return M.map((r, i) => r[i]).sort((x, y) => x - y);
  }
  // Mass-weighted dynamical matrix of the periodic N-cell chain
  function dynMatrix(N, p) {
    const D = Array.from({ length: 2 * N }, () => new Array(2 * N).fill(0)), m = i => (i % 2 === 0 ? p.m1 : p.m2);
    const K = Array.from({ length: 2 * N }, () => new Array(2 * N).fill(0));
    const spring = (i, j, c) => { K[i][i] += c; K[j][j] += c; K[i][j] -= c; K[j][i] -= c; };
    for (let n = 0; n < N; n++) { spring(2 * n, 2 * n + 1, p.C1); spring(2 * n + 1, 2 * ((n + 1) % N), p.C2); }
    for (let i = 0; i < 2 * N; i++) for (let j = 0; j < 2 * N; j++) D[i][j] = K[i][j] / Math.sqrt(m(i) * m(j));
    return D;
  }

  return { omega2, omega, eigvec, ratioB, disp, dispPhase, observables, soundSpeed, residual, accel, energy, verlet, jacobiEig, dynMatrix, cmul, cabs, cdiv };
})();
if (typeof module !== 'undefined') module.exports = Core;
