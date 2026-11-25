/* Square-lattice scalar mass-and-spring model, units a = omega_D = 1.
   u'' = omega_D^2 (u[i+1,j]+u[i-1,j]+u[i,j+1]+u[i,j-1]-4u[i,j])  ->  Eq. 4.12 */
const Core = (function () {
  const PI = Math.PI, SQ2 = Math.SQRT2;
  const s2h = x => { const s = Math.sin(x / 2); return s * s; };

  // Eq. 4.12: omega^2 = 4[sin^2(kx/2)+sin^2(ky/2)]
  const omega2 = (kx, ky) => 4 * (s2h(kx) + s2h(ky));
  const omega = (kx, ky) => Math.sqrt(omega2(kx, ky));

  // v_g = grad_k omega, analytic: d omega/d kx = sin(kx)/omega
  function grad(kx, ky) {
    const w = omega(kx, ky);
    if (w < 1e-9) { const K = Math.hypot(kx, ky); return K < 1e-12 ? [0, 0] : [kx / K, ky / K]; }
    return [Math.sin(kx) / w, Math.sin(ky) / w];
  }

  const kvec = (K, th) => [K * Math.cos(th), K * Math.sin(th)];
  const omegaDir = (K, th) => { const k = kvec(K, th); return omega(k[0], k[1]); };
  const phaseSpeed = (K, th) => omegaDir(K, th) / K;           // omega/|k|, units omega_D a
  const ratio1110 = K => omega(K / SQ2, K / SQ2) / omega(K, 0);  // omega([11])/omega([10]) at equal |k|
  // largest |k| in direction th (rad) that stays inside the square FBZ
  const kMax = th => PI / Math.max(Math.abs(Math.cos(th)), Math.abs(Math.sin(th)));
  // small-k expansion: (v_p / omega_D a)^2 = 1 - K^2 (cos^4+sin^4)/12 + O(K^4)
  const vpSmallK = (K, th) => Math.sqrt(1 - K * K * (Math.pow(Math.cos(th), 4) + Math.pow(Math.sin(th), 4)) / 12);

  // plane-wave mode u(i,j,t) = cos(kx i + ky j - omega t)
  function mode(kx, ky, i, j, t) { return Math.cos(kx * i + ky * j - omega(kx, ky) * t); }
  // residual of the scalar equation of motion for the mode (finite-difference u_tt)
  function modeResidualFD(kx, ky, i, j, t, h) {
    h = h || 1e-3;
    const f = (a, b, tt) => mode(kx, ky, a, b, tt);
    const utt = (f(i, j, t + h) - 2 * f(i, j, t) + f(i, j, t - h)) / (h * h);
    const rhs = f(i + 1, j, t) + f(i - 1, j, t) + f(i, j + 1, t) + f(i, j - 1, t) - 4 * f(i, j, t);
    return utt - rhs;
  }
  // residual using analytic u_tt = -omega^2 u
  function modeResidualExact(kx, ky, i, j, t) {
    const w2 = omega2(kx, ky), f = (a, b) => mode(kx, ky, a, b, t);
    return -w2 * f(i, j) - (f(i + 1, j) + f(i - 1, j) + f(i, j + 1) + f(i, j - 1) - 4 * f(i, j));
  }

  // periodic N x N lattice integrated with velocity Verlet (used for testing only)
  function makeLattice(N) { return { N, u: new Float64Array(N * N), v: new Float64Array(N * N), a: new Float64Array(N * N) }; }
  function accel(L) {
    const N = L.N, u = L.u, a = L.a;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const ip = (i + 1) % N, im = (i + N - 1) % N, jp = (j + 1) % N, jm = (j + N - 1) % N;
      a[j * N + i] = u[j * N + ip] + u[j * N + im] + u[jp * N + i] + u[jm * N + i] - 4 * u[j * N + i];
    }
  }
  function initMode(L, kx, ky) {
    const w = omega(kx, ky);
    for (let j = 0; j < L.N; j++) for (let i = 0; i < L.N; i++) {
      const ph = kx * i + ky * j; L.u[j * L.N + i] = Math.cos(ph); L.v[j * L.N + i] = w * Math.sin(ph);
    }
    accel(L);
  }
  function verlet(L, dt) {
    const n = L.u.length;
    for (let m = 0; m < n; m++) { L.v[m] += 0.5 * dt * L.a[m]; L.u[m] += dt * L.v[m]; }
    accel(L);
    for (let m = 0; m < n; m++) L.v[m] += 0.5 * dt * L.a[m];
  }
  function energy(L) {
    const N = L.N; let E = 0;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const m = j * N + i, du1 = L.u[m] - L.u[j * N + (i + 1) % N], du2 = L.u[m] - L.u[((j + 1) % N) * N + i];
      E += 0.5 * L.v[m] * L.v[m] + 0.5 * du1 * du1 + 0.5 * du2 * du2;
    }
    return E;
  }
  return { omega2, omega, grad, kvec, omegaDir, phaseSpeed, ratio1110, kMax, vpSmallK, mode,
           modeResidualFD, modeResidualExact, makeLattice, initMode, verlet, energy, accel };
})();
if (typeof module !== 'undefined') module.exports = Core;
