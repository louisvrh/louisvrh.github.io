// Coupled Quantum Wells — faithful JS port of the PHY60003 notebook (R. Oulton).
// Complex numbers are [re, im]; 2x2 matrices are [a, b, c, d] for [[a,b],[c,d]].
const Core = (() => {
  const m = 9.11e-31, hbar = 6.63e-34 / (2 * Math.PI), e = 1.6e-19;
  const K = Math.sqrt(2 * e * m) / hbar;

  const cmul = (x, y) => [x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]];
  const cadd = (x, y) => [x[0] + y[0], x[1] + y[1]];
  const csub = (x, y) => [x[0] - y[0], x[1] - y[1]];
  const cdiv = (x, y) => {
    const den = y[0] * y[0] + y[1] * y[1];
    return [(x[0] * y[0] + x[1] * y[1]) / den, (x[1] * y[0] - x[0] * y[1]) / den];
  };
  const cabs = (x) => Math.hypot(x[0], x[1]);
  const cexp = (x) => { const r = Math.exp(x[0]); return [r * Math.cos(x[1]), r * Math.sin(x[1])]; };

  // M2 * M1 (matrix product, M2 on the left, as np.matmul(M2, M1))
  function mm(A, B) {
    return [
      cadd(cmul(A[0], B[0]), cmul(A[1], B[2])),
      cadd(cmul(A[0], B[1]), cmul(A[1], B[3])),
      cadd(cmul(A[2], B[0]), cmul(A[3], B[2])),
      cadd(cmul(A[2], B[1]), cmul(A[3], B[3])),
    ];
  }
  const I2 = () => [[1, 0], [0, 0], [0, 0], [1, 0]];

  // generate_struc(N, V0, dW, delta, dB)  -> V [eV], d [m]
  function generateStruc(N, V0, dW, delta, dB) {
    const n = 2 * (N + 1) + 1;
    const V = new Array(n).fill(0), d = new Array(n).fill(0);
    for (let i = 1; i < 2 * (N + 1); i += 2) V[i] = V0;
    for (let i = 0; i < n; i += 2) d[i] = dW * 1e-9;
    for (let i = 0; i < n; i += 4) d[i] = (dW + delta) * 1e-9;
    for (let i = 1; i < 2 * (N + 1); i += 2) d[i] = dB * 1e-9;
    d[0] = 5e-9; d[2 * (N + 1)] = 5e-9;
    return { V, d };
  }

  // TM_T: interface matrix
  function TM_T(E, V1, V2) {
    const c1 = E >= V1 ? [Math.sqrt(E - V1), 0] : [0, Math.sqrt(V1 - E)];
    const c2 = E >= V2 ? [Math.sqrt(E - V2), 0] : [0, Math.sqrt(V2 - E)];
    const c = cdiv(c1, c2);
    const p = [0.5 * (1 + c[0]), 0.5 * c[1]], q = [0.5 * (1 - c[0]), -0.5 * c[1]];
    return [p, q, q, p];
  }

  // TM_P: propagation matrix
  function TM_P(E, V, d) {
    const kd = E >= V ? [0, d * K * Math.sqrt(E - V)] : [-d * K * Math.sqrt(V - E), 0];
    return [cexp(kd), [0, 0], [0, 0], cexp([-kd[0], -kd[1]])];
  }

  // numpy.linalg.det: LU with partial pivoting (pivot by |re|+|im|), then sign*exp(logdet)
  function detLU(T) {
    let a00 = T[0], a01 = T[1], a10 = T[2], a11 = T[3], sgn = 1;
    if (Math.abs(a10[0]) + Math.abs(a10[1]) > Math.abs(a00[0]) + Math.abs(a00[1])) {
      [a00, a10] = [a10, a00]; [a01, a11] = [a11, a01]; sgn = -1;
    }
    let u11 = a11;
    if (a00[0] !== 0 || a00[1] !== 0) {
      const l = cmul(a10, cdiv([1, 0], a00));
      u11 = csub(a11, cmul(l, a01));
    }
    const m0 = Math.hypot(a00[0], a00[1]), m1 = Math.hypot(u11[0], u11[1]);
    if (m0 === 0 || m1 === 0) return [0, 0];
    const s = cmul([a00[0] / m0, a00[1] / m0], [u11[0] / m1, u11[1] / m1]);
    const mag = Math.exp(Math.log(m0) + Math.log(m1));
    return [sgn * s[0] * mag, sgn * s[1] * mag];
  }

  function totalT(E, V, d) {
    let T = I2();
    for (let j = 0; j < V.length - 1; j++) {
      if (j === 0) T = mm(TM_P(E, V[0], d[0]), T);
      T = mm(TM_T(E, V[j], V[j + 1]), T);
      T = mm(TM_P(E, V[j + 1], d[j + 1]), T);
    }
    return T;
  }

  function linspace(a, b, n) {
    const out = new Float64Array(n);
    if (n === 1) { out[0] = a; return out; }
    const step = (b - a) / (n - 1);
    for (let i = 0; i < n; i++) out[i] = a + i * step;
    out[n - 1] = b;
    return out;
  }

  // transferMatrix(EnergyMaxMin, d, V) -> E, T(E)=|t|^2 (NaN where det check fails)
  function transferMatrix(Emin, Emax, d, V) {
    const Eres = 0.0005;
    const val = Math.floor((Emax - Emin) / Eres);
    const npoints = Math.max(val, 1000);
    const E = linspace(Emin, Emax, npoints);
    const T2 = new Float64Array(npoints);
    for (let i = 0; i < npoints; i++) {
      const T = totalT(E[i], V, d);
      const det = detLU(T);
      if (cabs([det[0] - 1, det[1]]) > 1e-9) { T2[i] = NaN; continue; }
      const r = cdiv([-T[2][0], -T[2][1]], T[3]);
      const t = cadd(T[0], cmul(T[1], r));
      T2[i] = t[0] * t[0] + t[1] * t[1];
    }
    return { E, T2 };
  }

  // scipy.signal.find_peaks with default arguments (_local_maxima_1d)
  function findPeaks(x) {
    const peaks = [];
    const n = x.length;
    let i = 1;
    const iMax = n - 1;
    while (i < iMax) {
      if (x[i - 1] < x[i]) {
        let ia = i + 1;
        while (ia < iMax && x[ia] === x[i]) ia++;
        if (x[ia] < x[i]) {
          const left = i, right = ia - 1;
          peaks.push(Math.floor((left + right) / 2));
          i = ia;
        }
      }
      i++;
    }
    return peaks;
  }

  // transferMatrix_z(E, d, V) -> z [m], amp (real), prob, VV, |t|^2
  function transferMatrixZ(E, d, V) {
    const Ttot = totalT(E, V, d);
    const r = cdiv([-Ttot[2][0], -Ttot[2][1]], Ttot[3]);
    const tt = cadd(Ttot[0], cmul(Ttot[1], r));

    const npoints = 1000;
    let L = 0; for (const x of d) L += x;
    const z = linspace(0, L, npoints);
    const amp = new Float64Array(npoints), prob = new Float64Array(npoints), VV = new Float64Array(npoints);
    let T = I2();
    let layer = 0;
    const layerZ = []; let acc = 0; for (const x of d) { acc += x; layerZ.push(acc); }
    const dz = z[1] - z[0];
    const v0 = [1, 0], v1 = r;
    prob[0] = 1 + (r[0] * r[0] + r[1] * r[1]);
    amp[0] = 1 + r[0];
    VV[0] = V[0];
    for (let k = 0; k < npoints - 2; k++) {
      if (z[k + 1] > layerZ[layer]) {
        T = mm(TM_P(E, V[layer], layerZ[layer] - z[k]), T);
        T = mm(TM_T(E, V[layer], V[layer + 1]), T);
        T = mm(TM_P(E, V[layer + 1], z[k + 1] - layerZ[layer]), T);
        layer += 1;
      } else if (z[k + 1] < layerZ[layer]) {
        T = mm(TM_P(E, V[layer], dz), T);
      } else if (z[k + 1] === layerZ[layer]) {
        T = mm(TM_P(E, V[layer], dz), T);
        T = mm(TM_T(E, V[layer], V[layer + 1]), T);
        layer += 1;
      } else break;
      VV[k + 1] = V[layer];
      const a1 = cadd(cmul(T[0], v0), cmul(T[1], v1));
      const a2 = cadd(cmul(T[2], v0), cmul(T[3], v1));
      const s = cadd(a1, a2);
      amp[k + 1] = s[0];
      prob[k + 1] = s[0] * s[0] + s[1] * s[1];
    }
    // numpy keeps amp[0] = 1 + r complex, so max|amp| uses |1 + r| there
    let maxAbs = Math.hypot(1 + r[0], r[1]);
    for (let i = 1; i < npoints; i++) maxAbs = Math.max(maxAbs, Math.abs(amp[i]));
    return { z, amp, prob, VV, maxAbs, T2: tt[0] * tt[0] + tt[1] * tt[1] };
  }

  // Full pipeline of plotIt (minus the drawing)
  function run(p) {
    const { V, d } = generateStruc(p.N, p.V0, p.wellWidth, p.delta || 0, p.barrierWidth);
    const Emin = p.Emin + 1e-10, Emax = p.Emax + 1e-10;
    const spec = transferMatrix(Emin, Emax, d, V);
    const peaks = findPeaks(spec.T2);
    return { V, d, spec, peaks, peakE: peaks.map((i) => spec.E[i]) };
  }

  return { generateStruc, transferMatrix, transferMatrixZ, findPeaks, run, TM_T, TM_P };
})();

if (typeof module !== 'undefined') module.exports = Core;
