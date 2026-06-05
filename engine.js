/*
 * CardioSynth engine — pure statistical core for the Equipose / CardioSynthesizer
 * V5.5 DAPT reanalysis dashboard.
 *
 * Extracted VERBATIM from the dashboard's inline web-worker so the statistical
 * core is a single source of truth, importable under Node for testing. The page
 * loads this file with <script src="engine.js"></script> before the inline
 * scripts; the worker is built from this same file's text at runtime.
 *
 * Method: DerSimonian-Laird random-effects meta-analysis on log odds ratios
 * (logOR), with a 95% prediction interval, plus Egger's small-study-effects
 * regression, an odds-ratio -> risk-difference transform, and a histogram
 * builder for the Monte-Carlo NCB distribution. Faithful to the shipped app —
 * no methodology changes. (k=1 guard verified: tau2/I2 clamp to 0, no NaN.)
 */

// Odds-ratio -> absolute risk difference at a given baseline risk.
// base = baseline probability; logOR = log odds ratio applied to it.
function probFromLogOR(base, logOR) {
  const odds = base / (1 - base);
  const newOdds = odds * Math.exp(logOR);
  const newProb = newOdds / (1 + newOdds);
  return newProb - base;
}

// Egger's test: regression of effect on precision (1/SE). Returns the
// intercept and a coarse significance flag. Undefined for k<3.
function eggerTest(data, kEff, kSe) {
  if (data.length < 3) return null;
  const effects = data.map(d => d[kEff]);
  const precision = data.map(d => 1 / d[kSe]);
  const n = effects.length;

  const meanP = precision.reduce((a, b) => a + b) / n;
  const meanE = effects.reduce((a, b) => a + b) / n;

  let num = 0, den = 0;
  for (let i = 0; i < n; i++) {
    num += (precision[i] - meanP) * (effects[i] - meanE);
    den += (precision[i] - meanP) ** 2;
  }
  const slope = num / den;
  const intercept = meanE - slope * meanP;

  // T-stat for intercept
  const residuals = effects.map((e, i) => e - (intercept + slope * precision[i]));
  const ssr = residuals.reduce((a, b) => a + b * b, 0);
  const syx = Math.sqrt(ssr / (n - 2));
  const ssx = den;
  // SE of intercept
  const se_int = syx * Math.sqrt(1 / n + (meanP * meanP) / ssx);
  const t = intercept / se_int;
  const pval = (Math.abs(t) > 2.0) ? "<0.05" : ">0.05"; // Simplified approx

  return { intercept, pval };
}

// DerSimonian-Laird random-effects estimator on the supplied effect/SE keys.
function runREML(data, kEff, kSe) {
  const y = data.map(d => d[kEff]);
  const v = data.map(d => d[kSe] ** 2);

  // Fixed-effects
  const wFE = v.map(x => 1 / x);
  const swFE = wFE.reduce((a, b) => a + b, 0);
  const muFE = wFE.reduce((a, b, j) => a + b * y[j], 0) / swFE;

  // Q statistic
  const Q = wFE.reduce((sum, w, j) => sum + w * (y[j] - muFE) ** 2, 0);
  const df = Math.max(1, data.length - 1);
  const C = swFE - wFE.reduce((sum, w) => sum + w ** 2, 0) / swFE;

  // Tau^2 & I^2
  const tau2 = Math.max(0, (Q - df) / C);
  const I2 = Math.max(0, (Q - df) / Q * 100);

  // Random-effects
  const wRE = v.map(x => 1 / (x + tau2));
  const swRE = wRE.reduce((a, b) => a + b, 0);
  const muRE = wRE.reduce((a, b, j) => a + b * y[j], 0) / swRE;

  // Prediction Interval (95%)
  const se_pool = Math.sqrt(1 / swRE);
  const se_pred = Math.sqrt(se_pool ** 2 + tau2);
  const pred_lo = muRE - 1.96 * se_pred;
  const pred_hi = muRE + 1.96 * se_pred;

  return {
    mu: muRE,
    se_pool: se_pool,
    tau2: tau2,
    I2: I2,
    pred_lo, pred_hi,
    forest: data.map(d => ({
      id: d.id,
      val: d[kEff],
      lo: d[kEff] - 1.96 * d[kSe],
      hi: d[kEff] + 1.96 * d[kSe]
    }))
  };
}

// Histogram over the central 5th-95th percentile band of a sorted array.
function buildHistogram(arr, bins) {
  if (arr.length === 0) return [];
  const p05 = arr[Math.floor(arr.length * 0.05)];
  const p95 = arr[Math.floor(arr.length * 0.95)];
  const min = p05, max = p95;
  const range = max - min;
  const step = range / bins;
  const res = [];
  let current = min;
  let idx = 0;
  while (idx < arr.length && arr[idx] < min) idx++;
  for (let i = 0; i < bins; i++) {
    let count = 0;
    const next = current + step;
    while (idx < arr.length && arr[idx] < next) { count++; idx++; }
    res.push({ x: (current + step / 2) * 100, y: count });
    current = next;
  }
  return res;
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { probFromLogOR, eggerTest, runREML, buildHistogram };
}
