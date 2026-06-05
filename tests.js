/*
 * Node tests for the CardioSynth engine. Run: node tests.js
 * Every expected value is hand-computed independently of the engine code.
 */
const { runREML, eggerTest, probFromLogOR, buildHistogram } = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  ok  ' + name); }
  else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

// --- HAND-WORKED 2-STUDY DerSimonian-Laird POOLING ---------------------------
// t1: y1=-0.16, s1=0.06 -> w1=1/0.0036 = 277.77778
// t2: y2=-0.21, s2=0.05 -> w2=1/0.0025 = 400
//   SumW = 677.77778
//   muFE = (277.77778*-0.16 + 400*-0.21)/677.77778
//        = (-44.44444 - 84)/677.77778 = -128.44444/677.77778 = -0.1895082
//   Q = 277.77778*(-0.16 - muFE)^2 + 400*(-0.21 - muFE)^2
//     = 277.77778*(0.0295082)^2 + 400*(-0.0204918)^2
//     = 277.77778*0.00087073 + 400*0.00041991 = 0.241870 + 0.167964 = 0.409834
//   df = k-1 = 1.  Q < df  ->  tau2 = max(0,(Q-1)/C) = 0,  I2 = 0.
//   tau2=0  => muRE = muFE = -0.1895082
//   se_pool = sqrt(1/SumW) = sqrt(1/677.77778) = sqrt(0.00147541) = 0.0384111
//   se_pred = sqrt(se_pool^2 + 0) = 0.0384111
//   pred_lo = -0.1895082 - 1.96*0.0384111 = -0.2647938
//   pred_hi = -0.1895082 + 1.96*0.0384111 = -0.1142225
const two = runREML(
  [{ id: "A", e: -0.16, s: 0.06 }, { id: "B", e: -0.21, s: 0.05 }], "e", "s");
ok('2study: muRE ~ -0.1895082', close(two.mu, -0.1895082, 1e-6), 'got ' + two.mu);
ok('2study: se_pool ~ 0.0384111', close(two.se_pool, 0.0384111, 1e-6), 'got ' + two.se_pool);
ok('2study: tau2 == 0 (Q<df)', two.tau2 === 0, 'got ' + two.tau2);
ok('2study: I2 == 0', two.I2 === 0, 'got ' + two.I2);
ok('2study: pred_lo ~ -0.2647938', close(two.pred_lo, -0.2647938, 1e-6), 'got ' + two.pred_lo);
ok('2study: pred_hi ~ -0.1142225', close(two.pred_hi, -0.1142225, 1e-6), 'got ' + two.pred_hi);
ok('2study: forest first lo = val-1.96*se', close(two.forest[0].lo, -0.16 - 1.96 * 0.06, 1e-9));
ok('2study: forest second hi = val+1.96*se', close(two.forest[1].hi, -0.21 + 1.96 * 0.05, 1e-9));

// --- HAND-WORKED HETEROGENEOUS 2-STUDY (Q > df, tau2 > 0) --------------------
// t1: y1=-0.50, s1=0.10 -> w1=100 ; t2: y2=0.10, s2=0.10 -> w2=100
//   SumW=200 ; muFE=(100*-0.5 + 100*0.1)/200 = -40/200 = -0.20
//   Q = 100*(-0.5+0.2)^2 + 100*(0.1+0.2)^2 = 100*0.09 + 100*0.09 = 18
//   df=1 ; C = SumW - (w1^2+w2^2)/SumW = 200 - (10000+10000)/200 = 200 - 100 = 100
//   tau2 = (18-1)/100 = 0.17
//   I2 = (18-1)/18 *100 = 94.4444%
//   wRE = 1/(0.01+0.17)=1/0.18=5.55556 each ; muRE = (5.55556*-0.5 + 5.55556*0.1)/11.11111 = -0.20
//   se_pool = sqrt(1/11.11111) = sqrt(0.09) = 0.30
const het = runREML(
  [{ id: "A", e: -0.50, s: 0.10 }, { id: "B", e: 0.10, s: 0.10 }], "e", "s");
ok('het: tau2 ~ 0.17', close(het.tau2, 0.17, 1e-9), 'got ' + het.tau2);
ok('het: I2 ~ 94.4444', close(het.I2, 94.44444, 1e-4), 'got ' + het.I2);
ok('het: muRE ~ -0.20', close(het.mu, -0.20, 1e-9), 'got ' + het.mu);
ok('het: se_pool ~ 0.30', close(het.se_pool, 0.30, 1e-9), 'got ' + het.se_pool);

// --- k=1 PASSTHROUGH (must NOT produce NaN) ---------------------------------
// single study: muRE=y, tau2=0, I2=0, se_pool=s
const one = runREML([{ id: "X", e: -0.20, s: 0.10 }], "e", "s");
ok('k1: mu == -0.20', close(one.mu, -0.20, 1e-12), 'got ' + one.mu);
ok('k1: tau2 == 0 (no NaN)', one.tau2 === 0 && !Number.isNaN(one.tau2), 'got ' + one.tau2);
ok('k1: I2 == 0 (no NaN)', one.I2 === 0 && !Number.isNaN(one.I2), 'got ' + one.I2);
ok('k1: se_pool == 0.10', close(one.se_pool, 0.10, 1e-12), 'got ' + one.se_pool);
ok('k1: nothing NaN', !Number.isNaN(one.mu) && !Number.isNaN(one.pred_lo) && !Number.isNaN(one.pred_hi));

// --- TWO IDENTICAL STUDIES => tau2=0, I2=0 -----------------------------------
// y=-0.20,s=0.10 each -> w=100 ; SumW=200 ; muFE=-0.20 ; Q=0 -> tau2=0, I2=0
//   se_pool = sqrt(1/200) = 0.0707107
const id2 = runREML([{ id: "A", e: -0.20, s: 0.10 }, { id: "B", e: -0.20, s: 0.10 }], "e", "s");
ok('identical: mu == -0.20', close(id2.mu, -0.20, 1e-12), 'got ' + id2.mu);
ok('identical: tau2 == 0', id2.tau2 === 0, 'got ' + id2.tau2);
ok('identical: I2 == 0', id2.I2 === 0, 'got ' + id2.I2);
ok('identical: se_pool == sqrt(1/200)', close(id2.se_pool, Math.sqrt(1 / 200), 1e-12), 'got ' + id2.se_pool);

// --- probFromLogOR ----------------------------------------------------------
// base=0.04, logOR=0 -> no change -> risk difference 0 exactly
ok('probFromLogOR(0.04,0) == 0', probFromLogOR(0.04, 0) === 0, 'got ' + probFromLogOR(0.04, 0));
// base=0.04, OR=0.5: odds0=0.04/0.96=0.0416667 ; newOdds=0.0208333 ;
//   newProb=0.0208333/1.0208333=0.0204082 ; diff=0.0204082-0.04=-0.0195918
ok('probFromLogOR(0.04, ln0.5) ~ -0.0195918',
   close(probFromLogOR(0.04, Math.log(0.5)), -0.0195918, 1e-6),
   'got ' + probFromLogOR(0.04, Math.log(0.5)));
// monotone: a protective logOR (<0) lowers risk -> negative diff
ok('probFromLogOR protective < 0', probFromLogOR(0.10, Math.log(0.7)) < 0);
ok('probFromLogOR harmful > 0', probFromLogOR(0.10, Math.log(1.5)) > 0);

// --- eggerTest --------------------------------------------------------------
// k<3 -> null (undefined)
ok('egger: k=2 -> null', eggerTest([{ e: 0.1, s: 0.1 }, { e: 0.2, s: 0.2 }], 'e', 's') === null);
// Perfectly symmetric funnel: effect independent of precision -> intercept ~ 0.
//   3 studies all y=-0.20, varying SE. Regression of y on 1/SE has slope 0,
//   intercept = mean(y) = -0.20 (since precision varies but y constant).
const eg = eggerTest(
  [{ e: -0.20, s: 0.05 }, { e: -0.20, s: 0.10 }, { e: -0.20, s: 0.20 }], 'e', 's');
ok('egger: constant effect -> intercept ~ -0.20', close(eg.intercept, -0.20, 1e-9), 'got ' + eg.intercept);

// --- buildHistogram ---------------------------------------------------------
// empty guard
ok('histogram: empty -> []', Array.isArray(buildHistogram([], 10)) && buildHistogram([], 10).length === 0);
// bins count: a uniformly spread sorted array yields exactly `bins` buckets
const sorted = Array.from({ length: 1000 }, (_, i) => i / 1000); // 0..0.999
const hist = buildHistogram(sorted, 20);
ok('histogram: returns requested bin count', hist.length === 20, 'got ' + hist.length);
ok('histogram: total counted <= input length', hist.reduce((a, b) => a + b.y, 0) <= sorted.length);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
