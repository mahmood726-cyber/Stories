/*
 * Minimal smoke test for the meta-analysis engine.
 * Run with:  node Evidence_Reversal_40/smoke.test.js
 * Exit code 0 = all checks passed, non-zero = a check failed.
 *
 * These checks pin a few hand-verifiable numerical anchors so an accidental
 * regression in the pooling / effect-size code is caught before publishing.
 */
'use strict';

const M = require('./meta-analysis.js');

let failures = 0;
function approx(label, got, want, tol) {
  tol = tol === undefined ? 1e-6 : tol;
  if (typeof got !== 'number' || !isFinite(got) || Math.abs(got - want) > tol) {
    console.error('FAIL ' + label + ': got ' + got + ', want ' + want + ' (tol ' + tol + ')');
    failures++;
  } else {
    console.log('ok   ' + label + ' = ' + got);
  }
}
function ok(label, cond) {
  if (cond) { console.log('ok   ' + label); } else { console.error('FAIL ' + label); failures++; }
}

// 1. Engine loads and exposes the public API.
ok('exports MetaAnalysis constructor', typeof M === 'function');
ok('exports EffectSizes', M.EffectSizes && typeof M.EffectSizes.logOR === 'function');

// 2. logOR is exact for a 2x2 table with no zero cells.
//    log((10*80)/(90*20)) = log(0.444444...) = -0.81093022
var es = M.EffectSizes.logOR(10, 90, 20, 80);
approx('logOR yi (10/90/20/80)', es.yi, -0.8109302162163288, 1e-9);
approx('logOR vi (10/90/20/80)', es.vi, 1 / 10 + 1 / 90 + 1 / 20 + 1 / 80, 1e-12);

// 3. Zero-cell continuity correction is applied ONLY when a cell is 0.
//    log((0.5*5.5)/(10.5*5.5)) = log(0.5/10.5) = -3.04452244
var zero = M.EffectSizes.logOR(0, 10, 5, 5);
approx('logOR zero-cell correction', zero.yi, Math.log(0.5 / 10.5), 1e-9);

// 4. Fixed-effects pooling matches the closed-form inverse-variance weighted mean.
var fe = new M({ method: 'DL', silent: true });
fe.addEffectSize({ yi: 0.2, vi: 0.04, study: 'A' });
fe.addEffectSize({ yi: 0.5, vi: 0.05, study: 'B' });
fe.addEffectSize({ yi: -0.1, vi: 0.03, study: 'C' });
var feRes = fe.runFixedEffectsModel();
var w = [1 / 0.04, 1 / 0.05, 1 / 0.03];
var sumW = w[0] + w[1] + w[2];
var wantFE = (w[0] * 0.2 + w[1] * 0.5 + w[2] * -0.1) / sumW;
approx('fixed-effects estimate', feRes.estimate, wantFE, 1e-9);
approx('fixed-effects se', feRes.se, Math.sqrt(1 / sumW), 1e-9);

// 5. Random-effects model runs and returns a finite estimate + non-negative tau^2.
var re = new M({ method: 'DL', silent: true });
re.addEffectSize({ yi: 0.2, vi: 0.04, study: 'A' });
re.addEffectSize({ yi: 0.5, vi: 0.05, study: 'B' });
re.addEffectSize({ yi: -0.1, vi: 0.03, study: 'C' });
var reRes = re.runRandomEffectsModel({ predictionInterval: true });
ok('random-effects estimate finite', isFinite(reRes.estimate));
ok('random-effects tau2 >= 0', reRes.tau2 >= 0);
ok('prediction interval lb < ub', reRes.pi_lb < reRes.pi_ub);

// 6. k=1 is guarded (pooling a single study must throw, not silently "pool").
var threw = false;
try {
  var single = new M({ silent: true });
  single.addEffectSize({ yi: 0.1, vi: 0.02, study: 'X' });
  single.runRandomEffectsModel();
} catch (e) { threw = true; }
ok('k=1 random-effects throws', threw);

// 7. addEffectSize must fail closed on a non-positive / non-finite variance.
//    vi=0 previously produced wi=1/0=Infinity -> pooled estimate=NaN, se=0
//    (silent corruption). vi<0 produced se=NaN. Both must now throw.
function throwsAdd(study) {
  try { new M({ silent: true }).addEffectSize(study); return false; }
  catch (e) { return true; }
}
ok('addEffectSize throws on vi = 0', throwsAdd({ yi: 0.2, vi: 0 }));
ok('addEffectSize throws on vi < 0', throwsAdd({ yi: 0.2, vi: -0.1 }));
ok('addEffectSize throws on vi = NaN', throwsAdd({ yi: 0.2, vi: NaN }));
ok('addEffectSize throws on vi = Infinity', throwsAdd({ yi: 0.2, vi: Infinity }));
ok('addEffectSize throws on yi = NaN', throwsAdd({ yi: NaN, vi: 0.04 }));
// A valid positive variance must still be accepted (guard is not over-broad).
ok('addEffectSize accepts vi > 0', !throwsAdd({ yi: 0.2, vi: 0.04 }));

// 8. Effect-size anchors (hand-derived exact values, no zero cells).
//    logRR(10,90,20,80): yi = log(10/100) - log(20/100) = log(0.5);
//    vi = (1/10 - 1/100) + (1/20 - 1/100) = 0.09 + 0.04 = 0.13.
var rr = M.EffectSizes.logRR(10, 90, 20, 80);
approx('logRR yi (10/90/20/80)', rr.yi, Math.log(0.5), 1e-12);
approx('logRR vi (10/90/20/80)', rr.vi, 0.13, 1e-12);
//    riskDiff(10,90,20,80): yi = 0.1 - 0.2 = -0.1;
//    vi = 0.1*0.9/100 + 0.2*0.8/100 = 0.0009 + 0.0016 = 0.0025.
var rd = M.EffectSizes.riskDiff(10, 90, 20, 80);
approx('riskDiff yi (10/90/20/80)', rd.yi, -0.1, 1e-12);
approx('riskDiff vi (10/90/20/80)', rd.vi, 0.0025, 1e-12);
//    fisherZ(r=0.5, n=28): yi = 0.5*log((1.5)/(0.5)) = 0.5*log(3); vi = 1/(28-3) = 0.04.
var fz = M.EffectSizes.fisherZ(0.5, 28);
approx('fisherZ yi (r=0.5,n=28)', fz.yi, 0.5 * Math.log(3), 1e-12);
approx('fisherZ vi (r=0.5,n=28)', fz.vi, 0.04, 1e-12);

// 9. SMD (Hedges' g): m1=10,sd1=2,n1=20, m2=8,sd2=2,n2=20.
//    pooledSD = sqrt((19*4 + 19*4)/38) = 2 exactly, so raw Cohen's d = 1 exactly.
//    Hedges' J (df=38) must lie in (0,1) and match the 1 - 3/(4*df-1) approximation
//    to ~3 decimals; g = J*d. This pins the exact J-correction snapshot.
var smd = M.EffectSizes.smd(10, 2, 20, 8, 2, 20);
approx('smd raw Cohen d = 1', smd.d, 1, 1e-12);
ok('smd J in (0,1)', smd.J > 0 && smd.J < 1);
approx('smd J ~ approximation', smd.J, 1 - 3 / (4 * 38 - 1), 1e-3);
approx('smd g = J*d', smd.yi, smd.J * smd.d, 1e-12);

// 10. DerSimonian-Laird tau^2 on the 3-study set is closed-form hand-verifiable.
//     Q = 4.595742, df = 2, C = sumW - sum(wi^2)/sumW = 51.0638 -> tau2 = (Q-2)/C.
var dl = new M({ method: 'DL', silent: true });
dl.addEffectSize({ yi: 0.2, vi: 0.04 });
dl.addEffectSize({ yi: 0.5, vi: 0.05 });
dl.addEffectSize({ yi: -0.1, vi: 0.03 });
var dlRes = dl.runRandomEffectsModel();
approx('DL tau2 (3-study)', dlRes.tau2, 0.05083333333333333, 1e-6);

// 11. Egger's test wrappers return NaN (not null) for k<3, and produce a stable
//     t/p snapshot on a fixed 5-study asymmetric dataset.
var eggerSmall = new M({ silent: true });
eggerSmall.addEffectSize({ yi: 0.3, vi: 0.02 });
eggerSmall.addEffectSize({ yi: 0.4, vi: 0.05 });
var esRes = eggerSmall.eggerTest();
ok('eggerTest returns NaN pval for k<3', typeof esRes.pval === 'number' && isNaN(esRes.pval));
var beggSmall = new M({ silent: true });
beggSmall.addEffectSize({ yi: 0.3, vi: 0.02 });
beggSmall.addEffectSize({ yi: 0.4, vi: 0.05 });
var bsRes = beggSmall.beggTest();
ok('beggTest returns NaN pval for k<3', typeof bsRes.pval === 'number' && isNaN(bsRes.pval));
var egger = new M({ silent: true });
egger.addEffectSize({ yi: 0.30, vi: 0.02 });
egger.addEffectSize({ yi: 0.45, vi: 0.05 });
egger.addEffectSize({ yi: 0.55, vi: 0.08 });
egger.addEffectSize({ yi: 0.70, vi: 0.12 });
egger.addEffectSize({ yi: 0.25, vi: 0.01 });
var egRes = egger.eggerTest();
approx('Egger t (fixed dataset)', egRes.tval, 17.772412685555636, 1e-6);
approx('Egger p (fixed dataset)', egRes.pval, 0.00038842192258958796, 1e-9);

if (failures > 0) {
  console.error('\n' + failures + ' check(s) failed');
  process.exit(1);
}
console.log('\nAll smoke checks passed');
