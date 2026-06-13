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

if (failures > 0) {
  console.error('\n' + failures + ' check(s) failed');
  process.exit(1);
}
console.log('\nAll smoke checks passed');
