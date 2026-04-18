/**
 * MetaJS - A JavaScript Library for Meta-Analysis
 *
 * Extracted from TruthCert-PairwisePro
 * Validated against R metafor package v4.8.0
 *
 * @version 1.3.0
 * @license MIT
 * @author TruthCert-PairwisePro Team
 *
 * Features:
 * - Effect size calculations (OR, RR, RD, SMD, MD, correlations, proportions)
 * - 8 tau² heterogeneity estimators (DL, REML, ML, PM, HS, SJ, HE, EB)
 * - Publication bias tests (Egger, Begg, Trim-and-Fill)
 * - Sensitivity analysis (Leave-one-out, influence diagnostics)
 * - HKSJ adjustment for small studies
 * - Prediction intervals
 * - Confidence intervals for tau² and I² (Q-profile method)
 * - Fixed-effects model
 * - Meta-regression with multiple moderators
 * - Subgroup analysis with Q-between test
 * - Convenience summary methods
 * - Export to JSON/CSV formats
 *
 * Usage:
 *   const meta = new MetaAnalysis();
 *   meta.calculateEffectSizes(studies, 'OR');
 *   const results = meta.runRandomEffectsModel({ method: 'REML' });
 */

(function(global) {
  'use strict';

  // ============================================
  // CONSTANTS
  // ============================================

  const TOLERANCE = Object.freeze({
    CONVERGENCE: 1e-10,
    ZERO_CHECK: 1e-10,
    SMALL_VALUE: 1e-30,
    BETA_EPS: 3e-7,
    PIVOT: 1e-12  // Minimum pivot value for matrix operations
  });

  const MAX_ITERATIONS = Object.freeze({
    REML: 200,
    ML: 100,
    PM: 100,
    EB: 100,
    BETA_CF: 100,
    BISECTION: 100
  });

  // Critical values for confidence intervals (frozen for safety)
  const CRITICAL_VALUES = Object.freeze({
    Z_95: 1.959963984540054,  // normalQuantile(0.975), exact value for reproducibility
    Z_99: 2.575829303548901,  // normalQuantile(0.995)
    Z_90: 1.644853626951472   // normalQuantile(0.95)
  });

  // ============================================
  // UTILITY FUNCTIONS
  // ============================================

  const sum = arr => arr.reduce((a, b) => a + b, 0);
  const mean = arr => sum(arr) / arr.length;
  const variance = arr => {
    const m = mean(arr);
    return sum(arr.map(x => Math.pow(x - m, 2))) / (arr.length - 1);
  };

  /**
   * Validate input data structure
   * @param {Array} data - Array of study objects
   * @param {String} context - Context for error messages
   */
  const validateData = (data, context = 'data') => {
    if (!Array.isArray(data)) {
      throw new Error(`${context} must be an array`);
    }
    if (data.length < 2) {
      throw new Error(`${context} must contain at least 2 studies`);
    }
  };

  /**
   * Validate numeric value
   * @param {Number} value - Value to validate
   * @param {String} name - Name for error messages
   */
  const validateNumeric = (value, name) => {
    if (typeof value !== 'number' || !isFinite(value)) {
      throw new Error(`${name} must be a finite number, got: ${value}`);
    }
  };

  /**
   * Validate non-negative value
   * @param {Number} value - Value to validate
   * @param {String} name - Name for error messages
   */
  const validateNonNegative = (value, name) => {
    validateNumeric(value, name);
    if (value < 0) {
      throw new Error(`${name} must be non-negative, got: ${value}`);
    }
  };

  /**
   * Validate positive value (> 0)
   * @param {Number} value - Value to validate
   * @param {String} name - Name for error messages
   */
  const validatePositive = (value, name) => {
    validateNumeric(value, name);
    if (value <= 0) {
      throw new Error(`${name} must be positive, got: ${value}`);
    }
  };

  /**
   * Validate correlation coefficient (-1 to 1)
   * @param {Number} value - Value to validate
   * @param {String} name - Name for error messages
   */
  const validateCorrelation = (value, name) => {
    validateNumeric(value, name);
    if (value < -1 || value > 1) {
      throw new Error(`${name} must be between -1 and 1, got: ${value}`);
    }
  };

  /**
   * Validate 2x2 table cell (non-negative integer-like)
   * @param {Number} value - Value to validate
   * @param {String} name - Name for error messages
   */
  const validateCell = (value, name) => {
    validateNonNegative(value, name);
  };

  // ============================================
  // STATISTICAL DISTRIBUTION FUNCTIONS
  // ============================================

  /**
   * Error function approximation (Abramowitz & Stegun)
   */
  const erf = x => {
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    const t = 1 / (1 + p * x);
    const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return sign * y;
  };

  /**
   * Standard normal CDF
   */
  const normalCDF = (x, mu = 0, sigma = 1) => {
    const z = (x - mu) / sigma;
    return 0.5 * (1 + erf(z / Math.sqrt(2)));
  };

  /**
   * Normal quantile function (inverse CDF)
   * Uses Abramowitz & Stegun approximation
   */
  const normalQuantile = (p, mu = 0, sigma = 1) => {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (p === 0.5) return mu;

    const a = [
      -3.969683028665376e1, 2.209460984245205e2,
      -2.759285104469687e2, 1.383577518672690e2,
      -3.066479806614716e1, 2.506628277459239e0
    ];
    const b = [
      -5.447609879822406e1, 1.615858368580409e2,
      -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1
    ];
    const c = [
      -7.784894002430293e-3, -3.223964580411365e-1,
      -2.400758277161838e0, -2.549732539343734e0,
      4.374664141464968e0, 2.938163982698783e0
    ];
    const d = [
      7.784695709041462e-3, 3.224671290700398e-1,
      2.445134137142996e0, 3.754408661907416e0
    ];

    const pLow = 0.02425;
    const pHigh = 1 - pLow;
    let q, r;

    if (p < pLow) {
      q = Math.sqrt(-2 * Math.log(p));
      return mu + sigma * (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p <= pHigh) {
      q = p - 0.5;
      r = q * q;
      return mu + sigma * (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    } else {
      q = Math.sqrt(-2 * Math.log(1 - p));
      return mu + sigma * -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
  };

  /**
   * Log gamma function (Lanczos approximation)
   */
  const lgamma = x => {
    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.001208650973866179, -0.000005395239384953];
    let y = x, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += c[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  };

  /**
   * Beta continued fraction for incomplete beta
   */
  const betacf = (x, a, b) => {
    const maxIter = MAX_ITERATIONS.BETA_CF;
    const eps = TOLERANCE.BETA_EPS;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < TOLERANCE.SMALL_VALUE) d = TOLERANCE.SMALL_VALUE;
    d = 1 / d;
    let h = d;

    for (let m = 1; m <= maxIter; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d;
      if (Math.abs(d) < TOLERANCE.SMALL_VALUE) d = TOLERANCE.SMALL_VALUE;
      c = 1 + aa / c;
      if (Math.abs(c) < TOLERANCE.SMALL_VALUE) c = TOLERANCE.SMALL_VALUE;
      d = 1 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d;
      if (Math.abs(d) < TOLERANCE.SMALL_VALUE) d = TOLERANCE.SMALL_VALUE;
      c = 1 + aa / c;
      if (Math.abs(c) < TOLERANCE.SMALL_VALUE) c = TOLERANCE.SMALL_VALUE;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < eps) break;
    }
    return h;
  };

  /**
   * Regularized incomplete beta function I_x(a,b)
   */
  const incompleteBeta = (x, a, b) => {
    if (x === 0) return 0;
    if (x === 1) return 1;

    const bt = Math.exp(
      lgamma(a + b) - lgamma(a) - lgamma(b) +
      a * Math.log(x) + b * Math.log(1 - x)
    );

    if (x < (a + 1) / (a + b + 2)) {
      return bt * betacf(x, a, b) / a;
    }
    return 1 - bt * betacf(1 - x, b, a) / b;
  };

  /**
   * Regularized lower incomplete gamma function P(a, x) = γ(a,x) / Γ(a)
   * This is the chi-square CDF when a = k/2 and x = chi2/2
   */
  const regularizedGammaP = (a, x) => {
    if (x < 0 || a <= 0) return 0;
    if (x === 0) return 0;

    if (x < a + 1) {
      // Series representation
      let sum = 1 / a;
      let term = 1 / a;
      for (let n = 1; n < 200; n++) {
        term *= x / (a + n);
        sum += term;
        if (Math.abs(term) < Math.abs(sum) * 1e-14) break;
      }
      return Math.exp(-x + a * Math.log(x) - lgamma(a)) * sum;
    } else {
      // Use relation P(a,x) = 1 - Q(a,x) where Q is upper regularized gamma
      let b = x + 1 - a;
      let c = 1 / TOLERANCE.SMALL_VALUE;
      let d = 1 / b;
      let h = d;

      for (let i = 1; i < 200; i++) {
        const an = -i * (i - a);
        b += 2;
        d = an * d + b;
        if (Math.abs(d) < TOLERANCE.SMALL_VALUE) d = TOLERANCE.SMALL_VALUE;
        c = b + an / c;
        if (Math.abs(c) < TOLERANCE.SMALL_VALUE) c = TOLERANCE.SMALL_VALUE;
        d = 1 / d;
        const del = d * c;
        h *= del;
        if (Math.abs(del - 1) < 1e-14) break;
      }

      return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
    }
  };

  /**
   * Chi-square CDF using regularized incomplete gamma
   * P(X ≤ x) = P(k/2, x/2) where P is regularized lower incomplete gamma
   */
  const chiSquareCDF = (x, df) => {
    if (x <= 0) return 0;
    if (df <= 0) return 0;
    return regularizedGammaP(df / 2, x / 2);
  };

  /**
   * Student's t-distribution CDF
   */
  const tCDF = (t, df) => {
    if (df <= 0) return NaN;
    const x = df / (df + t * t);
    const beta = 0.5 * incompleteBeta(x, df / 2, 0.5);
    return t > 0 ? 1 - beta : beta;
  };

  /**
   * Student's t-distribution quantile
   * Uses lookup table for small df, Cornish-Fisher expansion otherwise
   */
  const tQuantile = (p, df) => {
    if (df === Infinity) return normalQuantile(p);
    if (p <= 0 || p >= 1) return NaN;

    // Lookup table for exact values at common probability levels
    const tTable = {
      1: { 0.975: 12.7062, 0.95: 6.3138, 0.90: 3.0777, 0.025: -12.7062, 0.05: -6.3138, 0.10: -3.0777 },
      2: { 0.975: 4.3027, 0.95: 2.9200, 0.90: 1.8856, 0.025: -4.3027, 0.05: -2.9200, 0.10: -1.8856 },
      3: { 0.975: 3.1824, 0.95: 2.3534, 0.90: 1.6377, 0.025: -3.1824, 0.05: -2.3534, 0.10: -1.6377 },
      4: { 0.975: 2.7764, 0.95: 2.1318, 0.90: 1.5332, 0.025: -2.7764, 0.05: -2.1318, 0.10: -1.5332 },
      5: { 0.975: 2.5706, 0.95: 2.0150, 0.90: 1.4759, 0.025: -2.5706, 0.05: -2.0150, 0.10: -1.4759 }
    };

    const roundedDf = Math.round(df);
    if (roundedDf >= 1 && roundedDf <= 5 && tTable[roundedDf]) {
      // Check for exact match
      const pRound = Math.round(p * 1000) / 1000;
      if (tTable[roundedDf][pRound] !== undefined) {
        return tTable[roundedDf][pRound];
      }
      // Also check 1-p for symmetry
      const pComp = Math.round((1 - p) * 1000) / 1000;
      if (tTable[roundedDf][pComp] !== undefined) {
        return -tTable[roundedDf][pComp];
      }
    }

    // Cornish-Fisher expansion for larger df
    const x = normalQuantile(p);
    const g1 = (x * x * x + x) / 4;
    const g2 = (5 * Math.pow(x, 5) + 16 * x * x * x + 3 * x) / 96;
    const g3 = (3 * Math.pow(x, 7) + 19 * Math.pow(x, 5) + 17 * x * x * x - 15 * x) / 384;
    const g4 = (79 * Math.pow(x, 9) + 776 * Math.pow(x, 7) + 1482 * Math.pow(x, 5) - 1920 * x * x * x - 945 * x) / 92160;

    return x + g1 / df + g2 / (df * df) + g3 / (df * df * df) + g4 / (df * df * df * df);
  };

  // ============================================
  // EFFECT SIZE CALCULATIONS
  // ============================================

  const EffectSizes = {
    /**
     * Log Odds Ratio with continuity correction
     * @param {Number} a - Events in treatment group
     * @param {Number} b - Non-events in treatment group
     * @param {Number} c - Events in control group
     * @param {Number} d - Non-events in control group
     * @param {Number} cc - Continuity correction (default 0.5)
     */
    logOR: function(a, b, c, d, cc = 0.5) {
      // Input validation
      validateCell(a, 'a (treatment events)');
      validateCell(b, 'b (treatment non-events)');
      validateCell(c, 'c (control events)');
      validateCell(d, 'd (control non-events)');
      validateNonNegative(cc, 'continuity correction');

      // Apply continuity correction if any cell is 0
      if (a === 0 || b === 0 || c === 0 || d === 0) {
        a += cc; b += cc; c += cc; d += cc;
      }

      const yi = Math.log((a * d) / (b * c));
      const vi = 1 / a + 1 / b + 1 / c + 1 / d;

      return { yi, vi, se: Math.sqrt(vi), measure: 'OR' };
    },

    /**
     * Log Risk Ratio
     */
    logRR: function(a, b, c, d, cc = 0.5) {
      // Input validation
      validateCell(a, 'a (treatment events)');
      validateCell(b, 'b (treatment non-events)');
      validateCell(c, 'c (control events)');
      validateCell(d, 'd (control non-events)');
      validateNonNegative(cc, 'continuity correction');

      const n1 = a + b, n2 = c + d;
      if (n1 <= 0 || n2 <= 0) {
        throw new Error('Group sizes must be positive');
      }

      let a_adj = a, c_adj = c;

      if (a === 0 || c === 0) {
        a_adj = a + cc;
        c_adj = c + cc;
      }

      const yi = Math.log(a_adj / n1) - Math.log(c_adj / n2);
      const vi = (1 / a_adj - 1 / n1) + (1 / c_adj - 1 / n2);

      return { yi, vi, se: Math.sqrt(vi), measure: 'RR' };
    },

    /**
     * Risk Difference
     */
    riskDiff: function(a, b, c, d) {
      // Input validation
      validateCell(a, 'a (treatment events)');
      validateCell(b, 'b (treatment non-events)');
      validateCell(c, 'c (control events)');
      validateCell(d, 'd (control non-events)');

      const n1 = a + b, n2 = c + d;
      if (n1 <= 0 || n2 <= 0) {
        throw new Error('Group sizes must be positive');
      }

      const p1 = a / n1, p2 = c / n2;
      const yi = p1 - p2;
      const vi = (p1 * (1 - p1)) / n1 + (p2 * (1 - p2)) / n2;

      return { yi, vi, se: Math.sqrt(vi), measure: 'RD' };
    },

    /**
     * Standardized Mean Difference (Hedges' g)
     * Uses exact J correction factor
     */
    smd: function(m1, sd1, n1, m2, sd2, n2) {
      // Input validation
      validateNumeric(m1, 'mean1');
      validateNumeric(m2, 'mean2');
      validatePositive(sd1, 'sd1');
      validatePositive(sd2, 'sd2');
      validatePositive(n1, 'n1');
      validatePositive(n2, 'n2');

      if (n1 < 2 || n2 < 2) {
        throw new Error('Sample sizes must be at least 2 for SMD calculation');
      }

      const pooledSD = Math.sqrt(((n1 - 1) * sd1 * sd1 + (n2 - 1) * sd2 * sd2) / (n1 + n2 - 2));
      const d = (m1 - m2) / pooledSD;

      // Hedges' correction factor J (exact formula using gamma function)
      const df = n1 + n2 - 2;
      // Approximation: J ≈ 1 - 3/(4*df - 1) is accurate for df > 10
      // For small df, use: J = gamma(df/2) / (sqrt(df/2) * gamma((df-1)/2))
      let J;
      if (df > 50) {
        J = 1 - 3 / (4 * df - 1);
      } else {
        J = Math.exp(lgamma(df / 2) - Math.log(Math.sqrt(df / 2)) - lgamma((df - 1) / 2));
      }

      const g = J * d;
      const vi = (n1 + n2) / (n1 * n2) + (g * g) / (2 * (n1 + n2));

      return { yi: g, vi, se: Math.sqrt(vi), measure: 'SMD', d, J };
    },

    /**
     * Raw Mean Difference
     */
    meanDiff: function(m1, sd1, n1, m2, sd2, n2) {
      // Input validation
      validateNumeric(m1, 'mean1');
      validateNumeric(m2, 'mean2');
      validatePositive(sd1, 'sd1');
      validatePositive(sd2, 'sd2');
      validatePositive(n1, 'n1');
      validatePositive(n2, 'n2');

      const yi = m1 - m2;
      const vi = (sd1 * sd1) / n1 + (sd2 * sd2) / n2;

      return { yi, vi, se: Math.sqrt(vi), measure: 'MD' };
    },

    /**
     * Fisher's z transformation for correlations
     */
    fisherZ: function(r, n) {
      // Input validation
      validateCorrelation(r, 'correlation r');
      validatePositive(n, 'sample size n');
      if (n < 4) {
        throw new Error('Sample size must be at least 4 for Fisher z transformation');
      }
      if (Math.abs(r) === 1) {
        throw new Error('Correlation cannot be exactly -1 or 1 for Fisher z transformation');
      }

      const yi = 0.5 * Math.log((1 + r) / (1 - r));
      const vi = 1 / (n - 3);

      return { yi, vi, se: Math.sqrt(vi), measure: 'ZCOR', r };
    },

    /**
     * Raw correlation (untransformed)
     */
    correlation: function(r, n) {
      // Input validation
      validateCorrelation(r, 'correlation r');
      validatePositive(n, 'sample size n');
      if (n < 2) {
        throw new Error('Sample size must be at least 2');
      }

      const vi = Math.pow(1 - r * r, 2) / (n - 1);
      return { yi: r, vi, se: Math.sqrt(vi), measure: 'COR', r };
    },

    /**
     * Back-transform Fisher's z to correlation
     */
    zToR: function(z) {
      return (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
    },

    /**
     * Logit proportion
     */
    logitProp: function(events, n) {
      // Input validation
      validateNonNegative(events, 'events');
      validatePositive(n, 'sample size n');
      if (events > n) {
        throw new Error('Events cannot exceed sample size');
      }

      let p = events / n;

      // Continuity correction for 0 or 1
      if (p === 0) p = 0.5 / n;
      if (p === 1) p = (n - 0.5) / n;

      const yi = Math.log(p / (1 - p));
      const vi = 1 / (n * p * (1 - p));

      return { yi, vi, se: Math.sqrt(vi), measure: 'PLO', p };
    },

    /**
     * Raw proportion
     */
    proportion: function(events, n) {
      // Input validation
      validateNonNegative(events, 'events');
      validatePositive(n, 'sample size n');
      if (events > n) {
        throw new Error('Events cannot exceed sample size');
      }

      const p = events / n;
      const vi = (p * (1 - p)) / n;
      return { yi: p, vi, se: Math.sqrt(vi), measure: 'PR', p };
    },

    /**
     * Freeman-Tukey double arcsine transformation
     */
    freemanTukey: function(events, n) {
      // Input validation
      validateNonNegative(events, 'events');
      validatePositive(n, 'sample size n');
      if (events > n) {
        throw new Error('Events cannot exceed sample size');
      }

      const yi = Math.asin(Math.sqrt(events / (n + 1))) +
        Math.asin(Math.sqrt((events + 1) / (n + 1)));
      const vi = 1 / (n + 0.5);

      return { yi, vi, se: Math.sqrt(vi), measure: 'PFT', p: events / n };
    },

    /**
     * Back-transform Freeman-Tukey to proportion
     */
    ftToProp: function(z, n) {
      return 0.5 * (1 - Math.sign(Math.cos(z)) *
        Math.sqrt(1 - Math.pow(Math.sin(z) + (Math.sin(z) - 1 / Math.sin(z)) / n, 2)));
    },

    /**
     * Arcsine square root transformation
     */
    arcsineProp: function(events, n) {
      // Input validation
      validateNonNegative(events, 'events');
      validatePositive(n, 'sample size n');
      if (events > n) {
        throw new Error('Events cannot exceed sample size');
      }

      const p = events / n;
      const yi = Math.asin(Math.sqrt(p));
      const vi = 1 / (4 * n);

      return { yi, vi, se: Math.sqrt(vi), measure: 'PAS', p };
    },

    /**
     * Hazard Ratio from CI
     */
    hazardRatio: function(hr, ci_lower, ci_upper) {
      // Input validation
      validatePositive(hr, 'hazard ratio');
      validatePositive(ci_lower, 'CI lower bound');
      validatePositive(ci_upper, 'CI upper bound');
      if (ci_lower >= ci_upper) {
        throw new Error('CI lower bound must be less than CI upper bound');
      }
      if (hr < ci_lower || hr > ci_upper) {
        throw new Error('Hazard ratio must be within the confidence interval');
      }

      const yi = Math.log(hr);
      const se = (Math.log(ci_upper) - Math.log(ci_lower)) / 3.92;
      const vi = se * se;

      return { yi, vi, se, measure: 'HR', hr };
    }
  };

  // ============================================
  // TAU-SQUARED ESTIMATORS
  // ============================================

  const Tau2Estimators = {
    /**
     * DerSimonian-Laird moment-based estimator
     */
    DL: function(yi, vi) {
      const k = yi.length;
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const sumW2 = sum(wi.map(w => w * w));

      const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;
      const Q = sum(yi.map((y, i) => wi[i] * Math.pow(y - theta, 2)));

      const C = sumW - sumW2 / sumW;
      const tau2 = Math.max(0, (Q - (k - 1)) / C);

      return { tau2, Q, df: k - 1, method: 'DL', converged: true };
    },

    /**
     * Restricted Maximum Likelihood (REML)
     * Uses Fisher scoring with proper REML likelihood
     */
    REML: function(yi, vi, maxIter = MAX_ITERATIONS.REML, tol = TOLERANCE.CONVERGENCE) {
      const k = yi.length;
      let tau2 = this.DL(yi, vi).tau2;

      // If DL gives 0, check if it should be 0
      const Q_dl = this.DL(yi, vi).Q;
      if (Q_dl <= k - 1) {
        return { tau2: 0, converged: true, iterations: 0, method: 'REML' };
      }

      // Fisher scoring for REML
      for (let iter = 0; iter < maxIter; iter++) {
        const wi = vi.map(v => 1 / (v + tau2));
        const sumW = sum(wi);
        const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

        const wi2 = wi.map(w => w * w);
        const sumW2 = sum(wi2);

        // REML score: derivative of REML log-likelihood
        const residuals = yi.map(y => y - theta);
        const Q_tau = sum(residuals.map((r, i) => wi2[i] * r * r));
        const traceP = sumW - sumW2 / sumW;

        const score = -0.5 * traceP + 0.5 * Q_tau;

        // Fisher information
        const traceP2 = sum(wi.map(w => {
          const p = w - w * w / sumW;
          return p * p;
        }));
        const info = 0.5 * traceP2;

        if (info < TOLERANCE.ZERO_CHECK) break;

        // Newton-Raphson update with damping
        let delta = score / info;
        let tau2New = tau2 + delta;

        // Ensure non-negative
        if (tau2New < 0) {
          tau2New = tau2 / 2;
        }

        // Check convergence
        if (Math.abs(tau2New - tau2) < tol || Math.abs(delta) < tol) {
          return { tau2: Math.max(0, tau2New), converged: true, iterations: iter + 1, method: 'REML' };
        }

        tau2 = Math.max(0, tau2New);
      }

      // Bisection fallback
      let lo = 0, hi = Math.max(1, tau2 * 5, variance(yi));

      for (let iter = 0; iter < MAX_ITERATIONS.BISECTION; iter++) {
        const mid = (lo + hi) / 2;
        const wi = vi.map(v => 1 / (v + mid));
        const sumW = sum(wi);
        const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

        const wi2 = wi.map(w => w * w);
        const sumW2 = sum(wi2);
        const residuals = yi.map(y => y - theta);
        const Q_tau = sum(residuals.map((r, i) => wi2[i] * r * r));
        const traceP = sumW - sumW2 / sumW;

        const score = -0.5 * traceP + 0.5 * Q_tau;

        if (Math.abs(score) < tol || (hi - lo) < tol) {
          return { tau2: mid, converged: true, iterations: iter, method: 'REML' };
        }

        if (score > 0) lo = mid;
        else hi = mid;
      }

      return { tau2: (lo + hi) / 2, converged: false, method: 'REML' };
    },

    /**
     * Maximum Likelihood
     */
    ML: function(yi, vi, maxIter = MAX_ITERATIONS.ML, tol = TOLERANCE.CONVERGENCE) {
      const k = yi.length;
      let tau2 = this.DL(yi, vi).tau2;
      if (tau2 < TOLERANCE.ZERO_CHECK) tau2 = 0.01;

      for (let iter = 0; iter < maxIter; iter++) {
        const wi = vi.map(v => 1 / (v + tau2));
        const sumW = sum(wi);
        const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

        const wi2 = wi.map(w => w * w);
        const residuals = yi.map(y => Math.pow(y - theta, 2));

        const score = -0.5 * sum(wi) + 0.5 * sum(wi.map((w, i) => w * w * residuals[i]));
        const info = 0.5 * sum(wi2);

        const tau2New = Math.max(0, tau2 + score / info);

        if (Math.abs(tau2New - tau2) < tol) {
          return { tau2: tau2New, converged: true, iterations: iter + 1, method: 'ML' };
        }

        tau2 = tau2New;
      }

      return { tau2, converged: false, iterations: maxIter, method: 'ML' };
    },

    /**
     * Paule-Mandel iterative estimator
     */
    PM: function(yi, vi, maxIter = MAX_ITERATIONS.PM, tol = TOLERANCE.CONVERGENCE) {
      const k = yi.length;
      let tau2 = this.DL(yi, vi).tau2;

      const Q_of_tau2 = (t2) => {
        const wi = vi.map(v => 1 / (v + t2));
        const sumW = sum(wi);
        const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;
        return sum(yi.map((y, i) => wi[i] * Math.pow(y - theta, 2)));
      };

      if (Q_of_tau2(0) <= k - 1) {
        return { tau2: 0, converged: true, iterations: 0, method: 'PM' };
      }

      // Iterative approach
      for (let iter = 0; iter < 20; iter++) {
        const wi = vi.map(v => 1 / (v + tau2));
        const sumW = sum(wi);
        const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;
        const Q = sum(yi.map((y, i) => wi[i] * Math.pow(y - theta, 2)));

        if (Math.abs(Q - (k - 1)) < tol) {
          return { tau2, converged: true, iterations: iter + 1, method: 'PM' };
        }

        const sumW2 = sum(wi.map(w => w * w));
        const C = sumW - sumW2 / sumW;
        tau2 = Math.max(0, tau2 + (Q - (k - 1)) / C);
      }

      // Bisection fallback
      let lower = 0, upper = Math.max(tau2 * 10, 100);

      for (let iter = 0; iter < MAX_ITERATIONS.BISECTION; iter++) {
        const mid = (lower + upper) / 2;
        const Q_mid = Q_of_tau2(mid);

        if (Math.abs(Q_mid - (k - 1)) < tol || (upper - lower) < tol) {
          return { tau2: mid, converged: true, iterations: iter + 21, method: 'PM' };
        }

        if (Q_mid > k - 1) lower = mid;
        else upper = mid;
      }

      return { tau2: (lower + upper) / 2, converged: false, method: 'PM' };
    },

    /**
     * Hunter-Schmidt estimator
     */
    HS: function(yi, vi) {
      const k = yi.length;
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;
      const Q = sum(yi.map((y, i) => wi[i] * Math.pow(y - theta, 2)));
      const tau2 = Math.max(0, (Q - (k - 1)) / sumW);

      return { tau2, Q, method: 'HS', converged: true };
    },

    /**
     * Sidik-Jonkman two-step estimator
     */
    SJ: function(yi, vi) {
      const k = yi.length;
      const y_bar = mean(yi);
      const s2 = sum(yi.map(y => Math.pow(y - y_bar, 2))) / (k - 1);
      const v_bar = mean(vi);
      const tau2_init = Math.max(0, s2 - v_bar);

      const wi = vi.map(v => 1 / (v + tau2_init));
      const sumW = sum(wi);
      const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;
      const Q = sum(yi.map((y, i) => wi[i] * Math.pow(y - theta, 2)));

      const sumW2 = sum(wi.map(w => w * w));
      const C = sumW - sumW2 / sumW;
      const tau2 = Math.max(0, (Q - (k - 1)) / C);

      return { tau2, Q, method: 'SJ', converged: true };
    },

    /**
     * Hedges estimator (unweighted)
     */
    HE: function(yi, vi) {
      const k = yi.length;
      const y_bar = mean(yi);
      const Q_unweighted = sum(yi.map(y => Math.pow(y - y_bar, 2)));
      const tau2 = Math.max(0, Q_unweighted / (k - 1) - mean(vi));

      return { tau2, Q: Q_unweighted, method: 'HE', converged: true };
    },

    /**
     * Empirical Bayes estimator (Morris)
     */
    EB: function(yi, vi, maxIter = MAX_ITERATIONS.EB, tol = TOLERANCE.CONVERGENCE) {
      const k = yi.length;
      let tau2 = this.DL(yi, vi).tau2;

      for (let iter = 0; iter < maxIter; iter++) {
        const wi = vi.map(v => 1 / (v + tau2));
        const sumW = sum(wi);
        const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

        const num = sum(yi.map((y, i) => Math.pow(y - theta, 2) / Math.pow(vi[i] + tau2, 2)));
        const den = sum(vi.map(v => 1 / Math.pow(v + tau2, 2)));

        const tau2New = Math.max(0, num / den - 1 / sumW);

        if (Math.abs(tau2New - tau2) < tol) {
          return { tau2: tau2New, converged: true, iterations: iter + 1, method: 'EB' };
        }

        tau2 = tau2New;
      }

      return { tau2, converged: false, iterations: maxIter, method: 'EB' };
    }
  };

  // ============================================
  // HETEROGENEITY STATISTICS
  // ============================================

  const Heterogeneity = {
    /**
     * Calculate I², H², and Q statistics
     * Uses tau²-based I² formula to match metafor: I² = tau²/(tau² + s²)
     */
    calculate: function(yi, vi, tau2) {
      const k = yi.length;
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const sumW2 = sum(wi.map(w => w * w));
      const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;
      const Q = sum(yi.map((y, i) => wi[i] * Math.pow(y - theta, 2)));
      const df = k - 1;
      const p_Q = 1 - chiSquareCDF(Q, df);

      // Typical within-study variance: s² = (k-1) / C
      const C = sumW - sumW2 / sumW;
      const typicalV = df / C;

      // I² using tau²-based formula to match metafor
      // I² = tau² / (tau² + s²) where s² is typical sampling variance
      const I2 = tau2 < TOLERANCE.ZERO_CHECK ? 0 : tau2 / (tau2 + typicalV);

      // H² = Q / df
      const H2 = Q / df;
      const H = Math.sqrt(H2);

      return {
        Q, df, p_Q,
        I2, H2, H,
        tau2, tau: Math.sqrt(tau2),
        typicalV
      };
    },

    /**
     * Q-profile confidence interval for tau² (Viechtbauer, 2007)
     * Uses chi-square distribution to find bounds
     */
    tau2ConfidenceInterval: function(yi, vi, Q, level = 0.95) {
      const k = yi.length;
      const df = k - 1;
      const alpha = 1 - level;

      // Fixed-effects weights
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const sumW2 = sum(wi.map(w => w * w));
      const C = sumW - sumW2 / sumW;

      // Chi-square critical values
      const chi2_lb = chiSquareQuantile(1 - alpha / 2, df);
      const chi2_ub = chiSquareQuantile(alpha / 2, df);

      // tau² bounds from Q-profile method
      // tau² = (Q - chi²) / C
      const tau2_lb = Math.max(0, (Q - chi2_lb) / C);
      const tau2_ub = Math.max(0, (Q - chi2_ub) / C);

      return {
        tau2_lb: tau2_lb,
        tau2_ub: tau2_ub,
        level: level,
        method: 'Q-profile'
      };
    },

    /**
     * Profile Likelihood CI for tau²
     * More accurate for REML/ML estimates than Q-profile method
     *
     * @param {Array} yi - Effect sizes
     * @param {Array} vi - Variances
     * @param {Number} tau2_hat - Point estimate of tau²
     * @param {String} method - Estimation method ('REML' or 'ML')
     * @param {Number} level - Confidence level (default: 0.95)
     * @returns {Object} tau² CI via profile likelihood
     */
    profileLikelihoodCI: function(yi, vi, tau2_hat, method = 'REML', level = 0.95) {
      const k = yi.length;
      const alpha = 1 - level;
      const critVal = chiSquareQuantile(level, 1);  // 3.841 for 95%

      // Calculate log-likelihood at tau2_hat
      const logLik = function(tau2) {
        tau2 = Math.max(0, tau2);
        const wi = vi.map(v => 1 / (v + tau2));
        const sumW = sum(wi);
        const mu = sum(yi.map((y, i) => wi[i] * y)) / sumW;
        const Q_w = sum(yi.map((y, i) => wi[i] * Math.pow(y - mu, 2)));

        // Log-likelihood
        let ll = -0.5 * sum(vi.map((v, i) => Math.log(v + tau2)));
        ll -= 0.5 * Q_w;

        if (method === 'REML') {
          // REML adjustment
          ll -= 0.5 * Math.log(sumW);
        }

        return ll;
      };

      const ll_max = logLik(tau2_hat);
      const target = ll_max - critVal / 2;  // 2 * (ll_max - ll) = critVal

      // Find lower bound (search from 0 to tau2_hat)
      let tau2_lb = 0;
      const searchLower = function() {
        let lo = 0;
        let hi = tau2_hat;
        for (let iter = 0; iter < 100; iter++) {
          const mid = (lo + hi) / 2;
          if (logLik(mid) < target) {
            lo = mid;
          } else {
            hi = mid;
          }
          if (hi - lo < 1e-8) break;
        }
        return lo;
      };
      tau2_lb = logLik(0) >= target ? 0 : searchLower();

      // Find upper bound (search from tau2_hat upward)
      let tau2_ub = tau2_hat;
      const searchUpper = function() {
        // First find an upper bound where ll < target
        let hi = tau2_hat * 2 + 1;
        while (logLik(hi) >= target && hi < 1000) {
          hi *= 2;
        }
        if (hi >= 1000) return hi;  // Unbounded

        let lo = tau2_hat;
        for (let iter = 0; iter < 100; iter++) {
          const mid = (lo + hi) / 2;
          if (logLik(mid) >= target) {
            lo = mid;
          } else {
            hi = mid;
          }
          if (hi - lo < 1e-8) break;
        }
        return hi;
      };
      tau2_ub = searchUpper();

      return {
        tau2_lb: tau2_lb,
        tau2_ub: tau2_ub,
        level: level,
        method: 'Profile Likelihood (' + method + ')'
      };
    },

    /**
     * Confidence interval for I² derived from tau² CI
     * Uses the formula: I² = tau² / (tau² + typicalV)
     */
    I2ConfidenceInterval: function(yi, vi, tau2CI, level = 0.95) {
      const k = yi.length;
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const sumW2 = sum(wi.map(w => w * w));
      const df = k - 1;
      const C = sumW - sumW2 / sumW;
      const typicalV = df / C;

      // Transform tau² bounds to I² bounds
      const I2_lb = tau2CI.tau2_lb / (tau2CI.tau2_lb + typicalV);
      const I2_ub = tau2CI.tau2_ub / (tau2CI.tau2_ub + typicalV);

      return {
        I2_lb: I2_lb,
        I2_ub: I2_ub,
        level: level,
        method: 'Derived from Q-profile tau² CI'
      };
    }
  };

  /**
   * Chi-square quantile function (inverse CDF)
   * Uses Newton-Raphson iteration
   */
  const chiSquareQuantile = (p, df) => {
    if (p <= 0) return 0;
    if (p >= 1) return Infinity;

    // Initial guess using Wilson-Hilferty approximation
    const z = normalQuantile(p);
    let x = df * Math.pow(1 - 2 / (9 * df) + z * Math.sqrt(2 / (9 * df)), 3);
    x = Math.max(0.001, x);

    // Newton-Raphson iteration
    for (let i = 0; i < 50; i++) {
      const cdf = chiSquareCDF(x, df);
      const pdf = chiSquarePDF(x, df);
      if (pdf < 1e-20) break;

      const dx = (cdf - p) / pdf;
      x = Math.max(0.001, x - dx);

      if (Math.abs(dx) < 1e-10) break;
    }

    return x;
  };

  /**
   * Chi-square PDF
   */
  const chiSquarePDF = (x, df) => {
    if (x <= 0 || df <= 0) return 0;
    const k = df / 2;
    return Math.exp((k - 1) * Math.log(x) - x / 2 - k * Math.log(2) - lgamma(k));
  };

  // ============================================
  // POOLED ESTIMATE CALCULATION
  // ============================================

  const PooledEstimate = {
    /**
     * Calculate pooled effect with weights
     */
    calculate: function(yi, vi, tau2, useHKSJ = false) {
      const k = yi.length;

      // Random-effects weights
      const wi = vi.map(v => 1 / (v + tau2));
      const sumW = sum(wi);

      // Pooled estimate
      const estimate = sum(yi.map((y, i) => wi[i] * y)) / sumW;

      // Standard error
      let se = Math.sqrt(1 / sumW);

      // HKSJ adjustment
      let hksj_mult = 1;
      if (useHKSJ && k > 1) {
        const Q_re = sum(yi.map((y, i) => wi[i] * Math.pow(y - estimate, 2)));
        hksj_mult = Q_re / (k - 1);
        // Use raw multiplier to match metafor (not max(1, hksj_mult))
        se = se * Math.sqrt(hksj_mult);
      }

      // Confidence interval
      const df = useHKSJ ? k - 1 : Infinity;
      const tCrit = useHKSJ ? tQuantile(0.975, df) : CRITICAL_VALUES.Z_95;
      const ci_lb = estimate - tCrit * se;
      const ci_ub = estimate + tCrit * se;

      // P-value
      const testStat = estimate / se;
      const pval = useHKSJ
        ? 2 * (1 - tCDF(Math.abs(testStat), df))
        : 2 * (1 - normalCDF(Math.abs(testStat)));

      return {
        estimate, se, ci_lb, ci_ub, zval: testStat, pval,
        weights: wi,
        hksj: useHKSJ,
        hksj_mult,
        df
      };
    },

    /**
     * Prediction interval
     */
    predictionInterval: function(estimate, se, tau2, k, level = 0.95) {
      const df = Math.max(1, k - 2);
      const alpha = 1 - level;
      const tCrit = tQuantile(1 - alpha / 2, df);
      const pi_se = Math.sqrt(se * se + tau2);

      return {
        pi_lb: estimate - tCrit * pi_se,
        pi_ub: estimate + tCrit * pi_se,
        pi_se, df
      };
    }
  };

  // ============================================
  // PUBLICATION BIAS TESTS
  // ============================================

  const PublicationBias = {
    /**
     * Egger's regression test for funnel plot asymmetry
     * Matches metafor's regtest(x, predictor="sei")
     * Fits: yi = b0 + b1*sei using WLS with weights = 1/vi
     * Tests whether b1 (slope) differs from 0
     */
    egger: function(yi, vi) {
      const k = yi.length;
      const sei = vi.map(v => Math.sqrt(v));
      const wi = vi.map(v => 1 / v);  // WLS weights

      // Weighted least squares: yi ~ 1 + sei with weights wi
      const sumW = sum(wi);
      const sumWx = sum(sei.map((s, i) => wi[i] * s));
      const sumWy = sum(yi.map((y, i) => wi[i] * y));
      const sumWxy = sum(yi.map((y, i) => wi[i] * y * sei[i]));
      const sumWx2 = sum(sei.map((s, i) => wi[i] * s * s));

      // Normal equations for WLS: (X'WX) * beta = X'Wy
      const det = sumW * sumWx2 - sumWx * sumWx;
      const b0 = (sumWx2 * sumWy - sumWx * sumWxy) / det;  // intercept
      const b1 = (sumW * sumWxy - sumWx * sumWy) / det;    // slope

      // Weighted residuals and MSE
      const residuals = yi.map((y, i) => y - b0 - b1 * sei[i]);
      const sse = sum(residuals.map((r, i) => wi[i] * r * r));
      const mse = sse / (k - 2);

      // Variance-covariance of coefficients: mse * (X'WX)^-1
      const var_b0 = mse * sumWx2 / det;
      const var_b1 = mse * sumW / det;

      // Test for slope (b1) different from 0 - this is the Egger test
      const t_val = b1 / Math.sqrt(var_b1);
      const p_val = 2 * (1 - tCDF(Math.abs(t_val), k - 2));

      // Note: metafor returns the intercept (b0) as "zval" in the output
      // but the actual test is for the slope (b1)
      return {
        intercept: b0,  // Estimated effect at sei=0
        slope: b1,      // Coefficient on sei (tested for asymmetry)
        se: Math.sqrt(var_b1),
        t: t_val,
        p: p_val,
        df: k - 2,
        method: 'Egger'
      };
    },

    /**
     * Begg's rank correlation test
     * Calculates Kendall's tau between standardized residuals and variances
     */
    begg: function(yi, vi) {
      const k = yi.length;
      const sei = vi.map(v => Math.sqrt(v));

      // Calculate pooled estimate (fixed effects)
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

      // Standardized residuals (yi - theta) / sei
      const ri = yi.map((y, i) => (y - theta) / sei[i]);

      // Kendall's tau between residuals and variances
      let concordant = 0, discordant = 0;
      for (let i = 0; i < k - 1; i++) {
        for (let j = i + 1; j < k; j++) {
          const sign1 = Math.sign(ri[i] - ri[j]);
          const sign2 = Math.sign(vi[i] - vi[j]);
          if (sign1 * sign2 > 0) concordant++;
          else if (sign1 * sign2 < 0) discordant++;
        }
      }

      const n_pairs = k * (k - 1) / 2;
      const tau = (concordant - discordant) / n_pairs;

      // Approximate p-value using normal approximation
      const se_tau = Math.sqrt(2 * (2 * k + 5) / (9 * k * (k - 1)));
      const z = tau / se_tau;
      const p = 2 * (1 - normalCDF(Math.abs(z)));

      return { tau, z, p, method: 'Begg' };
    },

    /**
     * Trim and Fill method
     */
    trimFill: function(yi, vi, tau2, side = 'auto', maxIter = 100) {
      const k = yi.length;

      // Calculate pooled estimate
      const wi = vi.map(v => 1 / (v + tau2));
      const sumW = sum(wi);
      let theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

      // Determine side if auto
      let estimatedSide = side;
      if (side === 'auto') {
        const eggerResult = this.egger(yi, vi);
        estimatedSide = eggerResult.slope > 0 ? 'left' : 'right';
      }

      // Calculate R0 (number of missing studies) using L0 estimator
      const centered = yi.map(y => y - theta);
      const ranks = this._rank(centered.map(Math.abs));

      let T_plus = 0;
      for (let i = 0; i < k; i++) {
        if ((estimatedSide === 'right' && centered[i] > 0) ||
          (estimatedSide === 'left' && centered[i] < 0)) {
          T_plus += ranks[i];
        }
      }

      // L0 estimator for k0
      const k0 = Math.max(0, Math.round((4 * T_plus - k * (k + 1)) / (2 * k - 1)));

      if (k0 === 0) {
        return {
          k0: 0,
          estimate: theta,
          se: Math.sqrt(1 / sumW),
          imputed: [],
          method: 'Trim-Fill'
        };
      }

      // Impute missing studies by reflection
      const sorted = yi.map((y, i) => ({ y, v: vi[i], centered: y - theta }))
        .sort((a, b) => Math.abs(b.centered) - Math.abs(a.centered));

      const imputed = [];
      for (let i = 0; i < k0 && i < sorted.length; i++) {
        const reflected = 2 * theta - sorted[i].y;
        imputed.push({ yi: reflected, vi: sorted[i].v });
      }

      // Recalculate with imputed studies
      const yi_all = [...yi, ...imputed.map(s => s.yi)];
      const vi_all = [...vi, ...imputed.map(s => s.vi)];
      const wi_all = vi_all.map(v => 1 / (v + tau2));
      const sumW_all = sum(wi_all);
      const theta_adj = sum(yi_all.map((y, i) => wi_all[i] * y)) / sumW_all;

      return {
        k0,
        estimate: theta_adj,
        se: Math.sqrt(1 / sumW_all),
        imputed,
        side: estimatedSide,
        method: 'Trim-Fill'
      };
    },

    /**
     * Fail-Safe N calculations (File Drawer Analysis)
     * Three methods: Rosenthal, Orwin, Rosenberg
     *
     * @param {Array} yi - Effect sizes
     * @param {Array} vi - Variances
     * @param {Object} options - Configuration options
     * @returns {Object} Fail-safe N values for each method
     */
    failSafeN: function(yi, vi, options) {
      options = options || {};
      const k = yi.length;
      const alpha = options.alpha || 0.05;
      const targetES = options.targetES || 0.2;  // Small effect threshold for Orwin

      // Calculate z-values for each study
      const zi = yi.map((y, i) => y / Math.sqrt(vi[i]));

      // Pooled z and effect
      const sumZ = sum(zi);
      const wi = vi.map(v => 1 / v);
      const sumW = sum(wi);
      const pooledES = sum(yi.map((y, i) => wi[i] * y)) / sumW;

      // Critical z-value
      const z_crit = normalQuantile(1 - alpha / 2);

      // ==== Rosenthal's Fail-Safe N ====
      // N_fs = (sum(z) / z_crit)^2 - k
      // Number of null studies needed to make p > alpha
      const rosenthal = Math.max(0, Math.round(
        Math.pow(sumZ / z_crit, 2) - k
      ));

      // ==== Orwin's Fail-Safe N ====
      // N_fs = k * ((ES_pooled - ES_target) / ES_target)
      // Number of null studies needed to reduce effect to target
      let orwin = 0;
      if (Math.abs(pooledES) > Math.abs(targetES) && targetES !== 0) {
        orwin = Math.max(0, Math.round(
          k * ((Math.abs(pooledES) - Math.abs(targetES)) / Math.abs(targetES))
        ));
      }

      // ==== Rosenberg's Fail-Safe N (Weighted) ====
      // Uses weighted sum of z-scores
      // N_fs = (sum(z)^2 - k * z_crit^2) / z_crit^2
      const rosenberg = Math.max(0, Math.round(
        (Math.pow(sumZ, 2) - k * Math.pow(z_crit, 2)) / Math.pow(z_crit, 2)
      ));

      // Interpretation thresholds (5k + 10 rule from Rosenthal)
      const rosenthalThreshold = 5 * k + 10;

      return {
        rosenthal: {
          N: rosenthal,
          interpretation: rosenthal > rosenthalThreshold ?
            'Robust: FSN exceeds 5k+10 threshold' :
            'Potentially vulnerable to publication bias'
        },
        orwin: {
          N: orwin,
          targetES: targetES,
          interpretation: orwin > k ?
            'Robust: Many null studies needed to reduce effect' :
            'Effect may be sensitive to publication bias'
        },
        rosenberg: {
          N: rosenberg,
          interpretation: rosenberg > rosenthalThreshold ?
            'Robust: Weighted FSN exceeds threshold' :
            'Potentially vulnerable to publication bias'
        },
        pooledEffect: pooledES,
        k: k,
        alpha: alpha,
        method: 'Fail-Safe N'
      };
    },

    _rank: function(arr) {
      const sorted = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
      const ranks = new Array(arr.length);
      for (let i = 0; i < sorted.length; i++) {
        ranks[sorted[i].i] = i + 1;
      }
      return ranks;
    }
  };

  // ============================================
  // SENSITIVITY ANALYSIS
  // ============================================

  const Sensitivity = {
    /**
     * Leave-one-out analysis
     */
    leaveOneOut: function(yi, vi, method = 'REML') {
      const k = yi.length;
      const results = [];

      for (let i = 0; i < k; i++) {
        const yi_loo = yi.filter((_, j) => j !== i);
        const vi_loo = vi.filter((_, j) => j !== i);

        const tau2_result = Tau2Estimators[method](yi_loo, vi_loo);
        const pooled = PooledEstimate.calculate(yi_loo, vi_loo, tau2_result.tau2);
        const het = Heterogeneity.calculate(yi_loo, vi_loo, tau2_result.tau2);

        results.push({
          omitted: i,
          estimate: pooled.estimate,
          se: pooled.se,
          ci_lb: pooled.ci_lb,
          ci_ub: pooled.ci_ub,
          tau2: tau2_result.tau2,
          I2: het.I2
        });
      }

      return results;
    },

    /**
     * Influence diagnostics
     */
    influence: function(yi, vi, tau2) {
      const k = yi.length;
      const wi = vi.map(v => 1 / (v + tau2));
      const sumW = sum(wi);
      const theta = sum(yi.map((y, i) => wi[i] * y)) / sumW;

      const results = [];

      for (let i = 0; i < k; i++) {
        // Standardized residual
        const resid = yi[i] - theta;
        const hi = wi[i] / sumW;  // Hat value
        const resid_se = Math.sqrt((vi[i] + tau2) * (1 - hi));
        const rstudent = resid / resid_se;

        // Cook's distance
        const cooksd = (wi[i] * resid * resid) / (2 * (1 - hi) * (1 - hi));

        // DFBETAS
        const dfbetas = (wi[i] * resid) / (sumW * Math.sqrt(1 - hi));

        results.push({
          study: i,
          residual: resid,
          rstudent,
          hat: hi,
          cooksd,
          dfbetas
        });
      }

      return results;
    }
  };

  // ============================================
  // MAIN META-ANALYSIS CLASS
  // ============================================

  /**
   * MetaAnalysis class for conducting meta-analyses
   * @class
   */
  function MetaAnalysis(options) {
    if (!(this instanceof MetaAnalysis)) {
      return new MetaAnalysis(options);
    }

    options = options || {};
    this.options = {
      method: options.method || 'REML',
      measure: options.measure || 'OR',
      useHKSJ: options.useHKSJ !== undefined ? options.useHKSJ : true,
      continuityCorrection: options.continuityCorrection || 0.5,
      level: options.level || 0.95,
      silent: options.silent || false,  // Suppress console warnings
      onWarning: options.onWarning || null,  // Callback for warnings
      onProgress: options.onProgress || null  // Progress callback for long operations
    };
    this.data = [];
    this.model = null;
    this._warnings = [];  // Store warnings for later retrieval
  }

  /**
   * Internal warning handler - respects silent mode and callbacks
   * @private
   */
  MetaAnalysis.prototype._warn = function(message, context) {
    const warning = { message: message, context: context || null, timestamp: new Date() };
    this._warnings.push(warning);

    if (this.options.onWarning) {
      this.options.onWarning(warning);
    } else if (!this.options.silent) {
      console.warn('MetaAnalysis: ' + message);
    }
  };

  /**
   * Get all accumulated warnings
   * @returns {Array} Array of warning objects
   */
  MetaAnalysis.prototype.getWarnings = function() {
    return this._warnings.slice();
  };

  /**
   * Clear accumulated warnings
   * @returns {MetaAnalysis} Returns this for method chaining
   */
  MetaAnalysis.prototype.clearWarnings = function() {
    this._warnings = [];
    return this;
  };

  /**
   * Set analysis option (chainable)
   * @param {String} key - Option name
   * @param {*} value - Option value
   * @returns {MetaAnalysis} Returns this for method chaining
   */
  MetaAnalysis.prototype.setOption = function(key, value) {
    if (this.options.hasOwnProperty(key) || ['method', 'measure', 'level', 'useHKSJ', 'continuityCorrection', 'silent'].indexOf(key) >= 0) {
      this.options[key] = value;
    } else {
      throw new Error('Unknown option: ' + key);
    }
    return this;
  };

  /**
   * Add a pre-calculated effect size to the dataset (chainable)
   * @param {Object} study - Study object with yi, vi, and study name
   * @returns {MetaAnalysis} Returns this for method chaining
   */
  MetaAnalysis.prototype.addEffectSize = function(study) {
    if (typeof study.yi !== 'number' || typeof study.vi !== 'number') {
      throw new Error('Study must have numeric yi (effect size) and vi (variance)');
    }
    if (!study.study) {
      study.study = 'Study ' + (this.data.length + 1);
    }
    study.se = Math.sqrt(study.vi);
    this.data.push(study);
    this.model = null;  // Invalidate cached model
    return this;
  };

  /**
   * Remove a study from the dataset (chainable)
   * @param {Number|String} identifier - Index or study name
   * @returns {MetaAnalysis} Returns this for method chaining
   */
  MetaAnalysis.prototype.removeStudy = function(identifier) {
    const idx = typeof identifier === 'number' ? identifier :
      this.data.findIndex(function(d) { return d.study === identifier; });

    if (idx < 0 || idx >= this.data.length) {
      throw new Error('Study not found: ' + identifier);
    }

    this.data.splice(idx, 1);
    this.model = null;  // Invalidate cached model
    return this;
  };

  /**
   * Reset the analysis (clear all data and model)
   * @returns {MetaAnalysis} Returns this for method chaining
   */
  MetaAnalysis.prototype.reset = function() {
    this.data = [];
    this.model = null;
    this._warnings = [];
    return this;
  };

  /**
   * Calculate effect sizes from raw study data
   * @param {Array} studies - Array of study objects
   * @param {String} measure - Effect size measure (OR, RR, RD, SMD, MD, etc.)
   * @returns {MetaAnalysis} Returns this for method chaining
   */
  MetaAnalysis.prototype.calculateEffectSizes = function(studies, measure) {
    // Validate input
    if (!Array.isArray(studies)) {
      throw new Error('Studies must be an array');
    }
    if (studies.length === 0) {
      throw new Error('Studies array cannot be empty');
    }

    measure = measure || 'OR';
    this.data = [];
    const cc = this.options.continuityCorrection;
    const self = this;

    for (let idx = 0; idx < studies.length; idx++) {
      const study = studies[idx];
      let result;
      const studyName = study.study || study.name || ('Study ' + (this.data.length + 1));

      try {
        if (study.yi !== undefined && study.vi !== undefined) {
          result = { yi: study.yi, vi: study.vi, se: Math.sqrt(study.vi) };
        } else if (study.events_t !== undefined && study.events_c !== undefined) {
          const a = study.events_t;
          const b = study.n_t - study.events_t;
          const c = study.events_c;
          const d = study.n_c - study.events_c;

          switch (measure) {
            case 'OR': result = EffectSizes.logOR(a, b, c, d, cc); break;
            case 'RR': result = EffectSizes.logRR(a, b, c, d, cc); break;
            case 'RD': result = EffectSizes.riskDiff(a, b, c, d); break;
            default: result = EffectSizes.logOR(a, b, c, d, cc);
          }
        } else if (study.mean_t !== undefined && study.mean_c !== undefined) {
          switch (measure) {
            case 'SMD':
              result = EffectSizes.smd(
                study.mean_t, study.sd_t, study.n_t,
                study.mean_c, study.sd_c, study.n_c
              );
              break;
            case 'MD':
            default:
              result = EffectSizes.meanDiff(
                study.mean_t, study.sd_t, study.n_t,
                study.mean_c, study.sd_c, study.n_c
              );
          }
        } else if (study.r !== undefined && study.n !== undefined) {
          switch (measure) {
            case 'COR': result = EffectSizes.correlation(study.r, study.n); break;
            case 'ZCOR':
            default: result = EffectSizes.fisherZ(study.r, study.n);
          }
        } else if (study.events !== undefined && study.n !== undefined) {
          switch (measure) {
            case 'PFT': result = EffectSizes.freemanTukey(study.events, study.n); break;
            case 'PAS': result = EffectSizes.arcsineProp(study.events, study.n); break;
            case 'PR': result = EffectSizes.proportion(study.events, study.n); break;
            case 'PLO':
            default: result = EffectSizes.logitProp(study.events, study.n);
          }
        }

        if (result && isFinite(result.yi) && isFinite(result.vi) && result.vi > 0) {
          // Copy all study properties (for moderators, risk of bias, etc.)
          const studyData = {
            study: studyName,
            yi: result.yi,
            vi: result.vi,
            se: Math.sqrt(result.vi),
            n_t: study.n_t,
            n_c: study.n_c
          };

          // Preserve additional properties (moderators, risk of bias)
          for (const key in study) {
            if (study.hasOwnProperty(key) && studyData[key] === undefined) {
              studyData[key] = study[key];
            }
          }

          this.data.push(studyData);
        }
      } catch (e) {
        self._warn('Error calculating effect size for study ' + studyName + ': ' + e.message, { study: studyName, error: e });
      }
    }

    this.options.measure = measure;
    return this;  // Method chaining
  };

  /**
   * Run random-effects meta-analysis model
   * @param {Object} options - Model options
   * @returns {Object} Model results
   */
  MetaAnalysis.prototype.runRandomEffectsModel = function(options) {
    options = options || {};
    const method = options.method || this.options.method || 'REML';
    const hksj = options.hksj !== undefined ? options.hksj : false;
    const predictionInterval = options.predictionInterval !== undefined ? options.predictionInterval : false;
    const profileLikelihoodCI = options.profileLikelihoodCI !== undefined ? options.profileLikelihoodCI : false;

    if (this.data.length < 2) {
      throw new Error('At least 2 studies required for meta-analysis');
    }

    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const k = yi.length;

    // Estimate tau²
    const tau2Result = Tau2Estimators[method](yi, vi);
    const tau2 = tau2Result.tau2;

    // Calculate random-effects weights
    const wi = vi.map(function(v) { return 1 / (v + tau2); });
    const sumW = sum(wi);

    // Pooled estimate
    const estimate = sum(yi.map(function(y, i) { return wi[i] * y; })) / sumW;
    let se = Math.sqrt(1 / sumW);

    // Fixed-effects Q statistic (for heterogeneity test)
    const wi_fe = vi.map(function(v) { return 1 / v; });
    const sumW_fe = sum(wi_fe);
    const sumW2_fe = sum(wi_fe.map(function(w) { return w * w; }));
    const theta_fe = sum(yi.map(function(y, i) { return wi_fe[i] * y; })) / sumW_fe;
    const Q = sum(yi.map(function(y, i) { return wi_fe[i] * Math.pow(y - theta_fe, 2); }));
    const Q_df = k - 1;

    // I² calculation using tau²-based formula to match metafor
    // I² = tau² / (tau² + s²) where s² is the "typical" sampling variance
    // s² = (k-1) / C where C = sum(wi) - sum(wi²)/sum(wi)
    let I2;
    if (tau2 < TOLERANCE.ZERO_CHECK) {
      I2 = 0;
    } else {
      const C = sumW_fe - sumW2_fe / sumW_fe;
      const typicalV = Q_df / C;  // typical sampling variance
      I2 = tau2 / (tau2 + typicalV);
    }

    // H² from Q
    const H2 = Q / Q_df;

    // Chi-square p-value for Q using proper incomplete gamma
    const Q_pval = 1 - chiSquareCDF(Q, Q_df);

    // HKSJ adjustment (Hartung-Knapp-Sidik-Jonkman)
    // metafor uses the raw multiplier without max(1, ...) constraint
    let hksj_ci_lb, hksj_ci_ub, hksj_tval, hksj_pval;
    if (hksj) {
      const Q_re = sum(yi.map(function(y, i) { return wi[i] * Math.pow(y - estimate, 2); }));
      const hksj_mult = Q_re / Q_df;
      // Use raw multiplier to match metafor (not max(1, hksj_mult))
      const hksj_se = se * Math.sqrt(hksj_mult);
      const tCrit = tQuantile(0.975, Q_df);
      hksj_ci_lb = estimate - tCrit * hksj_se;
      hksj_ci_ub = estimate + tCrit * hksj_se;
      hksj_tval = estimate / hksj_se;
      hksj_pval = 2 * (1 - tCDF(Math.abs(hksj_tval), Q_df));
    }

    // Standard CI (z-based)
    const ci_lb = estimate - CRITICAL_VALUES.Z_95 * se;
    const ci_ub = estimate + CRITICAL_VALUES.Z_95 * se;
    const zval = estimate / se;
    const pval = 2 * (1 - normalCDF(Math.abs(zval)));

    // Prediction interval
    let pi_lb, pi_ub;
    if (predictionInterval) {
      if (tau2 < TOLERANCE.ZERO_CHECK) {
        // No heterogeneity: PI = CI
        pi_lb = estimate - CRITICAL_VALUES.Z_95 * se;
        pi_ub = estimate + CRITICAL_VALUES.Z_95 * se;
      } else {
        const pi_se = Math.sqrt(se * se + tau2);
        const pi_df = Math.max(1, k - 2);
        const pi_tCrit = tQuantile(0.975, pi_df);
        pi_lb = estimate - pi_tCrit * pi_se;
        pi_ub = estimate + pi_tCrit * pi_se;
      }
    }

    // Store weights in data
    for (let i = 0; i < this.data.length; i++) {
      this.data[i].weight = wi[i] / sumW * 100;
    }

    // Calculate confidence intervals for tau² and I²
    // Use profile likelihood if requested (more accurate for REML/ML), otherwise Q-profile
    let tau2CI;
    if (profileLikelihoodCI && (method === 'REML' || method === 'ML')) {
      tau2CI = Heterogeneity.profileLikelihoodCI(yi, vi, tau2, method, 0.95);
    } else {
      tau2CI = Heterogeneity.tau2ConfidenceInterval(yi, vi, Q, 0.95);
    }
    const I2CI = Heterogeneity.I2ConfidenceInterval(yi, vi, tau2CI, 0.95);

    this.model = {
      estimate: estimate,
      se: se,
      ci_lb: ci_lb,
      ci_ub: ci_ub,
      zval: zval,
      pval: pval,
      tau2: tau2,
      tau: Math.sqrt(tau2),
      tau2_ci_lb: tau2CI.tau2_lb,
      tau2_ci_ub: tau2CI.tau2_ub,
      tau2_ci_method: tau2CI.method,
      I2: I2,
      I2_ci_lb: I2CI.I2_lb,
      I2_ci_ub: I2CI.I2_ub,
      H2: H2,
      Q: Q,
      Q_df: Q_df,
      Q_pval: Q_pval,
      k: k,
      method: method,
      hksj: hksj,
      converged: tau2Result.converged
    };

    if (hksj) {
      this.model.hksj_ci_lb = hksj_ci_lb;
      this.model.hksj_ci_ub = hksj_ci_ub;
      this.model.hksj_tval = hksj_tval;
      this.model.hksj_pval = hksj_pval;
    }

    if (predictionInterval) {
      this.model.pi_lb = pi_lb;
      this.model.pi_ub = pi_ub;
    }

    this.model.predictionInterval = predictionInterval;

    return this.model;
  };

  /**
   * Run fixed-effects meta-analysis model
   * @param {Object} options - Model options
   * @returns {Object} Model results
   */
  MetaAnalysis.prototype.runFixedEffectsModel = function(options) {
    options = options || {};

    if (this.data.length < 2) {
      throw new Error('At least 2 studies required for meta-analysis');
    }

    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const k = yi.length;

    // Fixed-effects weights (inverse variance)
    const wi = vi.map(function(v) { return 1 / v; });
    const sumW = sum(wi);
    const sumW2 = sum(wi.map(function(w) { return w * w; }));

    // Pooled estimate
    const estimate = sum(yi.map(function(y, i) { return wi[i] * y; })) / sumW;

    // Standard error
    const se = Math.sqrt(1 / sumW);

    // Q statistic for heterogeneity
    const Q = sum(yi.map(function(y, i) { return wi[i] * Math.pow(y - estimate, 2); }));
    const Q_df = k - 1;
    const Q_pval = 1 - chiSquareCDF(Q, Q_df);

    // I² and H²
    const C = sumW - sumW2 / sumW;
    const typicalV = Q_df / C;
    const tau2_DL = Math.max(0, (Q - Q_df) / C);
    const I2 = tau2_DL < TOLERANCE.ZERO_CHECK ? 0 : tau2_DL / (tau2_DL + typicalV);
    const H2 = Q / Q_df;

    // Confidence interval (z-based)
    const ci_lb = estimate - CRITICAL_VALUES.Z_95 * se;
    const ci_ub = estimate + CRITICAL_VALUES.Z_95 * se;
    const zval = estimate / se;
    const pval = 2 * (1 - normalCDF(Math.abs(zval)));

    // Store weights in data
    for (let i = 0; i < this.data.length; i++) {
      this.data[i].weight = wi[i] / sumW * 100;
    }

    this.model = {
      estimate: estimate,
      se: se,
      ci_lb: ci_lb,
      ci_ub: ci_ub,
      zval: zval,
      pval: pval,
      tau2: 0,  // Fixed-effects assumes no between-study variance
      tau: 0,
      I2: I2,   // Still report for information
      H2: H2,
      Q: Q,
      Q_df: Q_df,
      Q_pval: Q_pval,
      k: k,
      method: 'FE',
      modelType: 'fixed'
    };

    return this.model;
  };

  /**
   * Meta-regression with one or more moderators
   * @param {String|Array} moderators - Moderator variable name(s) in study data
   * @param {Object} options - Model options
   * @returns {Object} Meta-regression results
   */
  MetaAnalysis.prototype.metaRegression = function(moderators, options) {
    options = options || {};
    const method = options.method || this.options.method || 'REML';

    if (this.data.length < 3) {
      throw new Error('At least 3 studies required for meta-regression');
    }

    // Ensure moderators is an array
    const mods = Array.isArray(moderators) ? moderators : [moderators];

    // Extract data
    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const k = yi.length;

    // Build design matrix X (intercept + moderators)
    const X = [];
    const self = this;
    for (let i = 0; i < k; i++) {
      const row = [1];  // Intercept
      for (let j = 0; j < mods.length; j++) {
        const modName = mods[j];
        if (self.data[i][modName] === undefined) {
          throw new Error('Moderator "' + modName + '" not found in study ' + self.data[i].study);
        }
        row.push(self.data[i][modName]);
      }
      X.push(row);
    }

    const p = mods.length + 1;  // Number of parameters (intercept + moderators)

    // Estimate tau² using method of moments with residual heterogeneity
    // Start with tau² from null model
    const wi_fe = vi.map(function(v) { return 1 / v; });
    const sumW = sum(wi_fe);
    const sumW2 = sum(wi_fe.map(function(w) { return w * w; }));
    const theta_fe = sum(yi.map(function(y, i) { return wi_fe[i] * y; })) / sumW;
    const Q_total = sum(yi.map(function(y, i) { return wi_fe[i] * Math.pow(y - theta_fe, 2); }));
    const C_total = sumW - sumW2 / sumW;

    // Initial tau² estimate (DL from null model)
    let tau2 = Math.max(0, (Q_total - (k - 1)) / C_total);

    // Iterative estimation (simplified weighted least squares)
    for (let iter = 0; iter < 50; iter++) {
      // Weights
      const wi = vi.map(function(v) { return 1 / (v + tau2); });

      // Weighted least squares: beta = (X'WX)^-1 X'Wy
      const XtWX = [];
      const XtWy = [];

      // Compute X'WX
      for (let a = 0; a < p; a++) {
        XtWX[a] = [];
        for (let b = 0; b < p; b++) {
          let s = 0;
          for (let i = 0; i < k; i++) {
            s += X[i][a] * wi[i] * X[i][b];
          }
          XtWX[a][b] = s;
        }
        // Compute X'Wy
        let s = 0;
        for (let i = 0; i < k; i++) {
          s += X[i][a] * wi[i] * yi[i];
        }
        XtWy[a] = s;
      }

      // Invert XtWX (for 2x2 or simple cases)
      const beta = solveLinearSystem(XtWX, XtWy);

      // Calculate residuals
      const predicted = X.map(function(row) {
        return sum(row.map(function(x, j) { return x * beta[j]; }));
      });
      const residuals = yi.map(function(y, i) { return y - predicted[i]; });

      // Update tau² via method of moments
      const Q_res = sum(residuals.map(function(r, i) { return wi[i] * r * r; }));
      const df_res = k - p;

      // Trace term for tau² update
      let trace = 0;
      for (let i = 0; i < k; i++) {
        trace += wi[i];
      }
      // This is simplified - full REML would need more complex trace calculation
      const C_res = trace - p;  // Approximate

      const tau2_new = Math.max(0, (Q_res - df_res) / Math.max(1, C_res));

      if (Math.abs(tau2_new - tau2) < TOLERANCE.CONVERGENCE) {
        tau2 = tau2_new;
        break;
      }
      tau2 = tau2_new;
    }

    // Final estimation with converged tau²
    const wi = vi.map(function(v) { return 1 / (v + tau2); });

    // Compute X'WX and X'Wy with final weights
    const XtWX = [];
    const XtWy = [];
    for (let a = 0; a < p; a++) {
      XtWX[a] = [];
      for (let b = 0; b < p; b++) {
        let s = 0;
        for (let i = 0; i < k; i++) {
          s += X[i][a] * wi[i] * X[i][b];
        }
        XtWX[a][b] = s;
      }
      let s = 0;
      for (let i = 0; i < k; i++) {
        s += X[i][a] * wi[i] * yi[i];
      }
      XtWy[a] = s;
    }

    const beta = solveLinearSystem(XtWX, XtWy);
    const vcov = invertMatrix(XtWX);

    // Calculate residuals and Q_res
    const predicted = X.map(function(row) {
      return sum(row.map(function(x, j) { return x * beta[j]; }));
    });
    const residuals = yi.map(function(y, i) { return y - predicted[i]; });
    const Q_res = sum(residuals.map(function(r, i) { return wi[i] * r * r; }));

    // Test of moderators (Q_M)
    const Q_M = Q_total - Q_res;
    const df_M = p - 1;
    const p_M = df_M > 0 ? 1 - chiSquareCDF(Q_M, df_M) : 1;

    // Build coefficient results
    const coefficients = [];
    const labels = ['intercept'].concat(mods);

    for (let j = 0; j < p; j++) {
      const se = Math.sqrt(vcov[j][j]);
      const zval = beta[j] / se;
      const pval = 2 * (1 - normalCDF(Math.abs(zval)));

      coefficients.push({
        name: labels[j],
        estimate: beta[j],
        se: se,
        zval: zval,
        pval: pval,
        ci_lb: beta[j] - CRITICAL_VALUES.Z_95 * se,
        ci_ub: beta[j] + CRITICAL_VALUES.Z_95 * se
      });
    }

    // Residual heterogeneity
    const df_res = k - p;
    const p_res = df_res > 0 ? 1 - chiSquareCDF(Q_res, df_res) : 1;

    // I² for residual heterogeneity
    const C = sum(wi) - sum(wi.map(function(w) { return w * w; })) / sum(wi);
    const I2_res = tau2 < TOLERANCE.ZERO_CHECK ? 0 : tau2 / (tau2 + (df_res / Math.max(1, C)));

    // R² (proportion of heterogeneity explained)
    // Note: R² can be negative when moderators increase residual heterogeneity
    // This is a known property of R² in meta-regression (Raudenbush, 1994)
    const tau2_null = Math.max(0, (Q_total - (k - 1)) / C_total);
    const R2_raw = tau2_null > 0 ? (tau2_null - tau2) / tau2_null : 0;
    const R2_truncated = Math.max(0, R2_raw);
    const R2_negative = R2_raw < 0;

    return {
      coefficients: coefficients,
      tau2: tau2,
      tau: Math.sqrt(tau2),
      I2_residual: I2_res,
      R2: R2_truncated,           // Truncated at 0 (conventional reporting)
      R2_raw: R2_raw,             // Raw value (can be negative)
      R2_warning: R2_negative ?
        'Negative R² indicates moderators may increase heterogeneity or model misspecification' : null,
      QM: {
        Q: Q_M,
        df: df_M,
        p: p_M
      },
      QE: {
        Q: Q_res,
        df: df_res,
        p: p_res
      },
      k: k,
      p: p,
      moderators: mods
    };
  };

  /**
   * Simple linear system solver (Gaussian elimination)
   * For small matrices used in meta-regression
   * Includes pivot tolerance check to prevent division by zero
   */
  function solveLinearSystem(A, b) {
    const n = A.length;
    const aug = A.map(function(row, i) { return row.slice().concat([b[i]]); });

    // Forward elimination with partial pivoting
    for (let i = 0; i < n; i++) {
      // Find pivot (row with largest absolute value in column i)
      let maxRow = i;
      let maxVal = Math.abs(aug[i][i]);
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(aug[k][i]) > maxVal) {
          maxVal = Math.abs(aug[k][i]);
          maxRow = k;
        }
      }

      // Check for near-singular matrix
      if (maxVal < TOLERANCE.PIVOT) {
        throw new Error('Matrix is singular or near-singular (pivot < ' + TOLERANCE.PIVOT + ')');
      }

      // Swap rows
      const tmp = aug[i];
      aug[i] = aug[maxRow];
      aug[maxRow] = tmp;

      // Eliminate column
      const pivot = aug[i][i];
      for (let k = i + 1; k < n; k++) {
        const c = aug[k][i] / pivot;
        for (let j = i; j <= n; j++) {
          aug[k][j] -= c * aug[i][j];
        }
      }
    }

    // Back substitution
    const x = new Array(n);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = aug[i][n];
      for (let j = i + 1; j < n; j++) {
        x[i] -= aug[i][j] * x[j];
      }
      // Check for near-zero diagonal
      if (Math.abs(aug[i][i]) < TOLERANCE.PIVOT) {
        throw new Error('Division by near-zero in back substitution');
      }
      x[i] /= aug[i][i];
    }

    return x;
  }

  /**
   * Matrix inversion (for variance-covariance matrix)
   * Uses Gauss-Jordan elimination with pivot tolerance checks
   */
  function invertMatrix(A) {
    const n = A.length;
    const aug = A.map(function(row, i) {
      const newRow = row.slice();
      for (let j = 0; j < n; j++) {
        newRow.push(i === j ? 1 : 0);
      }
      return newRow;
    });

    // Forward elimination with partial pivoting
    for (let i = 0; i < n; i++) {
      // Find best pivot
      let maxRow = i;
      let maxVal = Math.abs(aug[i][i]);
      for (let k = i + 1; k < n; k++) {
        if (Math.abs(aug[k][i]) > maxVal) {
          maxVal = Math.abs(aug[k][i]);
          maxRow = k;
        }
      }

      // Check for near-singular matrix
      if (maxVal < TOLERANCE.PIVOT) {
        throw new Error('Matrix is singular or near-singular, cannot invert');
      }

      // Swap rows
      const tmp = aug[i];
      aug[i] = aug[maxRow];
      aug[maxRow] = tmp;

      const pivot = aug[i][i];
      for (let j = 0; j < 2 * n; j++) {
        aug[i][j] /= pivot;
      }

      for (let k = 0; k < n; k++) {
        if (k !== i) {
          const c = aug[k][i];
          for (let j = 0; j < 2 * n; j++) {
            aug[k][j] -= c * aug[i][j];
          }
        }
      }
    }

    // Extract inverse
    return aug.map(function(row) { return row.slice(n); });
  }

  /**
   * Egger's regression test for publication bias
   * @returns {Object} Test results
   */
  MetaAnalysis.prototype.eggerTest = function() {
    if (this.data.length < 3) {
      return { intercept: NaN, se: NaN, tval: NaN, pval: NaN };
    }

    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const result = PublicationBias.egger(yi, vi);

    return {
      intercept: result.intercept,
      se: result.se,
      tval: result.t,
      pval: result.p,
      slope: result.slope
    };
  };

  /**
   * Begg's rank correlation test for publication bias
   * @returns {Object} Test results
   */
  MetaAnalysis.prototype.beggTest = function() {
    if (this.data.length < 3) {
      return { tau: NaN, pval: NaN };
    }

    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const result = PublicationBias.begg(yi, vi);

    return {
      tau: result.tau,
      z: result.z,
      pval: result.p
    };
  };

  /**
   * Trim-and-Fill method for publication bias
   * @param {String} side - Side to fill ('left', 'right', or 'auto')
   * @returns {Object} Adjusted results
   */
  MetaAnalysis.prototype.trimAndFill = function(side) {
    side = side || 'auto';
    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const result = PublicationBias.trimFill(yi, vi, this.model.tau2, side);

    return {
      k0: result.k0,
      estimate: result.estimate,
      se: result.se,
      ci_lb: result.estimate - CRITICAL_VALUES.Z_95 * result.se,
      ci_ub: result.estimate + CRITICAL_VALUES.Z_95 * result.se,
      side: result.side,
      filled: result.imputed
    };
  };

  /**
   * Fail-Safe N analysis (File Drawer Analysis)
   * Calculates number of null studies needed to nullify results
   * Three methods: Rosenthal, Orwin, Rosenberg
   *
   * @param {Object} options - Configuration options
   * @param {Number} options.alpha - Significance level (default: 0.05)
   * @param {Number} options.targetES - Target effect size for Orwin method (default: 0.2)
   * @returns {Object} Fail-safe N results for all three methods
   */
  MetaAnalysis.prototype.failSafeN = function(options) {
    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    return PublicationBias.failSafeN(yi, vi, options);
  };

  /**
   * Leave-one-out sensitivity analysis
   * @param {Object} options - Analysis options
   * @param {String} options.method - Tau² estimation method
   * @param {Function} options.onProgress - Progress callback(current, total, study)
   * @returns {Array} Results for each omitted study
   */
  MetaAnalysis.prototype.leaveOneOut = function(options) {
    options = options || {};
    const method = options.method || this.options.method || 'REML';
    const onProgress = options.onProgress || this.options.onProgress;
    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const self = this;
    const k = yi.length;
    const results = [];

    for (let i = 0; i < k; i++) {
      // Progress callback
      if (onProgress) {
        onProgress(i + 1, k, self.data[i].study);
      }

      const yi_loo = yi.filter(function(_, j) { return j !== i; });
      const vi_loo = vi.filter(function(_, j) { return j !== i; });

      const tau2_result = Tau2Estimators[method](yi_loo, vi_loo);
      const pooled = PooledEstimate.calculate(yi_loo, vi_loo, tau2_result.tau2);
      const het = Heterogeneity.calculate(yi_loo, vi_loo, tau2_result.tau2);

      results.push({
        omitted: self.data[i].study,
        estimate: pooled.estimate,
        se: pooled.se,
        ci_lb: pooled.ci_lb,
        ci_ub: pooled.ci_ub,
        tau2: tau2_result.tau2,
        I2: het.I2
      });
    }

    return results;
  };

  /**
   * Cumulative meta-analysis
   * @param {Object} options - Analysis options
   * @param {String} options.sortBy - Field to sort studies by
   * @param {Function} options.onProgress - Progress callback(current, total, study)
   * @returns {Array} Cumulative results
   */
  MetaAnalysis.prototype.cumulativeMeta = function(options) {
    // Support legacy string argument for sortBy
    if (typeof options === 'string') {
      options = { sortBy: options };
    }
    options = options || {};

    const sortBy = options.sortBy;
    const onProgress = options.onProgress || this.options.onProgress;
    const sortedData = this.data.slice();

    if (sortBy && this.data[0][sortBy] !== undefined) {
      sortedData.sort(function(a, b) { return a[sortBy] - b[sortBy]; });
    }

    const results = [];
    const tempMA = new MetaAnalysis(this.options);
    const self = this;
    const total = sortedData.length;

    for (let i = 1; i <= total; i++) {
      // Progress callback
      if (onProgress) {
        onProgress(i, total, sortedData[i - 1].study);
      }

      tempMA.data = sortedData.slice(0, i).map(function(d) {
        return { yi: d.yi, vi: d.vi, study: d.study };
      });
      if (i >= 2) {
        const model = tempMA.runRandomEffectsModel({ method: self.options.method });
        results.push({
          k: i,
          lastStudy: sortedData[i - 1].study,
          estimate: model.estimate,
          se: model.se,
          ci_lb: model.ci_lb,
          ci_ub: model.ci_ub,
          tau2: model.tau2,
          I2: model.I2
        });
      }
    }

    return results;
  };

  /**
   * Influence diagnostics
   * @returns {Array} Diagnostic measures for each study
   */
  MetaAnalysis.prototype.influenceDiagnostics = function() {
    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const yi = this.data.map(function(d) { return d.yi; });
    const vi = this.data.map(function(d) { return d.vi; });
    const results = Sensitivity.influence(yi, vi, this.model.tau2);
    const self = this;

    return results.map(function(r, i) {
      return {
        study: self.data[i].study,
        residual: r.residual,
        rstudent: r.rstudent,
        hat: r.hat,
        cooksd: r.cooksd,
        dfbetas: r.dfbetas
      };
    });
  };

  /**
   * Format result for display
   * @param {Object} options - Formatting options
   * @returns {Object} Formatted results
   */
  MetaAnalysis.prototype.formatResult = function(options) {
    options = options || {};
    const decimals = options.decimalPlaces !== undefined ? options.decimalPlaces : 3;
    const transform = options.transformBack !== undefined ? options.transformBack : false;
    const measure = this.options.measure;

    if (!this.model) {
      return null;
    }

    let est = this.model.estimate;
    let ci_lb = this.model.ci_lb;
    let ci_ub = this.model.ci_ub;

    // Back-transform for ratio measures
    if (transform && ['OR', 'RR', 'HR'].indexOf(measure) >= 0) {
      est = Math.exp(est);
      ci_lb = Math.exp(ci_lb);
      ci_ub = Math.exp(ci_ub);
    }

    return {
      estimate: est.toFixed(decimals),
      ci: '[' + ci_lb.toFixed(decimals) + ', ' + ci_ub.toFixed(decimals) + ']',
      formatted: measure + ' = ' + est.toFixed(decimals) + ' (95% CI: ' + ci_lb.toFixed(decimals) + ', ' + ci_ub.toFixed(decimals) + ')',
      pval: this.model.pval < 0.001 ? '< 0.001' : this.model.pval.toFixed(3),
      I2: (this.model.I2 * 100).toFixed(1) + '%',
      tau2: this.model.tau2.toFixed(4)
    };
  };

  /**
   * Get forest plot data
   * @returns {Array} Data ready for forest plot visualization
   */
  MetaAnalysis.prototype.getForestPlotData = function() {
    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const measure = this.options.measure;
    const isRatio = ['OR', 'RR', 'HR'].indexOf(measure) >= 0;
    const self = this;

    const studies = this.data.map(function(d, i) {
      let est = d.yi;
      let ci_lb = d.yi - CRITICAL_VALUES.Z_95 * d.se;
      let ci_ub = d.yi + CRITICAL_VALUES.Z_95 * d.se;

      if (isRatio) {
        est = Math.exp(est);
        ci_lb = Math.exp(ci_lb);
        ci_ub = Math.exp(ci_ub);
      }

      return {
        study: d.study,
        y: i,
        estimate: est,
        ci_lb: ci_lb,
        ci_ub: ci_ub,
        weight: d.weight,
        isPooled: false
      };
    });

    // Add pooled estimate
    let pooledEst = this.model.estimate;
    let pooledLb = this.model.ci_lb;
    let pooledUb = this.model.ci_ub;

    if (isRatio) {
      pooledEst = Math.exp(pooledEst);
      pooledLb = Math.exp(pooledLb);
      pooledUb = Math.exp(pooledUb);
    }

    // Calculate prediction interval if tau2 > 0
    let pi_lb = null;
    let pi_ub = null;
    if (this.model.tau2 > 0 && this.model.k >= 3) {
      // Use t-distribution for PI (k-2 df)
      const df = this.model.k - 2;
      const t_crit = tQuantile(0.975, df);
      const pi_se = Math.sqrt(this.model.se * this.model.se + this.model.tau2);
      pi_lb = this.model.estimate - t_crit * pi_se;
      pi_ub = this.model.estimate + t_crit * pi_se;

      if (isRatio) {
        pi_lb = Math.exp(pi_lb);
        pi_ub = Math.exp(pi_ub);
      }
    }

    studies.push({
      study: 'Pooled',
      y: this.data.length,
      estimate: pooledEst,
      ci_lb: pooledLb,
      ci_ub: pooledUb,
      weight: 100,
      isPooled: true,
      // Prediction interval (for random effects)
      pi_lb: pi_lb,
      pi_ub: pi_ub,
      hasPredictionInterval: pi_lb !== null
    });

    // Return object with studies and metadata for easier plotting
    return {
      studies: studies,
      pooled: {
        estimate: pooledEst,
        ci: [pooledLb, pooledUb],
        pi: pi_lb !== null ? [pi_lb, pi_ub] : null
      },
      meta: {
        k: this.model.k,
        measure: measure,
        isRatio: isRatio,
        nullValue: isRatio ? 1 : 0,
        tau2: this.model.tau2,
        I2: this.model.I2
      }
    };
  };

  /**
   * Get funnel plot data
   * @returns {Object} Data for funnel plot visualization
   */
  MetaAnalysis.prototype.getFunnelPlotData = function() {
    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const measure = this.options.measure;
    const isRatio = ['OR', 'RR', 'HR'].indexOf(measure) >= 0;
    const pooledEstimate = this.model.estimate;

    const points = this.data.map(function(d) {
      return {
        study: d.study,
        x: isRatio ? Math.exp(d.yi) : d.yi,  // Effect size
        y: d.se,                              // Standard error (precision)
        yi: d.yi,
        se: d.se
      };
    });

    // Calculate funnel boundaries (95% CI pseudo-confidence regions)
    const seMax = Math.max.apply(null, this.data.map(function(d) { return d.se; }));
    const seMin = 0;
    const seRange = [];
    for (let s = 0; s <= 1; s += 0.02) {
      seRange.push(seMin + s * seMax);
    }

    const lowerBound = seRange.map(function(se) {
      const y = pooledEstimate - CRITICAL_VALUES.Z_95 * se;
      return { se: se, x: isRatio ? Math.exp(y) : y };
    });

    const upperBound = seRange.map(function(se) {
      const y = pooledEstimate + CRITICAL_VALUES.Z_95 * se;
      return { se: se, x: isRatio ? Math.exp(y) : y };
    });

    return {
      points: points,
      pooledEstimate: isRatio ? Math.exp(pooledEstimate) : pooledEstimate,
      lowerBound: lowerBound,
      upperBound: upperBound,
      seMax: seMax,
      measure: measure
    };
  };

  /**
   * Subgroup analysis
   * @param {String} groupVar - Variable name to group by
   * @returns {Object} Results by subgroup and overall
   */
  MetaAnalysis.prototype.subgroupAnalysis = function(groupVar) {
    if (!groupVar) {
      throw new Error('Group variable is required for subgroup analysis');
    }

    // Check if group variable exists in data
    if (!this.data[0] || this.data[0][groupVar] === undefined) {
      throw new Error('Group variable "' + groupVar + '" not found in data');
    }

    // Group studies
    const groups = {};
    this.data.forEach(function(d) {
      const group = d[groupVar];
      if (!groups[group]) {
        groups[group] = [];
      }
      groups[group].push(d);
    });

    const groupNames = Object.keys(groups);
    const subgroupResults = {};
    const self = this;

    // Run meta-analysis for each subgroup
    groupNames.forEach(function(groupName) {
      const groupData = groups[groupName];
      if (groupData.length >= 2) {
        const tempMA = new MetaAnalysis(self.options);
        tempMA.data = groupData;
        subgroupResults[groupName] = tempMA.runRandomEffectsModel({
          method: self.options.method
        });
        subgroupResults[groupName].k = groupData.length;
      } else {
        // Single study - just report the study values
        subgroupResults[groupName] = {
          estimate: groupData[0].yi,
          se: groupData[0].se,
          ci_lb: groupData[0].yi - CRITICAL_VALUES.Z_95 * groupData[0].se,
          ci_ub: groupData[0].yi + CRITICAL_VALUES.Z_95 * groupData[0].se,
          k: 1,
          tau2: 0,
          I2: 0
        };
      }
    });

    // Test for subgroup differences (Q-between)
    const overallModel = this.model || this.runRandomEffectsModel();
    let Q_within = 0;
    let df_within = 0;

    groupNames.forEach(function(groupName) {
      if (subgroupResults[groupName].Q) {
        Q_within += subgroupResults[groupName].Q;
        df_within += subgroupResults[groupName].Q_df || 0;
      }
    });

    const Q_total = overallModel.Q;
    const Q_between = Q_total - Q_within;
    const df_between = groupNames.length - 1;
    const p_between = df_between > 0 ? 1 - chiSquareCDF(Q_between, df_between) : 1;

    return {
      subgroups: subgroupResults,
      test: {
        Q_between: Q_between,
        df: df_between,
        p: p_between,
        heterogeneityDifference: Q_between > 0 ? 'Significant' : 'Not significant'
      },
      overall: overallModel
    };
  };

  /**
   * Export data to JSON format
   * @param {Object} options - Export options
   * @returns {String} JSON string
   */
  MetaAnalysis.prototype.toJSON = function(options) {
    options = options || {};
    const includeModel = options.includeModel !== false;
    const includeStudies = options.includeStudies !== false;

    const result = {
      meta: {
        version: '1.2.0',
        measure: this.options.measure,
        method: this.options.method,
        generatedAt: new Date().toISOString()
      }
    };

    if (includeStudies) {
      result.studies = this.data.map(function(d) {
        return {
          study: d.study,
          yi: d.yi,
          vi: d.vi,
          se: d.se,
          weight: d.weight
        };
      });
    }

    if (includeModel && this.model) {
      result.model = {
        estimate: this.model.estimate,
        se: this.model.se,
        ci_lb: this.model.ci_lb,
        ci_ub: this.model.ci_ub,
        pval: this.model.pval,
        tau2: this.model.tau2,
        I2: this.model.I2,
        Q: this.model.Q,
        Q_pval: this.model.Q_pval,
        k: this.model.k
      };
    }

    return JSON.stringify(result, null, options.pretty ? 2 : 0);
  };

  /**
   * Export data to CSV format
   * @param {Object} options - Export options
   * @param {String} options.delimiter - Column delimiter (default: ',')
   * @param {Boolean} options.includeHeader - Include header row (default: true)
   * @param {Boolean} options.includeModel - Include model results section (default: false)
   * @param {Boolean} options.includeRoB - Include Risk of Bias data (default: false)
   * @returns {String} CSV string
   */
  MetaAnalysis.prototype.toCSV = function(options) {
    options = options || {};
    const delimiter = options.delimiter || ',';
    const includeHeader = options.includeHeader !== false;
    const includeModel = options.includeModel || false;
    const includeRoB = options.includeRoB || false;

    const rows = [];

    // Study data header
    if (includeHeader) {
      let header = ['study', 'yi', 'vi', 'se', 'ci_lb', 'ci_ub', 'weight'];
      if (includeRoB) {
        header.push('rob_overall');
      }
      rows.push(header.join(delimiter));
    }

    // Study data
    this.data.forEach(function(d) {
      const ci_lb = d.yi - CRITICAL_VALUES.Z_95 * d.se;
      const ci_ub = d.yi + CRITICAL_VALUES.Z_95 * d.se;
      let row = [
        '"' + (d.study || '').replace(/"/g, '""') + '"',
        d.yi.toFixed(6),
        d.vi.toFixed(6),
        d.se.toFixed(6),
        ci_lb.toFixed(6),
        ci_ub.toFixed(6),
        (d.weight || 0).toFixed(2)
      ];
      if (includeRoB) {
        row.push(d.rob && d.rob.overall ? d.rob.overall : '');
      }
      rows.push(row.join(delimiter));
    });

    // Model results section
    if (includeModel && this.model) {
      rows.push('');  // Empty row separator
      rows.push('# Model Results');
      rows.push(['Parameter', 'Value'].join(delimiter));
      rows.push(['Method', this.model.method].join(delimiter));
      rows.push(['k (studies)', this.model.k].join(delimiter));
      rows.push(['Pooled estimate', this.model.estimate.toFixed(6)].join(delimiter));
      rows.push(['SE', this.model.se.toFixed(6)].join(delimiter));
      rows.push(['95% CI lower', this.model.ci_lb.toFixed(6)].join(delimiter));
      rows.push(['95% CI upper', this.model.ci_ub.toFixed(6)].join(delimiter));
      rows.push(['z-value', this.model.zval.toFixed(4)].join(delimiter));
      rows.push(['p-value', this.model.pval < 0.0001 ? '< 0.0001' : this.model.pval.toFixed(6)].join(delimiter));
      rows.push(['tau²', this.model.tau2.toFixed(6)].join(delimiter));
      rows.push(['tau² 95% CI', '[' + this.model.tau2_ci_lb.toFixed(4) + ', ' + this.model.tau2_ci_ub.toFixed(4) + ']'].join(delimiter));
      rows.push(['I²', (this.model.I2 * 100).toFixed(2) + '%'].join(delimiter));
      rows.push(['Q statistic', this.model.Q.toFixed(4)].join(delimiter));
      rows.push(['Q p-value', this.model.Q_pval.toFixed(6)].join(delimiter));

      // HKSJ results if available
      if (this.model.hksj) {
        rows.push(['HKSJ CI lower', this.model.hksj_ci_lb.toFixed(6)].join(delimiter));
        rows.push(['HKSJ CI upper', this.model.hksj_ci_ub.toFixed(6)].join(delimiter));
        rows.push(['HKSJ p-value', this.model.hksj_pval.toFixed(6)].join(delimiter));
      }

      // Prediction interval if available
      if (this.model.predictionInterval) {
        rows.push(['PI lower', this.model.pi_lb.toFixed(6)].join(delimiter));
        rows.push(['PI upper', this.model.pi_ub.toFixed(6)].join(delimiter));
      }
    }

    return rows.join('\n');
  };

  /**
   * Validate data with detailed error reporting
   * @returns {Object} Validation result with errors and warnings
   */
  MetaAnalysis.prototype.validate = function() {
    const errors = [];
    const warnings = [];
    const self = this;

    // Check if data exists
    if (!this.data || this.data.length === 0) {
      errors.push({
        code: 'NO_DATA',
        message: 'No data loaded. Use calculateEffectSizes() first.'
      });
      return { valid: false, errors: errors, warnings: warnings };
    }

    // Check minimum studies
    if (this.data.length < 2) {
      errors.push({
        code: 'INSUFFICIENT_STUDIES',
        message: 'At least 2 studies are required for meta-analysis.'
      });
    }

    // Validate each study
    this.data.forEach(function(d, i) {
      const studyId = d.study || ('Study ' + (i + 1));

      // Check effect size
      if (d.yi === undefined || d.yi === null || !isFinite(d.yi)) {
        errors.push({
          code: 'INVALID_EFFECT_SIZE',
          study: studyId,
          message: 'Invalid effect size (yi) for study "' + studyId + '"'
        });
      }

      // Check variance
      if (d.vi === undefined || d.vi === null || !isFinite(d.vi)) {
        errors.push({
          code: 'INVALID_VARIANCE',
          study: studyId,
          message: 'Invalid variance (vi) for study "' + studyId + '"'
        });
      } else if (d.vi <= 0) {
        errors.push({
          code: 'NON_POSITIVE_VARIANCE',
          study: studyId,
          message: 'Variance must be positive for study "' + studyId + '"'
        });
      }

      // Warn about extreme effect sizes
      if (d.yi !== undefined && Math.abs(d.yi) > 10) {
        warnings.push({
          code: 'EXTREME_EFFECT',
          study: studyId,
          message: 'Effect size seems extreme for study "' + studyId + '" (yi = ' + d.yi.toFixed(3) + ')'
        });
      }

      // Warn about very small variance
      if (d.vi !== undefined && d.vi < 0.0001) {
        warnings.push({
          code: 'SMALL_VARIANCE',
          study: studyId,
          message: 'Very small variance for study "' + studyId + '", which may dominate the pooled estimate'
        });
      }
    });

    return {
      valid: errors.length === 0,
      errors: errors,
      warnings: warnings,
      summary: {
        studies: this.data.length,
        errorCount: errors.length,
        warningCount: warnings.length
      }
    };
  };

  // ============================================
  // CONVENIENCE METHODS
  // ============================================

  /**
   * Get a comprehensive summary of the analysis
   * @returns {Object} Summary with key statistics
   */
  MetaAnalysis.prototype.getSummary = function() {
    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const measure = this.options.measure;
    const isRatio = ['OR', 'RR', 'HR'].indexOf(measure) >= 0;

    let estimate = this.model.estimate;
    let ci_lb = this.model.ci_lb;
    let ci_ub = this.model.ci_ub;

    if (isRatio) {
      estimate = Math.exp(estimate);
      ci_lb = Math.exp(ci_lb);
      ci_ub = Math.exp(ci_ub);
    }

    return {
      // Effect estimate
      estimate: estimate,
      ci: [ci_lb, ci_ub],
      pval: this.model.pval,
      significant: this.model.pval < 0.05,

      // Heterogeneity
      heterogeneity: {
        tau2: this.model.tau2,
        tau2_ci: this.model.tau2_ci_lb !== undefined ? [this.model.tau2_ci_lb, this.model.tau2_ci_ub] : null,
        I2: this.model.I2,
        I2_ci: this.model.I2_ci_lb !== undefined ? [this.model.I2_ci_lb, this.model.I2_ci_ub] : null,
        I2_percent: (this.model.I2 * 100).toFixed(1) + '%',
        level: this.getHeterogeneityLevel(),
        Q: this.model.Q,
        Q_pval: this.model.Q_pval
      },

      // Study info
      studies: this.model.k,
      method: this.model.method,
      measure: measure,

      // Formatted strings
      formatted: {
        estimate: isRatio ? estimate.toFixed(3) : estimate.toFixed(4),
        ci: '[' + ci_lb.toFixed(3) + ', ' + ci_ub.toFixed(3) + ']',
        pval: this.model.pval < 0.001 ? '< 0.001' : this.model.pval.toFixed(4),
        I2: (this.model.I2 * 100).toFixed(1) + '%'
      }
    };
  };

  /**
   * GRADE (Grading of Recommendations Assessment, Development and Evaluation) summary
   * Assesses certainty of evidence across 5 domains
   *
   * @param {Object} options - GRADE assessment options
   * @param {String} options.riskOfBias - Overall RoB: 'low', 'some', 'high' (default: null)
   * @param {String} options.indirectness - 'none', 'some', 'serious' (default: null)
   * @param {Number} options.optimalInformationSize - OIS for imprecision (default: calculated)
   * @param {Boolean} options.publicationBiasDetected - From Egger/Begg tests (default: auto)
   * @returns {Object} GRADE assessment with certainty rating
   */
  MetaAnalysis.prototype.getGRADESummary = function(options) {
    options = options || {};

    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const downgrades = [];
    let certaintyScore = 4;  // Start at HIGH (4 = high, 3 = moderate, 2 = low, 1 = very low)

    // Domain 1: Risk of Bias
    let robAssessment = options.riskOfBias || null;
    if (robAssessment === 'high') {
      downgrades.push({ domain: 'Risk of bias', reason: 'High risk of bias', levels: 2 });
      certaintyScore -= 2;
    } else if (robAssessment === 'some') {
      downgrades.push({ domain: 'Risk of bias', reason: 'Some concerns', levels: 1 });
      certaintyScore -= 1;
    }

    // Domain 2: Inconsistency (based on I²)
    const I2 = this.model.I2;
    let inconsistencyAssessment = 'none';
    if (I2 > 0.75) {
      downgrades.push({ domain: 'Inconsistency', reason: 'I² > 75% (considerable)', levels: 2 });
      certaintyScore -= 2;
      inconsistencyAssessment = 'serious';
    } else if (I2 > 0.50) {
      downgrades.push({ domain: 'Inconsistency', reason: 'I² 50-75% (substantial)', levels: 1 });
      certaintyScore -= 1;
      inconsistencyAssessment = 'some';
    } else if (I2 > 0.25 && this.model.Q_pval < 0.10) {
      downgrades.push({ domain: 'Inconsistency', reason: 'I² > 25% with significant Q', levels: 1 });
      certaintyScore -= 1;
      inconsistencyAssessment = 'some';
    }

    // Domain 3: Indirectness (user-provided)
    const indirectness = options.indirectness || null;
    if (indirectness === 'serious') {
      downgrades.push({ domain: 'Indirectness', reason: 'Serious concerns', levels: 2 });
      certaintyScore -= 2;
    } else if (indirectness === 'some') {
      downgrades.push({ domain: 'Indirectness', reason: 'Some concerns', levels: 1 });
      certaintyScore -= 1;
    }

    // Domain 4: Imprecision
    // Check if CI crosses clinically important thresholds and/or OIS
    const measure = this.options.measure;
    const isRatio = ['OR', 'RR', 'HR'].indexOf(measure) >= 0;
    let ci_lb = this.model.ci_lb;
    let ci_ub = this.model.ci_ub;
    let estimate = this.model.estimate;

    if (isRatio) {
      ci_lb = Math.exp(ci_lb);
      ci_ub = Math.exp(ci_ub);
      estimate = Math.exp(estimate);
    }

    // Check for imprecision: CI crosses null or threshold, or wide CI
    const nullValue = isRatio ? 1 : 0;
    const crossesNull = (ci_lb <= nullValue && ci_ub >= nullValue);
    const relativeCI = isRatio ? (ci_ub / ci_lb) : (ci_ub - ci_lb);

    // Calculate total sample size
    const totalN = this.data.reduce(function(sum, d) {
      return sum + (d.n_t || 0) + (d.n_c || 0);
    }, 0);

    // Optimal Information Size (OIS) check - simplified rule
    const ois = options.optimalInformationSize || null;
    let imprecisionAssessment = 'none';

    if (crossesNull && relativeCI > (isRatio ? 2.5 : 0.5)) {
      downgrades.push({ domain: 'Imprecision', reason: 'Wide CI crossing null', levels: 2 });
      certaintyScore -= 2;
      imprecisionAssessment = 'serious';
    } else if (crossesNull) {
      downgrades.push({ domain: 'Imprecision', reason: 'CI crosses null value', levels: 1 });
      certaintyScore -= 1;
      imprecisionAssessment = 'some';
    } else if (ois && totalN < ois) {
      downgrades.push({ domain: 'Imprecision', reason: 'Sample size below OIS', levels: 1 });
      certaintyScore -= 1;
      imprecisionAssessment = 'some';
    }

    // Domain 5: Publication Bias
    let pubBiasDetected = options.publicationBiasDetected;
    if (pubBiasDetected === undefined) {
      // Auto-detect from Egger test if k >= 10
      if (this.model.k >= 10) {
        const egger = this.eggerTest();
        pubBiasDetected = egger.pval < 0.10;
      } else {
        pubBiasDetected = null;  // Cannot assess with few studies
      }
    }

    if (pubBiasDetected === true) {
      downgrades.push({ domain: 'Publication bias', reason: 'Funnel plot asymmetry detected', levels: 1 });
      certaintyScore -= 1;
    }

    // Determine final certainty level
    certaintyScore = Math.max(1, certaintyScore);  // Floor at 1
    const certaintyLevels = ['', 'Very Low', 'Low', 'Moderate', 'High'];
    const certaintyLevel = certaintyLevels[certaintyScore];

    // GRADE symbols
    const gradeSymbols = {
      'High': '⊕⊕⊕⊕',
      'Moderate': '⊕⊕⊕○',
      'Low': '⊕⊕○○',
      'Very Low': '⊕○○○'
    };

    return {
      certainty: certaintyLevel,
      certaintyScore: certaintyScore,
      symbol: gradeSymbols[certaintyLevel],
      downgrades: downgrades,
      domains: {
        riskOfBias: robAssessment || 'not assessed',
        inconsistency: inconsistencyAssessment,
        indirectness: indirectness || 'not assessed',
        imprecision: imprecisionAssessment,
        publicationBias: pubBiasDetected === null ? 'not assessable (k < 10)' :
          (pubBiasDetected ? 'detected' : 'not detected')
      },
      evidence: {
        k: this.model.k,
        totalN: totalN || 'N/A',
        estimate: estimate,
        ci: [ci_lb, ci_ub],
        I2: I2
      },
      interpretation: this._getGRADEInterpretation(certaintyLevel)
    };
  };

  /**
   * Get GRADE interpretation text
   * @private
   */
  MetaAnalysis.prototype._getGRADEInterpretation = function(level) {
    const interpretations = {
      'High': 'We are very confident that the true effect lies close to the estimate.',
      'Moderate': 'We are moderately confident; the true effect is likely close to the estimate but may be substantially different.',
      'Low': 'Our confidence is limited; the true effect may be substantially different from the estimate.',
      'Very Low': 'We have very little confidence; the true effect is likely substantially different from the estimate.'
    };
    return interpretations[level] || '';
  };

  /**
   * Risk of Bias (RoB) domains based on Cochrane RoB 2 tool
   * @constant
   */
  MetaAnalysis.ROB_DOMAINS = Object.freeze({
    // RoB 2 domains for randomized trials
    ROB2: [
      'randomization',      // D1: Risk of bias arising from the randomization process
      'deviations',         // D2: Risk of bias due to deviations from intended interventions
      'missingData',        // D3: Risk of bias due to missing outcome data
      'measurement',        // D4: Risk of bias in measurement of the outcome
      'selection'           // D5: Risk of bias in selection of the reported result
    ],
    // ROBINS-I domains for non-randomized studies
    ROBINS_I: [
      'confounding',        // Bias due to confounding
      'selection',          // Bias in selection of participants
      'classification',     // Bias in classification of interventions
      'deviations',         // Bias due to deviations from intended interventions
      'missingData',        // Bias due to missing data
      'measurement',        // Bias in measurement of outcomes
      'reportedResult'      // Bias in selection of the reported result
    ],
    // Judgment levels
    LEVELS: ['low', 'some', 'high', 'unclear']
  });

  /**
   * Set Risk of Bias assessment for a study
   * Supports both RoB 2 and ROBINS-I frameworks
   *
   * @param {Number|String} studyIdentifier - Index or study name
   * @param {Object} robData - RoB assessment data
   * @param {String} robData.overall - Overall judgment: 'low', 'some', 'high', 'unclear'
   * @param {Object} robData.domains - Domain-specific judgments
   * @param {String} robData.framework - 'ROB2' or 'ROBINS_I' (default: 'ROB2')
   * @returns {MetaAnalysis} this for method chaining
   */
  MetaAnalysis.prototype.setRiskOfBias = function(studyIdentifier, robData) {
    const idx = typeof studyIdentifier === 'number' ? studyIdentifier :
      this.data.findIndex(function(d) { return d.study === studyIdentifier; });

    if (idx < 0 || idx >= this.data.length) {
      throw new Error('Study not found: ' + studyIdentifier);
    }

    // Initialize RoB structure if not present
    if (!this.data[idx].rob) {
      this.data[idx].rob = {
        overall: null,
        domains: {},
        framework: 'ROB2',
        notes: ''
      };
    }

    // Update with provided data
    if (robData.overall) {
      this.data[idx].rob.overall = robData.overall;
    }
    if (robData.domains) {
      Object.assign(this.data[idx].rob.domains, robData.domains);
    }
    if (robData.framework) {
      this.data[idx].rob.framework = robData.framework;
    }
    if (robData.notes) {
      this.data[idx].rob.notes = robData.notes;
    }

    return this;  // Method chaining
  };

  /**
   * Get Risk of Bias summary across all studies
   *
   * @returns {Object} RoB summary with counts and percentages
   */
  MetaAnalysis.prototype.getRiskOfBiasSummary = function() {
    const self = this;
    const summary = {
      overall: { low: 0, some: 0, high: 0, unclear: 0, notAssessed: 0 },
      domains: {},
      studies: [],
      percentLowRisk: 0
    };

    // Count overall judgments
    this.data.forEach(function(d, i) {
      const studyRob = {
        study: d.study,
        overall: null,
        domains: {}
      };

      if (d.rob && d.rob.overall) {
        summary.overall[d.rob.overall]++;
        studyRob.overall = d.rob.overall;
        studyRob.domains = d.rob.domains || {};
      } else {
        summary.overall.notAssessed++;
      }

      summary.studies.push(studyRob);
    });

    // Calculate percentage of low risk studies
    const assessed = this.data.length - summary.overall.notAssessed;
    if (assessed > 0) {
      summary.percentLowRisk = (summary.overall.low / assessed) * 100;
    }

    // Aggregate domain-level assessments
    const frameworks = ['ROB2', 'ROBINS_I'];
    frameworks.forEach(function(fw) {
      const domains = MetaAnalysis.ROB_DOMAINS[fw];
      domains.forEach(function(domain) {
        if (!summary.domains[domain]) {
          summary.domains[domain] = { low: 0, some: 0, high: 0, unclear: 0 };
        }
        self.data.forEach(function(d) {
          if (d.rob && d.rob.domains && d.rob.domains[domain]) {
            summary.domains[domain][d.rob.domains[domain]]++;
          }
        });
      });
    });

    // Overall assessment for GRADE
    summary.gradeAssessment = summary.overall.high > 0 ? 'high' :
      (summary.overall.some > summary.overall.low ? 'some' : 'low');

    return summary;
  };

  /**
   * Get the number of studies in the analysis
   * @returns {Number} Number of studies
   */
  MetaAnalysis.prototype.getStudyCount = function() {
    return this.data.length;
  };

  /**
   * Check if the pooled effect is statistically significant
   * @param {Number} alpha - Significance level (default: 0.05)
   * @returns {Boolean} True if significant
   */
  MetaAnalysis.prototype.isSignificant = function(alpha) {
    alpha = alpha || 0.05;
    if (!this.model) {
      this.runRandomEffectsModel();
    }
    return this.model.pval < alpha;
  };

  /**
   * Get heterogeneity level interpretation
   * Based on Higgins & Thompson (2002) guidelines
   * @returns {String} Interpretation of I² level
   */
  MetaAnalysis.prototype.getHeterogeneityLevel = function() {
    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const I2 = this.model.I2;

    if (I2 < 0.25) {
      return 'low';
    } else if (I2 < 0.50) {
      return 'moderate';
    } else if (I2 < 0.75) {
      return 'substantial';
    } else {
      return 'considerable';
    }
  };

  /**
   * Get effect size direction based on statistical significance
   * Uses CI to determine if effect is significantly different from null
   *
   * @param {Object} options - Direction options
   * @param {Number} options.alpha - Significance level (default: 0.05)
   * @param {Boolean} options.favorLower - Whether lower values are favorable (default: true for OR/RR/HR)
   * @returns {Object} Effect direction with significance info
   */
  MetaAnalysis.prototype.getEffectDirection = function(options) {
    options = options || {};

    if (!this.model) {
      this.runRandomEffectsModel();
    }

    const measure = this.options.measure;
    const isRatio = ['OR', 'RR', 'HR'].indexOf(measure) >= 0;
    const estimate = this.model.estimate;
    const ci_lb = this.model.ci_lb;
    const ci_ub = this.model.ci_ub;
    const pval = this.model.pval;
    const alpha = options.alpha || 0.05;

    // Null value (0 for differences, 0 for log-ratios which is 1 on original scale)
    const nullValue = 0;

    // Check if CI excludes null (statistically significant)
    const ciExcludesNull = (ci_lb > nullValue || ci_ub < nullValue);
    const significant = pval < alpha;

    // Determine direction
    let direction;
    let directionMagnitude;
    let interpretation;

    // Default: lower is favorable for ratio measures (OR < 1 = protective)
    const favorLower = options.favorLower !== undefined ? options.favorLower :
      isRatio;

    if (!significant || !ciExcludesNull) {
      direction = 'not significant';
      directionMagnitude = 'none';
      interpretation = 'No statistically significant effect detected';
    } else if (estimate < nullValue) {
      // Effect is negative (on log scale for ratios, or directly for differences)
      direction = favorLower ? 'favorable' : 'unfavorable';
      directionMagnitude = isRatio ?
        (estimate < -0.69 ? 'large' : estimate < -0.22 ? 'moderate' : 'small') :  // Based on Cohen's d thresholds
        (Math.abs(estimate) > 0.8 ? 'large' : Math.abs(estimate) > 0.5 ? 'moderate' : 'small');
      interpretation = isRatio ?
        'Significant reduction in risk/odds (OR = ' + Math.exp(estimate).toFixed(3) + ')' :
        'Significant decrease (MD = ' + estimate.toFixed(3) + ')';
    } else {
      // Effect is positive
      direction = favorLower ? 'unfavorable' : 'favorable';
      directionMagnitude = isRatio ?
        (estimate > 0.69 ? 'large' : estimate > 0.22 ? 'moderate' : 'small') :
        (Math.abs(estimate) > 0.8 ? 'large' : Math.abs(estimate) > 0.5 ? 'moderate' : 'small');
      interpretation = isRatio ?
        'Significant increase in risk/odds (OR = ' + Math.exp(estimate).toFixed(3) + ')' :
        'Significant increase (MD = ' + estimate.toFixed(3) + ')';
    }

    return {
      direction: direction,
      magnitude: directionMagnitude,
      significant: significant,
      pval: pval,
      estimate: isRatio ? Math.exp(estimate) : estimate,
      ci: isRatio ? [Math.exp(ci_lb), Math.exp(ci_ub)] : [ci_lb, ci_ub],
      interpretation: interpretation,
      measure: measure
    };
  };

  /**
   * Compare effect from two models (e.g., fixed vs random)
   * @param {Object} model1 - First model result
   * @param {Object} model2 - Second model result
   * @returns {Object} Comparison metrics
   */
  MetaAnalysis.prototype.compareModels = function(model1, model2) {
    if (!model1 || !model2) {
      throw new Error('Two model results required for comparison');
    }

    return {
      estimateDiff: Math.abs(model1.estimate - model2.estimate),
      seDiff: Math.abs(model1.se - model2.se),
      model1: {
        estimate: model1.estimate,
        se: model1.se,
        method: model1.method
      },
      model2: {
        estimate: model2.estimate,
        se: model2.se,
        method: model2.method
      },
      agreement: Math.abs(model1.estimate - model2.estimate) < 0.1 * Math.max(Math.abs(model1.estimate), Math.abs(model2.estimate))
    };
  };

  // Static references to internal modules
  MetaAnalysis.EffectSizes = EffectSizes;
  MetaAnalysis.Tau2Estimators = Tau2Estimators;
  MetaAnalysis.Heterogeneity = Heterogeneity;
  MetaAnalysis.PublicationBias = PublicationBias;
  MetaAnalysis.Sensitivity = Sensitivity;
  MetaAnalysis.PooledEstimate = PooledEstimate;

  // Utility functions
  MetaAnalysis.utils = {
    sum: sum,
    mean: mean,
    variance: variance,
    normalCDF: normalCDF,
    normalQuantile: normalQuantile,
    tCDF: tCDF,
    tQuantile: tQuantile,
    chiSquareCDF: chiSquareCDF,
    chiSquareQuantile: chiSquareQuantile,
    chiSquarePDF: chiSquarePDF,
    lgamma: lgamma
  };

  // Export constants
  MetaAnalysis.TOLERANCE = TOLERANCE;
  MetaAnalysis.MAX_ITERATIONS = MAX_ITERATIONS;
  MetaAnalysis.CRITICAL_VALUES = CRITICAL_VALUES;

  // ============================================
  // EXPORT
  // ============================================

  // UMD export
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = MetaAnalysis;
  } else if (typeof define === 'function' && define.amd) {
    define([], function() { return MetaAnalysis; });
  } else {
    global.MetaAnalysis = MetaAnalysis;
  }

})(typeof window !== 'undefined' ? window : global);
