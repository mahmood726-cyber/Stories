# Stories

Interactive narrative evidence-communication pieces built around historical
clinical-trial reversals.

## Contents

- **`CAST_WhenCertaintyKills.html`** — a single-file interactive narrative on the
  Cardiac Arrhythmia Suppression Trial (CAST), where anti-arrhythmic drugs that
  suppressed premature ventricular contractions increased mortality. Readers
  predict outcomes at each decision point before the results are revealed.
- **`Evidence_Reversal_40/`** — a multi-story interactive app (40 evidence
  reversals) plus a self-contained meta-analysis engine (`meta-analysis.js`):
  effect sizes (log OR / RR / RD, Fisher-z), tau-squared estimators
  (DL, REML, ML, PM, HE), heterogeneity (I-squared, Q, Q-profile CI),
  fixed- and random-effects pooling, prediction intervals, and publication-bias
  tests (Egger). Plotly is vendored locally (`plotly-2.27.0.min.js`).
- **`e156-submission/`** — the E156 micro-paper submission artifacts.
- **`E156-PROTOCOL.md`** — the E156 protocol and abstract.

## Engine tests

The meta-analysis engine has a Node smoke test that pins hand-verifiable
numerical anchors:

```
node Evidence_Reversal_40/smoke.test.js
```

## Notes

The HTML apps are self-contained and run offline; the only external reference is
an optional Google Fonts stylesheet that degrades to system fonts if unavailable.
