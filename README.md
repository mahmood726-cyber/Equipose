# CardioSynthesizer — DAPT Reanalysis Engine (Equipose)

A single-file, **fully offline** dashboard that reanalyses dual antiplatelet
therapy (DAPT) trial evidence. It pools log odds ratios for ischemic benefit and
bleeding harm with a DerSimonian–Laird random-effects model, runs a Monte-Carlo
net-clinical-benefit (NCB) simulation, and renders a tipping-point curve, forest
plot, risk-benefit scatter, time-horizon twin, and NCB histogram. Effect
modification by anatomical complexity is applied through an explicit interaction
term (no hidden scalars).

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN.

## Layout

```
index.html   single-file UI (loads engine.js; builds the worker from it)
engine.js    pure statistical core — runs in Node and the browser
tests.js     Node test harness, 30 assertions
LICENSE      Apache-2.0
```

## Statistical core (`engine.js`)

| Function | What it does |
|---|---|
| `runREML(data, kEff, kSe)` | DerSimonian–Laird random-effects pooling of log odds ratios: τ² via method-of-moments `(Q−(k−1))/C`, I², pooled SE, and a 95% prediction interval; returns per-study forest CIs |
| `eggerTest(data, kEff, kSe)` | Egger small-study-effects regression of effect on precision (intercept + coarse significance flag); `null` for k<3 |
| `probFromLogOR(base, logOR)` | converts a log odds ratio to an absolute risk difference at a given baseline risk |
| `buildHistogram(arr, bins)` | bins the central 5th–95th-percentile band of the Monte-Carlo NCB distribution |

The page loads `engine.js` as the page-side source of truth and serialises the
same four functions (`.toString()`) into the Web Worker blob, so the worker
thread runs the identical, test-verified code with no second copy — and with no
network fetch (works from `file://`).

## Fixes applied during revival (2026-06-05)

- **Offline:** removed the Google Fonts `<link>`; the app now loads no external
  resource (system fonts fall back).
- **Single source of truth:** extracted the pure stat functions into `engine.js`
  and deleted the inline duplicates from the worker block; the worker is now
  built by concatenating `engine.js`'s functions into its blob.
- **Tests:** added `tests.js` (30 assertions, all passing), with two hand-worked
  DL pooling cases, a k=1 no-NaN check, a two-identical-study τ²=0/I²=0 check,
  and edge guards.
- Renamed `Equipose.html` → `index.html`; added Pages scaffold (`.nojekyll`,
  `.gitignore`).

The pooling math was independently re-derived and verified correct, so it was
left unchanged. The DerSimonian–Laird estimator already clamps τ² and I² to 0,
so the k=1 single-study case returns the study itself without producing `NaN`.

## Tests

```
node tests.js
# 30 passed, 0 failed
```

Checks include a hand-computed homogeneous 2-study case
(μ≈−0.18951, SE≈0.038411, τ²=0, prediction interval [−0.26479, −0.11422]),
a heterogeneous 2-study case (τ²=0.17, I²≈94.44%), a k=1 passthrough that must
not be `NaN`, a two-identical-study case (τ²=0, I²=0, SE=√(1/200)), the
odds-ratio→risk-difference transform, an Egger constant-effect intercept, and
histogram edge guards.

## Caveats

DerSimonian–Laird is known to under-estimate τ² for small *k* (REML or
Paule–Mandel are preferred for k<10); this dashboard preserves the original
method for continuity and reports τ² and I² alongside every estimate. The
Monte-Carlo NCB, time-horizon curves, and complexity interaction are explicit
modelling assumptions, not fitted quantities — treat outputs as
hypothesis-generating, not a clinical decision rule. Apache-2.0 licensed.
