# E156-PROTOCOL — CardioSynthesizer (DAPT Reanalysis Engine)

- **Project:** Equipose (GitHub repo `Equipose`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `Equipose.html` dump,
  "CardioSynthesizer V5.5: Reanalysis Engine")
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: removed the Google Fonts CDN `<link>` (replaced with
  an offline comment); the app now loads no external resource.
- Extracted the statistical core (DerSimonian–Laird pooling, Egger's test,
  odds-ratio→risk-difference transform, histogram builder) into a pure
  `engine.js` (single source of truth). The inline worker duplicates were
  deleted; the worker blob is now built by serialising `engine.js`'s functions
  into it, so the worker thread runs identical, test-verified code offline.
- Added `tests.js` (30 assertions, all passing) with hand-computed expectations.
- Verified pooling math correct and left unchanged; confirmed the k=1
  single-study path returns the study (τ²/I² clamp to 0) without `NaN`.
- Added Pages scaffold (`.nojekyll`, `.gitignore`, README); renamed
  `Equipose.html` → `index.html`.

## Body (E156 draft — CURRENT BODY)

When does extended versus short dual antiplatelet therapy deliver net clinical
benefit once anatomical complexity and bleeding risk are weighed together? This
offline dashboard reanalyses eleven DAPT trials, pooling log odds ratios for
ischemic benefit and bleeding harm. It fits a DerSimonian–Laird random-effects
model, reporting τ², I², a pooled estimate and a 95% prediction interval. It
then runs a forty-thousand-draw Monte-Carlo net-clinical-benefit simulation with
an explicit complexity interaction and an ischemia–bleeding correlation. Across
the bundled trials the net benefit stays narrow and the probability of benefit
hovers near the equipoise band, so no strategy dominates once a realistic
bleeding penalty applies. A revival audit extracted the statistical core into a
single tested source of truth and confirmed the pooling math and the k=1 no-NaN
behaviour against thirty hand-derived assertions. The honest read is
decision-dependent equipoise rather than a universal recommendation, and the
tool exists to make that trade-off transparent.

SUBMITTED: [ ]
