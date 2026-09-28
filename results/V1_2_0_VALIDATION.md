# Version 1.2.0 validation — 28 September 2026

This revision separates the report languages and improves arrowhead readability. Mathematical generators, coefficients, path calculations, invariant evaluation and stored certificates are unchanged from 1.1.0.

- All **16 tests passed**, including the original mathematical tests, the complete cube/tetra arc-label checks, and new language/arrowhead checks.
- Both editions contain **11 reports**. Their section structure, diagram geometry, numeric table cells and raw JSON agree. English report text contains no Chinese characters. Inline emphasis may move to follow translated sentence order.
- Report prose, headings, table labels, captions, controls and accessibility descriptions use paired English/Chinese text. Mathematical identifiers and machine-readable JSON field names remain shared. The bilingual root index only selects a report edition.
- Every generated arrowhead has a fixed display size (`markerUnits=userSpaceOnUse`, 11 by 9 SVG units); its shaft ends 7 SVG units before the endpoint center. Gray shafts are slightly thicker and muted arrows retain readable contrast.
- Inspected the bracket report's gray arrowheads in both editions and at mobile width. All **46 static browser checks** passed without document overflow or page-script errors.
- Checked **6 localhost routes in both languages**, including bracket, cube, tetrahedron, commutation and the reading guide. Switching languages preserves all mathematical case parameters.
- All offline links and anchors resolve across **23 pages** (22 reports plus the language chooser).
- Checked Chinese `show`, English `loop`, notebook code compilation and the 1.2.0 wheel build.
- Checked **32 mathematical source/data/certificate/provenance files** byte-for-byte against the 1.1.0 archive. No new full verification run is claimed for this presentation-only revision; the prior full run and its results remain documented in REVISION_2026-09-28.md.

Detailed records: v1.2.0-summary.json, v1.2.0-tests.txt, v1.2.0-browser-checks.json, v1.2.0-live-checks.json.

Open `reports/en/index.html` for English or `reports/zh/index.html` for Chinese. Run `python -m cocycle.examples` to regenerate both editions; use `--lang en` / `--lang zh` with `show`, `loop` or `serve`.
