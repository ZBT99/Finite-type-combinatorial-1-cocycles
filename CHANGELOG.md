# Revision 1.2.0 — 28 September 2026

- Separate every offline report into complete `reports/en/` and `reports/zh/` editions. Corresponding paragraphs, captions, controls and tables use paired text; mathematics and JSON identifiers are shared.
- Add report-language selection to the command line, localhost explorer and notebook.
- Give all arrowheads a fixed display size independent of shaft thickness, shorten shafts before endpoint dots, and improve gray-arrow contrast.
- Rebuild the formula catalogue from the same formula arrays and original manuscript images in both languages.
- Check language parity, English text purity, identical diagram geometry and raw data, and arrowhead clearance.

See `results/V1_2_0_VALIDATION.md`.

# Revision 1.1.0 — 28 September 2026

This revision makes the distinction between raw virtual residuals and classical identity certificates explicit. Formula arrays, path calculations, local generators and the original mathematical certificates are unchanged.

- Draw whole-core endpoint blocks as solid arcs: three for cube, four for tetrahedron, six for R3–R3 commutation. Draw complementary based intervals as dashed arcs and label them B_0,… from the basepoint. Outside endpoints remain on dashed arcs.
- Distinguish crossing labels, endpoint row numbers, gap labels, crossing signs, endpoint direction, R3 coorientation and path orientation.
- Replace the erroneous Unicode subscript beta in B_b. Explain X as a sum of crossing signs, and the constant as a signed core contribution rather than a count of core arrows.
- Show a named residual increment for the selected cocycle, plus its core and full one-arrow residual. Preserve the integer lift for β₃ while displaying its residue modulo two.
- Explain cube (1,1) cases 16 and 18 with the actual term-6 and term-10 contributions. Conditional cancellation in one complete diagram is not presented as a pairing of independent test inputs or as a constructed classical completion.
- Separate the computed residual, geometric linking identities and exact coefficient comparison. Display every coefficient, including the constant, rather than only a final zero.
- Add an offline reading guide, improve the Chinese guide, and rebuild examples through `python -m cocycle.examples`. The live browser and notebook share the revised renderer.
- Add checks of all 5,760 cube and 1,440 tetra outside-arrow labels against their residual polynomials, plus commutation block and case-16/18 smoothing checks.

See `results/REVISION_2026-09-28.md` for this revision's validation. The 24 September reports remain historical computation records; their native-matcher audit is not represented as a new native run.
