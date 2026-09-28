# Verification results — 24 September 2026

This is the original computation record. For the 28 September presentation revision and new validation, see [REVISION_2026-09-28.md](REVISION_2026-09-28.md).

All requested computational families passed, with the distinctions below. The β₃ combinatorial-order minimality computation was excluded.

| Calculation | Exhaustive cases / comparisons | Result |
|---|---:|---|
| Positive tetrahedron core | 24 | All 72 cocycle certificates reconstruct every coefficient |
| Tetrahedron, one signed outside arrow | 1,440 | Raw virtual nonzero counts: β₁ 122, β₂ 372, β₃ 56; all covered by the classical certificates |
| Cube core | 144 | All 432 cocycle certificates reconstruct every coefficient |
| Cube, one signed outside arrow | 5,760 | 16 raw virtual β₁ residuals; all resolved by the displayed smoothing identity |
| Cube β₁ terms 6 / 10 | 160 | All other terms have zero net sum; residual is exactly the stated pair |
| Signed R3–R3 commutations | 23,040 | Every four-event sum is zero; β₁ and β₃ opposite edges cancel |
| Python / original native matcher | 115,680 R3 events | Every term agrees |
| Loop/formula comparisons | 136 | All specified comparisons agree |
| Manuscript vector diagrams | 55 | All endpoint words, directions and sign conditions agree |
| Mathematical regression checks | 9 | All pass, including corrupted-certificate rejection |

A separate full run with the native matcher disabled reproduced **identical** local certificates, commutation event values and all 136 loop results. The package does not require the author's C extension or absolute paths. A distributable Python wheel was built successfully.

## Rolling table (paper orientation)

| Knot | β₁(Roll₀) | β₂(Roll₀) | β₃(Roll₀) | v₂ | v₃ | v₄,₁=J | v₄,₂=E |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3₁⁺ | 0 | −6 | 0 | 1 | 1 | 1 | 0 |
| 4₁ | 2 | 2 | 0 | −1 | 0 | 0 | 0 |
| 5₁ | 0 | −50 | 0 | 3 | 5 | 5 | 3 |
| 5₂ | 0 | −26 | 0 | 2 | 3 | 3 | 1 |
| 6₁ | 8 | 10 | 0 | −2 | −1 | 0 | −1 |

The raw notebook Fox–Hatcher sums have the opposite sign over Z. Every report records the raw sum and the explicit inverse-path sign. No β coefficient is changed.

The five rows, together with the unknot, give determinant **−1** for the evaluation matrix of `(1,v₂,v₃,J,E,binom(v₂,2))`.

## Independence table, modulo two

| Loop | α₃¹ | β₁ | β₂ | β₃ |
|---|---:|---:|---:|---:|
| Rot(4₁) | 1 | 0 | 0 | 0 |
| ½Roll(4₁) | 0 | 1 | 1 | 0 |
| Rot(6₁) | 0 | 0 | 1 | 0 |
| ½Roll(6₁) | 1 | 0 | 1 | 1 |

Half rolling is evaluated on explicitly periodic zero-framed diagrams, not by halving a mod-two value.

## Scope and source observations

- The exact equality `v₄,₁=v4_j` and `v₄,₂=v4_e` was checked diagram by diagram, not just inferred from sample values. The older functions with “primitive” in their names have different normalizations.
- Main Appendix B.2 says “4 arrows from the quadruple point”. A quadruple point has **four strands and six crossing arrows**, as the thesis and generator use. The arXiv package preserves the supplied manuscript; this textual slip was not silently edited.
- Cube identities are classical linking identities. Nonzero formal virtual residuals are deliberately preserved and displayed.
- The full signed commutation check extends the 360 positive cases of the reference notebook. β₂ has 448 cases with nonzero opposite-edge sums but zero full-meridian sum.
- The 136 loop comparisons comprise 24 named-knot rotation/rolling/half-bracket cases, 4 brackets, 4 push-table cases, 4 rotation/half-rolling cases on periodic diagrams, and 100 rotation/rolling calculations on the first 50 diagrams of the author's table.
- General loop identities still use the finite-type bounds and basis arguments in main. The program supplies their computational input and additional checks, rather than claiming a sample-only proof.
- No code to classify all lower-order mod-two cocycles or prove β₃'s minimal combinatorial order is included.

The local case explorer and its cube/half-rolling endpoints were tested, including opening a nonzero cube occurrence and its identity in the browser. The HTML examples work offline. Source hashes, explicit Gauss matrices and exact coefficients are included.
