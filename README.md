# Finite-type-combinatorial-1-cocycles
This project is aimed for the calculation or the proof of the results concerning finite type combinatorial 1-cocycles.

## Order-four 1-cocycles: exact verification

Computational companion to **Butian Zhang, Pairings of combinatorial 1-cocycles with loops in knot spaces** (manuscript dated 21 September 2026; sources supplied 24 September 2026).

[中文使用说明](docs/使用说明.md) · [Mathematical conventions and coverage](docs/MATHEMATICS.md) · [Results](results/REPORT.md) · [Example reports](reports/index.html)

This is a mathematical case inspector. Every requested case can display every R3 event, every occurrence of a formula term, its selected arrows, coefficient, crossing signs, coorientation and contribution. Tetrahedron and cube reports also display exact linking-number identities, the crossings to smooth, their components and the coefficient-by-coefficient certificate. Opposite contributions within a term are retained separately.

## Choose a report language (revision 1.2.0)

- [English reports](reports/en/index.html)
- [中文版报告](reports/zh/index.html)

Every report has a complete English edition and a corresponding Chinese edition. Explanations, captions, tables and controls are translated; mathematical identifiers and raw JSON keys remain shared. Arrowheads have a fixed size and stop before endpoint dots, including gray arrows in loop reports.

Use `--lang en` or `--lang zh` with `show`, `loop` or `serve`. The notebook has a language selector.

## Reading the diagrams

Open [the offline reading guide](reports/en/reading-guide.html), then the [worked cube cases 16 and 18](reports/en/cube-1-1-beta1.html#case-16-18).

- **Solid circle arcs:** fixed core endpoint blocks — cube 3, tetrahedron 4, R3–R3 commute 6.
- **Dashed circle arcs:** complementary intervals B_0,… where outside endpoints may be inserted. Reading counterclockwise from the top basepoint, B_0 is immediately to its left and the last B is immediately to its right.
- **X[a,b,d]:** sum of crossing signs of outside arrows in gaps B_a,B_b. The endpoint direction d is distinct from the crossing sign being summed.
- **Residual increment:** R(core + one outside arrow) − R(core), shown only for the selected cocycle. A nonzero formal one-arrow residual is not hidden.
- **Certificate:** the program compares every coefficient of R with a geometrically constructed combination of linking identities. It does not claim cancellation within each outside-arrow case or a one-to-one pairing of test cases.

See [language-edition validation](results/V1_2_0_VALIDATION.md), [revision notes](CHANGELOG.md) and [new validation](results/REVISION_2026-09-28.md).

## Start here

Python 3.10 or later is required. From this directory:

```sh
python -m pip install .
python -m cocycle serve
```

Open **http://127.0.0.1:8765** in a browser. Select a tetrahedron, cube, commutation case or a loop from the form. No Jupyter or native C extension is needed. The bundled HTML examples can be opened directly without installing anything.

For an interactive mathematical notebook, open `Explore.ipynb` after installing `python -m pip install '.[notebook]'`.

## Reproduce the calculations

```sh
python -m unittest discover -s tests -v
python -m cocycle verify --output new-results
python -m cocycle replay --certificates new-results/local-certificates.json
python -m cocycle.examples  # rebuild both language editions from bundled records
```

The full run covers:

- **Tetrahedron:** 24 positive global cores, with 0 and 1 outside arrow; 1,440 signed one-arrow diagrams. Every residual is certified using geometrically constructed linking identities.
- **Cube:** 144 global cores (48 templates × 3 basepoint positions), with 0 and 1 outside arrow; 5,760 signed one-arrow diagrams. The 16 nonzero virtual β₁ residuals are retained and certified, including the term 6 / term 10 identity.
- **R3–R3:** 23,040 signed cases (48 coorientation-positive global R3 types for each disk × 10 interleavings). Reversing an event handles the other coorientations. β₁ and β₃ cancel on opposite edges; β₂ is checked on the full meridian.
- **Loops:** 136 comparisons of actual R3 sums against the paper's formulas, including the paper's rolling, push and mod-2 independence tables and 50 additional knot-table diagrams for rotation/rolling.

The original C matcher is optional. When it is installed, the local verification cross-checks all terms against it. The supplied run records **115,680 successful event comparisons**. The portable Python implementation supplies all actual computations and all detailed occurrences.

## Request a particular case

```sh
python -m cocycle show tetra --case 1 --beta beta2 --output tetra1-beta2.html
python -m cocycle show cube --case 1,1 --beta beta1 --output cube1-beta1.html
python -m cocycle show commute --case 0,8,7 --beta beta2 --output commute.html
python -m cocycle loop rolling 4_1 --output rolling-4_1.html
python -m cocycle loop half-rolling 6_1_periodic --output half-rolling-6_1.html
python -m cocycle loop bracket 4_1 --other 6_1 --output bracket.html
```

Without `--outside`, tetrahedron/cube reports include **all** signed outside placements and every contribution. To inspect one formal virtual diagram, add `--outside N` (tetrahedron: 0–59; cube: 0–39). These formal diagrams need not be classically realizable: a nonzero raw value is not suppressed.

Tetrahedron `s` is 1-based, as in the reference notebook (`template=48*(s-1)`). Cube templates, rotations, outside indices and commutation triples are 0-based. Commutation types are specified in `docs/case-index.json`; `--positive` selects the original six-type/360-case notebook convention.

## Important conventions

- **`v41 = v4_j = J15`, `v42 = v4_e = E34`.** All 49 arrow diagrams agree with the supplied manuscript PDFs, including partial sign conditions. The notebook's older `v_4_1_primitive` and `v_4_2_primitive` are different normalizations and must not be substituted.
- **Rolling direction:** the notebook's Fox–Hatcher routine traverses the inverse of the paper's rolling convention. We show its raw sum and apply an explicit path sign −1. See the caption “The reverse of the rolling loop” in main and the detailed orientation discussion.
- **β₃ is over F₂.** Its integer lift is shown for auditing, then reduced modulo two.
- **Half rolling requires a genuinely periodic, zero-framed input diagram.** It is not computed by dividing a mod-2 number by two.
- Loop examples and interpolation data are computational inputs to the paper's mathematical arguments; a finite sample is not advertised as a new proof of the general formula.
- The proof that **β₃ has combinatorial order exactly four is explicitly excluded**. No lower-order cocycle search is performed here.

## Files and provenance

`cocycle/data/` contains the formula arrays and explicit diagrams (no pickle, no large external database). `cocycle/geometry/` adapts the author's supplied generators; `paths.py` adapts the push/rolling paths from `For_these.ipynb`. `local.py` preserves and audits the notebook's transparent local verification. `SOURCE_HASHES.json` identifies the exact source files; `RELEASE_SHA256.json` identifies the delivered files.

`results/local-certificates.json` contains all residual vectors, geometrically generated identities, exact certificate coefficients and commutation event values. `results/loop-checks.json` contains the independently calculated invariant values, actual loop sums and predicted formulas. Reports are self-contained and require no remote scripts.

GitHub Actions runs the tests and the complete verification. Publishing to GitHub is separate from preparing this repository; no remote repository is created by this package.
