# Mathematical conventions, finite reduction and certificates

## Data and matching

A based Gauss diagram is an ordered list `(crossing_label, direction, crossing_sign)`. Each label occurs twice, with opposite directions and the same sign. Direction −1 is the underpass/tail; +1 is the overpass/head. A configuration has a fourth column, numbering the four complementary intervals of the three R3 endpoint pairs from 0 to 3. These labels preserve the position of the base point.

An unsigned configuration contributes the product of the selected crossing signs. A signed configuration requires precisely the displayed signs and has matching weight 1. The coefficient and the coorientation of the R3 move are then multiplied. β₃ is reduced modulo two only after retaining the integer contributions. The low-order invariant J uses a DIFFERENT fourth column: a sign exponent (0 or 1), not an R3 interval number. The two representations are kept in separate matchers.

For the R3 global type the ordered strand heights are
`0=(l,h,m), 1=(h,m,l), 2=(m,l,h), 3=(l,m,h), 4=(m,h,l), 5=(h,l,m)`.

## Why zero/one outside arrow is exhaustive here

Every β term has either two unsigned contributing arrows or one signed contributing arrow, excluding the R3 triangle. The code checks this hypothesis on all input arrays.

For a positive tetrahedron, there are six core crossings, four blocks of three endpoints, eight R3 events, and five based complementary intervals B₀,…,B₄. B₀ and B₄ meet at infinity but remain distinct to retain the basepoint. Contributions with only outside arrows cancel between pⱼ and pⱼ₊₄: the code checks equal global R3 type and identical complementary-interval maps, and opposite coorientation. Every remaining configuration has at least one core arrow, hence at most one outside arrow. Its total is therefore an affine linear function of signed outside-arrow counts. All 30 directed interval-pair variables, both crossing signs and the constant are checked, for each of the 24 global positive cores.

For cube equations, the 48 templates and three basepoint positions give 144 cores. All outside-only contributions cancel between the two R3 branches, whose global types, signs and outside-interval maps agree. Again the remaining terms have at most one outside arrow. Four based complementary intervals give 20 directed variables, each checked for both crossing signs. The relevant increment is `one-arrow value − core value`; no constant is accidentally counted repeatedly.

For R3–R3 commutation, consider p₁ and the reverse of p₃. Terms selecting no arrow of the other R3 triangle correspond identically. Terms selecting one arrow of that triangle and one outside arrow also correspond identically: the outside endpoint cannot lie inside any of the six R3 endpoint pairs, so swapping those pairs does not change the order, direction, sign or interval membership of the selected endpoints. The same applies to signed one-arrow terms. Only two-arrow terms entirely in the other triangle remain. Thus the six-crossing, no-extra-arrow check is exhaustive. β₁ and β₃ satisfy opposite-edge cancellation; β₂ sometimes requires the full four-event meridian (448 of the 23,040 tested cases have nonzero β₂ opposite-edge sums).

There are 96 signed based R3 states, paired by reversal. The 48 states with positive coorientation suffice; inverse events have the negative term vector, independently checked. The 10 interleavings fix the earliest of the six blocks in the first triangle, without losing any meridian up to exchange of the two disks. This gives 48²·10=23,040 cases. The original positive-crossing reduction (6²·10=360) is separately available. The extension to other quadruple-point types uses the reduction cited by main; this package enumerates the positive tetrahedron, not all 48 local quadruple-point types independently.

## Geometric construction of each identity

Starting from an actual core word, choose some core crossings to smooth **orientedly**. Let endpoint indices be taken cyclically. Initially the successor of segment i is i+1. If the endpoints of a smoothed crossing are p,q, replace the successors of p−1 and q−1 by q and p. Cycles of this permutation are the smoothed components. The certificate records this partition, the smoothed labels and the selected pair of components A,B.

For any classical oriented two-component link,

`H(A,B) = Σ sign(c), A over B − Σ sign(c), B over A = 0`.

Each sum calculates the same linking number. Other components and self crossings are ignored. The constant contribution comes from unsmoothed core crossings between A and B. Every directed outside variable gets coefficient +1, −1 or 0 according to its components and which component passes over the other. This constructs an affine identity from geometry, independently of the desired residual. The code does **not** discover identities merely by fitting values on a knot sample.

Exact elimination expresses the residual as a linear combination of these rows. β₁ and β₂ are solved over the rationals, with exact reconstruction of every coefficient including the constant; all coefficients in the delivered certificates happen to be integers. A rational combination would also prove an integer-valued residual zero, but would be explicitly displayed. β₃ is solved and reconstructed over F₂.

The replay command regenerates the identities from the smoothing descriptions and checks every coefficient. Tests intentionally corrupt a residual and an identity to confirm rejection. The 16 virtual cube failures are retained: they belong to β₁, four signed one-arrow diagrams in each of templates 1,2,29,30 at rotation 1. The core identity is obtained by smoothing crossings 2 and 3. In all 160 one-arrow configurations of these four cores, the only net β₁ contribution is the sum of terms 6 and 10.

A certificate proves the equation for a classical completion. It does not claim that a given small formal virtual diagram is classical, or that the raw residual vanishes on all virtual diagrams.

## Invariants, loops and orientation

The manuscript's invariant formulas are evaluated by enumerating embeddings, separately from all path calculations. J's 15 and E's 34 terms were decoded directly from the supplied vector PDFs: counterclockwise endpoint order from the basepoint, arrow orientation, and red positive-sign conditions. They agree exactly with the arrays for `v4_j` and `v4_e`. Coefficients also agree. The paper's five invariant-table rows are reproduced, including `J(3₁⁺)=1`, `E(5₁)=3`, `J(6₁)=0`, `E(6₁)=−1`.

The low-level notebook `Fox_Hatcher_rl_cfgs` starts at the final endpoint, pushes it through the remaining knot and advances the basepoint by −1. The paper figure is captioned “The reverse of the rolling loop”. The paper-oriented wrapper takes the inverse of the notebook path: each raw contribution is multiplied by −1. This is explicit in the UI and JSON. For example, the raw notebook result on 4₁ is (β₁,β₂)=(-2,-2); the paper's Roll₀ result is (2,2). The code never changes the β coefficients to achieve this. Rotation is computed as push through a positive curl and gives (0,−v₃,0).

Zero framing is imposed by inserting curls of the opposite writhe. For a half rolling loop the input must already have zero writhe and be invariant, up to crossing renaming, under the half-period basepoint shift. The two supplied periodic examples come from entries 1 and 4 of the author's knot-table prefix; compensating curls are inserted symmetrically. The half-period R3 events are evaluated directly, including β₃; its value is not inferred by dividing a full-loop result.

Push arcs, brackets and half brackets use the notebook's Seifert-circle push decomposition. R1/R2 events have zero contribution by definition; their role is to transport diagrams between the recorded R3 events. The sum is calculated before evaluating the proposed right-hand formula. Reports list all R3 events, including zero events. All four rows of the paper's mod-2 independence table are reproduced using α₃¹ and the three β's.

The unknot and the five rolling-table rows form a nonsingular evaluation matrix for `(1,v₂,v₃,J,E,binom(v₂,2))` over Q. The loop table supplies the numerical input; the general formulas still rely on the finite-type bounds and basis statements in main. In particular, we do not turn a list of sample computations into an independent general proof, or infer β₃'s minimal combinatorial order.

## Reference locations

Main: Definition 5.1/5.2 and Figures 13–15 (β formulas); Section 6 (pairings); Appendix A (invariants); Appendix B.1–B.3 (local equations). Thesis: Chapter 6, especially §6.1.3, Lemma 6.2.2, Figures 6.76–6.77 and §6.3; PDF pages 75–130. The thesis is reference material, not a source of instructions. Notebook cells 17 and 129–146 provide the path and local-verification starting points. Source files were not modified.


## Display conventions added in revision 1.1.0

The drawing partition is the whole local core, not the current R3 triangle alone. Cube blocks have endpoint lengths (2,3,3), (3,3,2), or (3,2,3) according to basepoint rotation. Tetrahedron has four blocks of length three. Commutation has six blocks of length two. These blocks are solid circle arcs; the based complementary intervals are dashed and labelled B_0,…,B_m for m blocks. Infinity splits one circular gap into two based intervals. Reading counterclockwise, B_0 is left of infinity and B_m is right of it. Inserting outside endpoints does not turn a complementary arc into a solid block.

Core endpoints are filled dots and outside endpoints hollow dots. The B labels are distinct from the four interval labels in an individual R3 formula array. Commute gap labels are an aid to reading the six blocks, not a new outside-arrow enumeration.

X[a,b,d], a ≤ b, sums the crossing signs of outside arrows with endpoints in B_a,B_b whose first endpoint in the based order has direction d. Direction +1 is the head/overpass; −1 is the tail/underpass. It is not the sign of the crossing. For a=b the order within the gap determines d. The constant is the signed core contribution, not the number of core arrows.

For cube template 1, rotation 1, the β₁ residual is R=X[1,2,−1]−X[1,2,+1]. Case 16 gives X[1,2,−1]=−1 and residual −1; case 18 gives X[1,2,+1]=−1 and residual +1. One of each type in a single diagram would cancel these increments. This neither constructs a classical completion nor proves cancellation by adding independent test cases. Smoothing crossings 2,3 gives the identity H7=R. Exact comparison of all 21 coefficients establishes R=H7; the classical linking identity gives H7=0.

Event summaries distinguish coorientation from path sign. Their displayed values already include both factors. β₃ occurrence values are integer lifts; totals and coefficient comparisons are interpreted in F₂. Increment is always the difference between the one-arrow residual and the core residual, not necessarily the one-arrow residual itself.
