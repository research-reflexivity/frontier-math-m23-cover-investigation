# Magma calculator rerun — 7 September 2026

**All 16 public inputs have successful Magma 2.29-10 executions.**
They were run sequentially on the University of Sydney's public calculator,
using the same XML endpoint as its web form. The limits were 60 seconds and
50,000 input bytes per request. No mathematical assertions were removed or
weakened. The exact input files match those published at public GitHub commit
`158c9afc6c06025197ecde74b8e8e0bbe8c96c92` before this rerun.

## The newly executed harmonic-model certificate

`verification/verify_harmonic_model_generated.m` passed unchanged in
**4.009 seconds** (39,831 input bytes). The input SHA-256 is
`ba42c7ecd378db402be2c0a92a3cf1d22a38bac14b40fd9ef63693aa0f1b88d0`.
Its recorded calculator seed is `3859956797`.

The successful assertions cover:

- smoothness in all four projective charts;
- exact binary quadric/cubic divisibility and nonzero discriminant for the
  two-point incidence on the osculating line;
- the canonical vanishing sequences at the two marks and the supplied
  quintic sections' order-23 base divisors;
- a one-dimensional multiplier kernel, its invertible witness, all six
  cross-products, and the identity with the supplied degree-23 function;
- a radical eight-dimensional residual critical algebra with a common
  nonzero critical value, giving the stated three-point passport together
  with the marked fibres and Riemann–Hurwitz.

It does **not** perform the analytic Hensel audit, the pi^33 comparison,
explicit rational descent, a new Arb continuation, or a mechanical review
of the geometric proof. The group identification still uses the separate
exact map-isomorphism and existing branch-cycle certificate. The finite
orientation/Bockstein checks below do not supply the unresolved relative
Fano–affine correspondence.

## Results

Times are calculator-reported seconds, not network wall time. See the
[machine-readable index](verification/calculator_2026_09_07_magma_summary.json)
for input hashes, versions, seeds and the exact request/response records.

| Input | Result | Seconds |
| --- | --- | ---: |
| `verification/verify_harmonic_model_generated.m` | PASS | 4.009 |
| `notes/compute_q0_tail_function_field_galois.m` | PASS after explicit-seed retry | 9.490 |
| `notes/compute_rational_e8_tail_function_field_galois.m` | PASS | 12.000 |
| `verification/verify_optimal_23_4.m` | PASS | 13.349 |
| `verification/certify_gauss_prolongation_obstruction.m` | PASS | 0.050 |
| `verification/certify_geometric_reflection_obstruction.m` | PASS | 0.070 |
| `verification/verify_canonical_quadric.m` | PASS | 36.200 |
| `verification/verify_hurwitz_degree23_branch.m` | PASS | 0.130 |
| `verification/verify_hurwitz_degree23_geometry.m` | PASS | 0.110 |
| `verification/verify_hurwitz_galois_closure.m` | PASS | 0.100 |
| `verification/verify_hurwitz_local_23.m` | PASS | 0.030 |
| `notes/certify_fano_affine_odd_fixed_point_lemma.m` | PASS | 0.110 |
| `notes/certify_pinched_tag_finite_identities.m` | PASS | 0.360 |
| `notes/certify_wild_parameter_orientation.m` | PASS | 0.020 |
| `notes/audit_pointed_relative_bockstein.m` | PASS | 1.159 |
| `notes/audit_log_quadratic_orientation_line.m` | PASS | 0.120 |

## Attempts retained, not hidden

There are **18 request/response records** for the 16 inputs.

1. The conductor-seven tail's first randomly initialized run exceeded the
   calculator's 60-second limit. This was an incomplete computation, not an
   assertion failure. Prepending only `SetSeed(3275039696);`, the seed in the
   earlier successful public record, gave a complete result in 9.490 seconds.
   The record preserves both the original source and submitted wrapper hashes.
   Its server-start seed differs from the effective seed set in the program.
2. The initial transcript reader failed to recognize one long geometry
   completion message because Magma wrapped it over two lines. The raw response
   already showed success with no error. The parser now normalizes whitespace
   when matching fixed messages, and the geometry script was rerun successfully.
   Both responses and the initial reader classification are retained.

The parser rejects a result containing an assertion/runtime/syntax error or
server warning even if a later line says PASS. Seven local regression tests
cover these cases, including timeouts, missing completion and line wrapping.
The archived XML is an execution/provenance record, not a digitally signed
formal proof certificate.

## Reproduction

Offline validation of the stored evidence:

```sh
make verify-magma-calculator-records
```

This checks deterministic generation of the harmonic input, current input
hashes and byte counts, exact submitted wrappers, raw response hashes,
completion markers and the index. It is now included in `verify-magma-record`
and hence in the default `verify-all`, but it does not contact Sydney.

An explicit fresh run that transmits the listed public scripts to Sydney:

```sh
make run-magma-calculator
```

To use the recorded successful conductor-seven seed:

```sh
python3 scripts/run_magma_calculator.py \
  --only notes/compute_q0_tail_function_field_galois.m --seed 3275039696
```

New runs retain new records rather than overwriting old ones. After reviewing
them, regenerate the index with
`python3 verification/verify_magma_calculator_runs.py --write-summary`.
The old Magma 2.29-9 summary files remain intact and continue to describe their
original executions.
