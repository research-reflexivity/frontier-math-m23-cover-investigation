# Source notes

The three principal canonical JSON inputs are:

- `data/Fint_coefficients_Z.json`
- `data/optimal_23_4_Z.json`
- `data/optimal_degree4_pencil.json`

Their SHA-256 fingerprints in this repository are:

```text
27dd1fd0de4f1c350f5a07ead6d5b747d7a8f34e4146ec42fa953f1536a65102  data/Fint_coefficients_Z.json
cfb4185b2800744d524d87b3f08f5459cd38d42afb3787c11fa1e2aa3f25ee84  data/optimal_23_4_Z.json
fdf571e4f47effc527e1726bdd99472be2c5e3390e45133fb3d7a3aab592619d  data/optimal_degree4_pencil.json
```

Their mathematical provenance and derivation are recorded in the present
paper.  The integral source equation comes from the model of Huang, Jackson,
Lee, Poonen, Pries, and Zhang cited as `HJLPPZ` in `paper/main.tex`.  They
prove that a degree-four function exists and is optimal but do not compute
it; the explicit equation and pencil realizing that function, together with
the arithmetic of the fibre above `T=0`, are derived and established here.

The repository includes the complete certificate suites for construction of
the minimal model, its arithmetic specializations, the Fano and affine
structures in the fibre above `T=0`, and their exact transport under the
change of generator.  The historical filename
`verification/verify_special_fibre_bridge.py` is retained as a stable
certificate interface; it reconstructs this comparison from the three JSON
inputs above.
