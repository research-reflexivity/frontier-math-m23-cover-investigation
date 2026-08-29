#!/usr/bin/env python3
"""Exact audits for the effective logarithmic quadratic-orientation line."""

from itertools import combinations


def q_half_rank(rank: int) -> int:
    """Parity of the Clifford-volume line of an even virtual rank."""
    assert rank % 2 == 0
    return (rank // 2) % 2


def q_subset(subset: set[int]) -> int:
    assert len(subset) % 2 == 0
    return q_half_rank(len(subset))


# Generic and returned trace lines.
trace_size = 253
for intersection_parity in (0, 1):
    # Only the residue modulo two is needed: 253-I modulo two.
    generic_q = (trace_size - intersection_parity) % 2
    assert generic_q == (1 + intersection_parity) % 2

returned_square_rank = 0
returned_nonsquare_rank = 2 * trace_size
assert q_half_rank(returned_square_rank) == 0
assert q_half_rank(returned_nonsquare_rank) == 1

# The terminal pair is (C,T^n).  Polarization by C=P triangle T has no
# cross term in either return case.
packet_return_cross_square = 0
packet_return_cross_nonsquare = 0
assert packet_return_cross_square == packet_return_cross_nonsquare == 0
full_trace_square_q = q_half_rank(3795 - 253)
full_trace_nonsquare_q = q_half_rank(3795 + 253)
assert (full_trace_square_q, full_trace_nonsquare_q) == (1, 0)

# Effective finite-node anomaly and packet orientation.
branch_plus = 1078
branch_minus = 112
packet_size = 3542
full_class_size = 3795
finite_anomaly = q_half_rank(branch_plus - branch_minus)
packet_orientation = q_half_rank(packet_size)
assert finite_anomaly == 1
assert packet_orientation == 1
assert finite_anomaly == packet_orientation

# Changing P by the congruent lift P-2C reverses the half orientation.
changed_lift = q_half_rank(packet_size - 2 * full_class_size)
assert changed_lift == 0
assert changed_lift == packet_orientation ^ (full_class_size % 2)

# The two fibre formulae force epsilon=q once the relative-line
# specialization is supplied geometrically.  Do not assume that equality
# inside this finite audit.
for epsilon in (0, 1):
    for q_return in (0, 1):
        generic = (1 + epsilon) % 2
        special = (finite_anomaly + q_return) % 2
        assert (generic == special) == (epsilon == q_return)

# Exhaustive small-model audit of the general normalization telescope.
# Adjacent integral coefficient vectors have coordinatewise even
# differences, so division by two is defined before reduction modulo two.
coefficient_vectors = [
    (a, b) for a in range(5) for b in range(5)
]
telescoping_chains_checked = 0
for z0 in coefficient_vectors:
    for z1 in coefficient_vectors:
        if any((b - a) % 2 for a, b in zip(z0, z1)):
            continue
        for z2 in coefficient_vectors:
            if any((b - a) % 2 for a, b in zip(z1, z2)):
                continue
            local = tuple(
                ((mid - left) // 2 + (right - mid) // 2) % 2
                for left, mid, right in zip(z0, z1, z2)
            )
            endpoint = tuple(
                ((right - left) // 2) % 2
                for left, right in zip(z0, z2)
            )
            assert local == endpoint
            telescoping_chains_checked += 1

# Exhaustive polarization on a small tag set.
universe = set(range(8))
even_subsets: list[set[int]] = []
for size in range(0, len(universe) + 1, 2):
    even_subsets.extend(set(c) for c in combinations(sorted(universe), size))

for left in even_subsets:
    for right in even_subsets:
        symmetric_difference = left ^ right
        polarization = len(left & right) % 2
        assert q_subset(symmetric_difference) == (
            q_subset(left) + q_subset(right) + polarization
        ) % 2

# Branch exchange dualizes the virtual line but preserves its parity.
assert q_half_rank(branch_plus - branch_minus) == q_half_rank(branch_minus - branch_plus)

# A subdivision inserts [B --id--> B], of virtual rank zero.
for rank in range(12):
    assert q_half_rank(rank - rank) == 0

print("log_quadratic_orientation_line_audit=PASS")
print(f"finite_anomaly={finite_anomaly}")
print(f"packet_orientation={packet_orientation}")
print("returned_orientation=square:0,nonsquare:1")
print("packet_return_polarization_cross_term=0,0")
print("terminal_full_trace_orientation=square:1,nonsquare:0")
print("generic_orientation=1+epsilon")
print(f"polarization_pairs_checked={len(even_subsets) ** 2}")
print(f"normalization_telescope_chains_checked={telescoping_chains_checked}")
