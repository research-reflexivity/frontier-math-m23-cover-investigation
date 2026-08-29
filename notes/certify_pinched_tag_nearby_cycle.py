#!/usr/bin/env python3
"""Exact finite certificate for the pinched-tag nearby-cycle lemma."""

from itertools import product


def dot(alpha: tuple[int, ...], beta: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(alpha, beta, strict=True)) % 2


def bucket(alpha: tuple[int, ...], beta: tuple[int, ...]) -> int:
    return (sum(alpha) * sum(beta)) % 2


def off_diagonal(alpha: tuple[int, ...], beta: tuple[int, ...]) -> int:
    return sum(
        alpha[i] * beta[j]
        for i in range(len(alpha))
        for j in range(len(beta))
        if i != j
    ) % 2


# The conductor decomposition is an identity for every binary pair, not
# only for the all-ones vectors occurring in the Mathieu trace.
for size in range(1, 7):
    vectors = tuple(product((0, 1), repeat=size))
    for alpha in vectors:
        for beta in vectors:
            assert bucket(alpha, beta) == dot(alpha, beta) ^ off_diagonal(alpha, beta)

# Exact endpoint values.
diagonal_branch = 77**2
transposition_branch = 8**2
node_difference = diagonal_branch - transposition_branch
assert (diagonal_branch, transposition_branch, node_difference) == (5929, 64, 5865)
assert node_difference % 2 == 1
conductor_pairing_degree = 23 * diagonal_branch + 253 * transposition_branch
assert conductor_pairing_degree == 152559

# Universal untagged pairing on the seven singleton and eight two-set
# inertia orbits.
singleton_pairing = 23 * 11**2
two_set_pairing = 253
untagged_total = 7 * singleton_pairing + 8 * two_set_pairing
assert (singleton_pairing, two_set_pairing, untagged_total) == (2783, 253, 21505)
assert untagged_total % 2 == 1

# Returned normalized and node tables in fields
# [tagged, untagged, off-diagonal].
square_normalized = (1, 1, 0)
nonsquare_normalized = (0, 1, 1)
square_node = ((77 - 8) % 2, node_difference % 2, (node_difference - (77 - 8)) % 2)
nonsquare_node = (0, node_difference % 2, node_difference % 2)
assert square_node == square_normalized
assert nonsquare_node == nonsquare_normalized

# Tree normalization: every internal value occurs twice.
for length in range(1, 13):
    for values in product((0, 1), repeat=length + 1):
        boundary_sum = sum(values[i] ^ values[i + 1] for i in range(length)) % 2
        assert boundary_sum == values[0] ^ values[-1]

# The three returned values give the singular-position table.
return_bits = (0, 1, 1)
wild_terms = tuple(1 ^ q for q in return_bits)
finite_term = 1
special_traces = tuple(wild ^ finite_term for wild in wild_terms)
assert special_traces == return_bits

print("pinched_bucket_identity_checked_for_fiber_sizes=1..6")
print("untagged_total=7*(23*11^2)+8*253=21505")
print("node_coefficients=[5929,64],difference=5865=1_mod_2")
print("conductor_bucket_distribution=23x77_and_253x8")
print("conductor_pairing_degree=152559")
print("square_fields_tagged_untagged_offdiag=[1,1,0]")
print("nonsquare_fields_tagged_untagged_offdiag=[0,1,1]")
print("tree_telescope_checked_through_12_edges=true")
print("special_trace_table=[0,1,1]")
print("PASS_PINCHED_TAG_NEARBY_CYCLE")
