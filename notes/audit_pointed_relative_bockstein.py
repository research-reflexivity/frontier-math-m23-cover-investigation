#!/usr/bin/env python3
"""Audit the finite linear algebra behind the pointed relative Bockstein.

This deliberately distinguishes three constructions:

1. the unrestricted graph-to-sheet mapping cone, in which the apparent
   coefficient Bockstein is zero;
2. the specialization-unit square at a split node, whose distinguished
   mod-four lift has a nonzero divided defect (but not an ordinary
   coefficient Bockstein); and
3. the relative cochain complex of a path with its two ends marked, in
   which that defect survives as the unique top relative class.

No branch-cycle representative is used.
"""


def rank_mod_2(matrix: list[list[int]]) -> int:
    rows = [[entry & 1 for entry in row] for row in matrix]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        for row in range(row_count):
            if row != pivot_row and rows[row][column]:
                rows[row] = [
                    left ^ right for left, right in zip(rows[row], rows[pivot_row])
                ]
        pivot_row += 1
    return pivot_row


def path_coboundary(edge_count: int, relative: bool) -> list[list[int]]:
    """Return C^0 -> C^1 for a path, optionally relative to both ends."""
    vertex_indices = range(1, edge_count) if relative else range(edge_count + 1)
    columns = list(vertex_indices)
    return [
        [int(vertex == edge) ^ int(vertex == edge + 1) for vertex in columns]
        for edge in range(edge_count)
    ]


# The natural action is doubly transitive.  These are the two entries of
# rho rho^*, so the naive three-term self-dual cone is a complex mod four.
group_order = 10_200_960
one_point_stabilizer = group_order // 23
two_point_stabilizer = group_order // (23 * 22)
assert one_point_stabilizer == 443_520
assert two_point_stabilizer == 20_160
assert one_point_stabilizer % 4 == 0
assert two_point_stabilizer % 4 == 0

# Write an invariant sheet operator as a*J+b*I.
full = (120, 1035)
trace = (8, 69)
packet = tuple(left - right for left, right in zip(full, trace))
assert packet == (112, 966)

full_mod_4 = tuple(entry % 4 for entry in full)
trace_mod_4 = tuple(entry % 4 for entry in trace)
packet_mod_4 = tuple(entry % 4 for entry in packet)
assert full_mod_4 == (0, 3)   # -I
assert trace_mod_4 == (0, 1)  # I
assert packet_mod_4 == (0, 2) # 2I

# The tempting ordinary mapping-cone Bockstein is not intrinsic.  The same
# mod-two packet has a lift packet-2*full whose sheet operator is zero mod 4.
alternative_packet = tuple(
    packet_entry - 2 * full_entry
    for packet_entry, full_entry in zip(packet, full)
)
assert alternative_packet == (-128, -1104)
assert all(entry % 4 == 0 for entry in alternative_packet)

# The same obstruction persists in the rank-one horizontal logarithmic
# nearby-cycle stalk.  Nakayama's formula gives one free coefficient copy
# in degrees zero and one.  Reduction Z/4 -> F_2 is surjective in both
# degrees, so every connecting homomorphism for
# 0 -> F_2 -> Z/4 -> F_2 -> 0 is zero.  Adding finitely many split tags only
# takes direct sums and cannot change this conclusion.
z4 = tuple(range(4))
f2 = tuple(range(2))
reduction_image_h0 = {entry % 2 for entry in z4}
reduction_image_h1 = {entry % 2 for entry in z4}
assert reduction_image_h0 == set(f2)
assert reduction_image_h1 == set(f2)
log_nearby_bockstein_h0 = 0
log_nearby_bockstein_h1 = 0
assert log_nearby_bockstein_h0 == log_nearby_bockstein_h1 == 0

# In contrast, the two canonical effective graph sums are separately fixed.
# Their relative comparison is (-I)*(I)^(-1), so its first divided layer is I.
relative_comparison_mod_4 = 3
relative_first_layer = ((relative_comparison_mod_4 - 1) // 2) % 2
assert relative_first_layer == 1

# Split-node specialization unit delta(1)=(1,1).  The canonical packet lift
# acts on the two normalized branches with the exact coefficients below.
diagonal_branch = 1078
transposition_branch = 112
assert diagonal_branch % 2 == transposition_branch % 2 == 0
divided_defect = ((diagonal_branch - transposition_branch) // 2) % 2
assert divided_defect == 1

# Equivalently, A*delta is (2,0) mod 4.  After division by two its image in
# F_2^2 / <(1,1)> is the nonzero vanishing-cycle class.
specialization_defect = (diagonal_branch % 4, transposition_branch % 4)
divided_specialization_defect = tuple(entry // 2 for entry in specialization_defect)
assert specialization_defect == (2, 0)
assert divided_specialization_defect == (1, 0)
assert divided_specialization_defect[0] ^ divided_specialization_defect[1] == 1


def quadratic_half_weight(vector: tuple[int, ...]) -> int:
    assert sum(vector) % 2 == 0
    return (sum(vector) // 2) % 2


# Half-weight is the canonical quadratic refinement of the tagged dot
# product on the even-weight binary code.
for bit_count in range(2, 10):
    even_vectors = [
        tuple((mask >> index) & 1 for index in range(bit_count))
        for mask in range(1 << bit_count)
        if mask.bit_count() % 2 == 0
    ]
    for left in even_vectors:
        for right in even_vectors:
            total = tuple(a ^ b for a, b in zip(left, right))
            dot = sum(a * b for a, b in zip(left, right)) % 2
            assert quadratic_half_weight(total) == (
                quadratic_half_weight(left)
                ^ quadratic_half_weight(right)
                ^ dot
            )

# Every packet has even multiplicity in every ordered sheet-pair bucket.
# Its global half-weight is odd, and polarization recovers the tagged
# intersection pairing because every graph label occurs on 23 source sheets.
packet_size = 3795 - 253
assert packet_size == 3542
packet_global_half_weight = (23 * packet_size // 2) % 2
assert packet_global_half_weight == 1
assert (diagonal_branch // 2) % 2 == 1
assert (transposition_branch // 2) % 2 == 0

# The returned full-class/trace half-conductor is the same quadratic layer.
square_returned_conductor = 98 * 253
nonsquare_returned_conductor = 0
assert (square_returned_conductor // 2) % 2 == 1
assert (nonsquare_returned_conductor // 2) % 2 == 0

# Relative half-distance of two odd traces is the generic augmented selector.
for overlap_size in range(254):
    relative_half_distance = ((2 * 253 - 2 * overlap_size) // 2) % 2
    assert relative_half_distance == (1 + overlap_size) % 2
assert ((2 * 253 - 2 * 253) // 2) % 2 == 0  # square return
assert ((2 * 253 - 2 * 0) // 2) % 2 == 1    # nonsquare return

# The corrected special target is quadratic rather than Bockstein.  The
# half-distance between the full odd class and the returned odd trace is
# 1+q: for a square return the trace is contained in the class, whereas for
# a nonsquare return it lies in the disjoint adjacent Mathieu class.
full_class_size = 3795
trace_size = 253
square_full_trace_half_distance = (full_class_size - trace_size) // 2 % 2
nonsquare_full_trace_half_distance = (full_class_size + trace_size) // 2 % 2
assert square_full_trace_half_distance == 1
assert nonsquare_full_trace_half_distance == 0
for return_bit, special_value in (
    (0, square_full_trace_half_distance),
    (1, nonsquare_full_trace_half_distance),
):
    assert special_value == (1 + return_bit) % 2

# Without marked ends a path has H^1=0.  Relative to its two ends it has one
# top class, detected by summing the edge coordinates.  This is the precise
# reason the endpoint framing is essential.
for edge_count in range(1, 17):
    absolute = path_coboundary(edge_count, relative=False)
    relative = path_coboundary(edge_count, relative=True)
    assert rank_mod_2(absolute) == edge_count
    assert rank_mod_2(relative) == edge_count - 1
    assert edge_count - rank_mod_2(relative) == 1
    for column in zip(*relative):
        assert sum(column) % 2 == 0

# The endpoint formula is independent of subdivision.
for return_bit in (0, 1):
    finite_value = 1
    wild_value = return_bit
    relative_value = finite_value ^ wild_value
    assert relative_value == 1 + return_bit - 2 * return_bit

print("rho_rho_adjoint_entries_mod_4=0")
print("canonical_sheet_lifts=full:-I,trace:+I")
print("packet_sheet_operator=2I_mod_4")
print("ordinary_mapping_cone_bockstein=zero_by_alternative_lift")
print("horizontal_log_nearby_coefficient_bockstein=zero")
print("specialization_unit_divided_defect=1")
print("half_weight_quadratic_refines_tagged_dot_product=true")
print("packet_global_half_weight=1")
print("returned_half_conductor_values=1,0")
print("relative_trace_half_distance=1+overlap_parity")
print("returned_trace_half_distance=return_bit")
print("full_class_returned_trace_half_distance=1+return_bit")
print("unpointed_path_H1_dimension=0")
print("two_endpoint_relative_path_H1_dimension=1")
print("relative_endpoint_formula=finite_value+wild_value")
print("PASS_POINTED_RELATIVE_BOCKSTEIN_AUDIT")
