#!/usr/bin/env sage-python
"""Verify the recorded Arb branch cycles and their Nielsen-class labels."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from sage.all import Permutation, PermutationGroup


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATH = ROOT / "verification" / "hurwitz_branch_cycle_summary.json"


def compose(first, second):
    return [second[first[index] - 1] for index in range(23)]


def inverse(permutation):
    result = [0] * 23
    for source, target in enumerate(permutation, 1):
        result[target - 1] = source
    return result


def conjugate(permutation, relabeling):
    return compose(compose(inverse(relabeling), permutation), relabeling)


def cycle_permutation(*cycles):
    result = list(range(1, 24))
    for cycle in cycles:
        for index, value in enumerate(cycle):
            result[value - 1] = cycle[(index + 1) % len(cycle)]
    return result


FIXED_Y = cycle_permutation(
    (1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,13,7,12)
)
REPRESENTATIVES = [
    cycle_permutation((4,16),(5,10),(6,21),(8,19),(9,18),(11,23),(13,14),(17,22)),
    cycle_permutation((4,19),(5,9),(6,17),(8,16),(10,18),(11,14),(13,23),(21,22)),
    cycle_permutation((4,6),(5,9),(7,21),(8,15),(11,19),(12,18),(13,20),(14,17)),
    cycle_permutation((3,17),(5,11),(7,18),(8,16),(9,21),(12,19),(14,22),(15,23)),
    cycle_permutation((3,11),(5,17),(6,20),(7,18),(8,19),(10,13),(12,16),(14,22)),
    cycle_permutation((3,21),(4,16),(9,22),(10,15),(11,20),(12,19),(13,18),(14,17)),
    cycle_permutation((3,8),(5,12),(6,15),(9,10),(11,16),(13,21),(17,19),(20,23)),
]


def cycle_lengths(permutation):
    seen = set()
    result = []
    for start in range(1, 24):
        if start in seen:
            continue
        value = start
        length = 0
        while value not in seen:
            seen.add(value)
            length += 1
            value = permutation[value - 1]
        result.append(length)
    return sorted(result)


def nielsen_matches(x, y):
    matches = set()
    for target_start in range(1, 24):
        relabeling = [0] * 23
        source = 1
        target = target_start
        for _ in range(23):
            relabeling[source - 1] = target
            source = y[source - 1]
            target = FIXED_Y[target - 1]
        transformed_x = conjugate(x, relabeling)
        power = list(range(1, 24))
        for class_id, representative in enumerate(REPRESENTATIVES, 1):
            for _ in range(23):
                if conjugate(representative, power) == transformed_x:
                    matches.add(class_id)
                power = compose(power, FIXED_Y)
    return matches


def identify_nielsen_class(x, y):
    matches = nielsen_matches(x, y)
    if len(matches) != 1:
        raise AssertionError(f"ambiguous Nielsen match: {sorted(matches)}")
    return matches.pop()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary = json.loads(SUMMARY_PATH.read_text())
    assert summary["schema"] == "m23.cover-investigation.branch-cycle-arb.v1"
    assert summary["status"] == "PASS_CERTIFIED_ARB_ROOT_CONTINUATION_ALL_SEVEN"
    assert [record["class_id"] for record in summary["records"]] == list(range(1, 8))
    assert [
        record["class_id"]
        for record in summary["records"]
        if record["component"] == "sextic"
    ] == [1, 2, 3, 4, 5, 7]
    assert [
        (record["embedding_index"], record["class_id"])
        for record in summary["records"]
        if record["component"] == "sextic"
    ] == [(3,1),(6,2),(5,3),(2,4),(4,5),(1,7)]

    identity = list(range(1, 24))
    total = 0
    maximum_depth = 0
    for record in summary["records"]:
        zero = record["loops"]["zero"]["permutation_images"]
        one = record["loops"]["one"]["permutation_images"]
        infinity = record["loops"]["infinity"]["permutation_images"]
        assert cycle_lengths(zero) == [23]
        assert cycle_lengths(one) == [1] * 7 + [2] * 8
        assert cycle_lengths(infinity) == [23]
        assert compose(compose(zero, one), infinity) == identity
        x = one
        y = zero if record["component"] == "sextic" else inverse(zero)
        assert identify_nielsen_class(x, y) == record["class_id"]
        group = PermutationGroup([Permutation(x), Permutation(y)])
        assert group.order() == 10200960
        for loop_name, item in record["loops"].items():
            assert item["attempted_tubes"] >= item["certified_tubes"] > 0
            original_segments = 4 if loop_name == "infinity" else 6
            assert (
                item["attempted_tubes"]
                == 2 * item["certified_tubes"] - original_segments
            )
            total += item["certified_tubes"]
            maximum_depth = max(maximum_depth, item["maximum_dyadic_depth"])

    assert total == summary["total_certified_tubes"] == 150145
    assert maximum_depth == summary["maximum_dyadic_depth"] == 18
    for relative, expected in summary["source_sha256"].items():
        assert sha256(ROOT / relative) == expected

    print("PASS Arb branch products and cycle shapes for all seven exact maps")
    print("PASS exact match with Nielsen classes 1 through 7")
    print("PASS every branch-cycle pair generates M23 of order 10200960")


if __name__ == "__main__":
    main()
