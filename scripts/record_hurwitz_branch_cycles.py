#!/usr/bin/env python3
"""Record the completed Arb branch-cycle continuation runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "verification" / "hurwitz_branch_cycle_summary.json"


def loop(permutation, attempted, accepted, depth):
    return {
        "permutation_images": permutation,
        "attempted_tubes": attempted,
        "certified_tubes": accepted,
        "maximum_dyadic_depth": depth,
    }


RECORDS = [
    {
        "class_id": 1,
        "component": "sextic",
        "embedding_index": 3,
        "precision_bits": 256,
        "root_scale": "1/2",
        "loops": {
            "zero": loop([4,21,7,11,22,9,17,20,5,3,2,10,8,12,14,15,23,6,13,16,18,19,1], 23614, 11810, 14),
            "one": loop([20,7,3,9,13,15,2,17,4,10,11,12,5,14,6,16,8,18,21,1,19,23,22], 24880, 12443, 18),
            "infinity": loop([8,3,10,6,19,16,11,7,1,12,4,14,9,15,18,20,13,21,2,23,22,17,5], 21480, 10742, 18),
        },
    },
    {
        "class_id": 2,
        "component": "sextic",
        "embedding_index": 6,
        "precision_bits": 256,
        "root_scale": "1/2",
        "loops": {
            "zero": loop([22,19,2,9,4,3,1,6,17,23,8,11,20,7,12,14,16,15,13,5,10,21,18], 30972, 15489, 16),
            "one": loop([5,14,3,4,1,6,10,8,12,7,11,9,21,2,15,16,22,19,18,23,13,17,20], 14118, 7062, 13),
            "infinity": loop([20,16,6,5,7,8,21,11,15,14,12,4,22,3,18,17,1,2,23,10,19,9,13], 13382, 6693, 13),
        },
    },
    {
        "class_id": 3,
        "component": "sextic",
        "embedding_index": 5,
        "precision_bits": 256,
        "root_scale": "2",
        "loops": {
            "zero": loop([4,6,7,13,2,23,15,19,3,1,16,11,22,12,8,5,9,20,21,14,18,17,10], 10842, 5424, 12),
            "one": loop([13,2,16,19,22,6,21,10,9,8,11,12,1,23,17,3,15,18,4,20,7,5,14], 12344, 6175, 11),
            "infinity": loop([4,5,11,8,13,2,19,23,17,15,12,14,10,6,22,9,7,21,1,18,3,16,20], 13340, 6672, 12),
        },
    },
    {
        "class_id": 4,
        "component": "sextic",
        "embedding_index": 2,
        "precision_bits": 384,
        "root_scale": "1/4",
        "loops": {
            "zero": loop([18,15,4,1,6,7,22,3,5,20,12,9,17,10,13,19,16,14,8,11,2,23,21], 19654, 9830, 14),
            "one": loop([5,14,7,17,1,6,3,8,9,10,15,12,13,2,11,22,4,20,21,18,19,16,23], 6570, 3288, 10),
            "infinity": loop([9,18,6,13,4,5,8,19,12,14,2,11,15,21,20,7,3,10,23,1,16,17,22], 16918, 8461, 14),
        },
    },
    {
        "class_id": 5,
        "component": "sextic",
        "embedding_index": 4,
        "precision_bits": 256,
        "root_scale": "2",
        "loops": {
            "zero": loop([13,9,20,6,4,10,17,1,3,21,19,11,2,7,22,14,23,12,15,18,16,5,8], 6718, 3362, 10),
            "one": loop([3,2,1,4,10,17,13,8,19,5,11,22,7,14,20,16,6,18,9,15,23,12,21], 4900, 2453, 10),
            "infinity": loop([9,13,8,5,6,7,1,23,11,22,12,15,14,16,3,21,4,20,2,19,17,18,10], 17066, 8535, 14),
        },
    },
    {
        "class_id": 6,
        "component": "degree_one",
        "embedding_index": None,
        "precision_bits": 256,
        "root_scale": "32",
        "loops": {
            "zero": loop([13,7,8,2,22,3,10,21,5,14,4,6,11,20,12,19,1,9,17,23,18,16,15], 4730, 2368, 10),
            "one": loop([6,2,3,5,4,1,7,13,21,22,23,12,8,14,20,16,17,19,18,15,9,10,11], 16320, 8163, 15),
            "infinity": loop([12,4,6,9,11,17,2,1,8,5,20,15,3,10,14,22,19,16,21,23,18,7,13], 8860, 4432, 13),
        },
    },
    {
        "class_id": 7,
        "component": "sextic",
        "embedding_index": 1,
        "precision_bits": 256,
        "root_scale": "1",
        "loops": {
            "zero": loop([7,21,16,13,2,10,5,17,6,19,8,15,1,11,14,18,3,4,22,12,23,20,9], 7936, 3971, 11),
            "one": loop([10,3,2,16,5,17,18,14,9,1,15,12,13,8,11,4,6,7,19,20,21,23,22], 19748, 9877, 16),
            "infinity": loop([6,17,5,3,7,8,16,15,23,13,12,20,4,11,14,18,9,1,10,22,2,21,19], 5786, 2895, 12),
        },
    },
]


def compose(first, second):
    return [second[first[index] - 1] for index in range(23)]


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


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    identity = list(range(1, 24))
    for record in RECORDS:
        zero = record["loops"]["zero"]["permutation_images"]
        one = record["loops"]["one"]["permutation_images"]
        infinity = record["loops"]["infinity"]["permutation_images"]
        if cycle_lengths(zero) != [23]:
            raise AssertionError(f"class {record['class_id']} zero cycle shape")
        if cycle_lengths(one) != [1] * 7 + [2] * 8:
            raise AssertionError(f"class {record['class_id']} one cycle shape")
        if cycle_lengths(infinity) != [23]:
            raise AssertionError(f"class {record['class_id']} infinity cycle shape")
        if compose(compose(zero, one), infinity) != identity:
            raise AssertionError(f"class {record['class_id']} branch product")

    source_paths = [
        ROOT / "data" / "optimal_23_4_Z.json",
        ROOT / "data" / "hurwitz_monodromy_eliminant_candidate.json",
        ROOT / "scripts" / "compute_hurwitz_monodromy_resultants.py",
        ROOT / "scripts" / "assemble_hurwitz_monodromy_eliminant.py",
        ROOT / "scripts" / "certify_degree_one_branch_cycles.py",
        ROOT / "scripts" / "certify_hurwitz_branch_cycles.py",
        Path(__file__).resolve(),
    ]
    total_certified = sum(
        item["certified_tubes"]
        for record in RECORDS
        for item in record["loops"].values()
    )
    payload = {
        "schema": "m23.cover-investigation.branch-cycle-arb.v1",
        "status": "PASS_CERTIFIED_ARB_ROOT_CONTINUATION_ALL_SEVEN",
        "method": (
            "pairwise-disjoint complex interval-Newton contraction disks, "
            "uniform on each dyadic target-path segment"
        ),
        "base_target": ["1/2", "2"],
        "finite_loop_orientation": "counterclockwise",
        "infinity_loop_orientation": "clockwise",
        "sextic_finite_loop_radius": "2/5",
        "degree_one_finite_loop_radius": "1/4",
        "nielsen_convention": (
            "x=monodromy(one); y=monodromy(zero) on beta=N/(lambda*D), "
            "while the independently normalized degree-one optimal model "
            "uses y=inverse(monodromy(zero)) because its two order-23 "
            "endpoints are ordered oppositely"
        ),
        "records": RECORDS,
        "total_certified_tubes": total_certified,
        "maximum_dyadic_depth": max(
            item["maximum_dyadic_depth"]
            for record in RECORDS
            for item in record["loops"].values()
        ),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in source_paths
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("PASS recorded branch products and cycle shapes for all seven maps")
    print(f"PASS {total_certified} uniform Arb continuation tubes")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
