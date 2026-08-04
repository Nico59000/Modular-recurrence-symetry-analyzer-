#!/usr/bin/env python3
"""Exact period and symmetry analysis for modular linear recurrences.

This module replaces finite-prefix/Mathematica-wrapper period guessing with
state-space cycle detection.  Coefficients are ordered from the oldest state
coordinate to the newest:

    u[n+k] = c[0] u[n] + c[1] u[n+1] + ... + c[k-1] u[n+k-1] (mod m).

The default modulus is 3, but the implementation works over Z/mZ.
"""

from __future__ import annotations

import argparse
import cmath
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable, Iterator, Sequence


IntTuple = tuple[int, ...]


def parse_int_tuple(text: str) -> IntTuple:
    text = text.strip()
    if not text:
        return ()
    if "," in text:
        parts = text.split(",")
    elif " " in text:
        parts = text.split()
    else:
        parts = list(text)
    try:
        return tuple(int(part.strip()) for part in parts if part.strip() != "")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid integer sequence: {text!r}") from exc


def canonical_residue(x: int, modulus: int) -> int:
    if modulus < 2:
        raise ValueError("modulus must be at least 2")
    return x % modulus


def negate_word(word: Sequence[int], modulus: int) -> IntTuple:
    return tuple((-x) % modulus for x in word)


def rotate_word(word: Sequence[int], shift: int) -> IntTuple:
    n = len(word)
    if n == 0:
        return ()
    shift %= n
    return tuple(word[shift:]) + tuple(word[:shift])


def canonical_rotation(word: Sequence[int]) -> IntTuple:
    """Lexicographically least cyclic rotation (Booth, O(n))."""
    s = tuple(word)
    n = len(s)
    if n == 0:
        return ()
    doubled = s + s
    i, j, k = 0, 1, 0
    while i < n and j < n and k < n:
        a, b = doubled[i + k], doubled[j + k]
        if a == b:
            k += 1
            continue
        if a > b:
            i = i + k + 1
            if i <= j:
                i = j + 1
        else:
            j = j + k + 1
            if j <= i:
                j = i + 1
        k = 0
    start = min(i, j)
    return doubled[start : start + n]


def minimal_period(word: Sequence[int]) -> IntTuple:
    """Return the primitive repeating block of an exactly periodic word."""
    s = tuple(word)
    n = len(s)
    if n <= 1:
        return s
    prefix = [0] * n
    j = 0
    for i in range(1, n):
        while j and s[i] != s[j]:
            j = prefix[j - 1]
        if s[i] == s[j]:
            j += 1
            prefix[i] = j
    p = n - prefix[-1]
    return s[:p] if n % p == 0 else s


def word_to_text(word: Sequence[int]) -> str:
    if all(0 <= x <= 9 for x in word):
        return "".join(str(x) for x in word)
    return "[" + ",".join(map(str, word)) + "]"


@dataclass(frozen=True)
class RecurrenceSpec:
    coefficients: IntTuple
    modulus: int = 3
    name: str = "custom"

    def __post_init__(self) -> None:
        if self.modulus < 2:
            raise ValueError("modulus must be at least 2")
        if not self.coefficients:
            raise ValueError("at least one coefficient is required")
        object.__setattr__(
            self,
            "coefficients",
            tuple(c % self.modulus for c in self.coefficients),
        )

    @property
    def order(self) -> int:
        return len(self.coefficients)

    @property
    def state_count(self) -> int:
        return self.modulus ** self.order

    @property
    def invertible(self) -> bool:
        # Determinant of the companion map is +/- the oldest coefficient.
        return math.gcd(self.coefficients[0], self.modulus) == 1


class PackedTransition:
    """Packed base-m state transition; oldest coordinate is least significant."""

    __slots__ = ("spec", "m", "k", "top_power")

    def __init__(self, spec: RecurrenceSpec) -> None:
        self.spec = spec
        self.m = spec.modulus
        self.k = spec.order
        self.top_power = self.m ** (self.k - 1)

    def encode(self, state: Sequence[int]) -> int:
        if len(state) != self.k:
            raise ValueError(f"seed must contain exactly {self.k} terms")
        value = 0
        power = 1
        for digit in state:
            value += (digit % self.m) * power
            power *= self.m
        return value

    def decode(self, packed: int) -> IntTuple:
        digits: list[int] = []
        x = packed
        for _ in range(self.k):
            digits.append(x % self.m)
            x //= self.m
        return tuple(digits)

    def first(self, packed: int) -> int:
        return packed % self.m

    def __call__(self, packed: int) -> int:
        x = packed
        total = 0
        for coefficient in self.spec.coefficients:
            digit = x % self.m
            x //= self.m
            total += coefficient * digit
        new_digit = total % self.m
        return packed // self.m + new_digit * self.top_power


def brent_cycle(
    transition: Callable[[int], int],
    initial_state: int,
    max_steps: int | None = None,
) -> tuple[int, int]:
    """Return (preperiod mu, cycle length lambda) with O(1) auxiliary memory."""
    power = lam = 1
    tortoise = initial_state
    hare = transition(initial_state)
    steps = 1
    while tortoise != hare:
        if max_steps is not None and steps > max_steps:
            raise RuntimeError("cycle search exceeded max_steps")
        if power == lam:
            tortoise = hare
            power *= 2
            lam = 0
        hare = transition(hare)
        lam += 1
        steps += 1

    mu = 0
    tortoise = hare = initial_state
    for _ in range(lam):
        hare = transition(hare)
    while tortoise != hare:
        if max_steps is not None and steps > 2 * max_steps:
            raise RuntimeError("preperiod search exceeded max_steps")
        tortoise = transition(tortoise)
        hare = transition(hare)
        mu += 1
        steps += 1
    return mu, lam


@dataclass
class CycleInfo:
    preperiod: int
    state_period: int
    output_period: int
    preperiod_word: IntTuple
    period_word: IntTuple
    cycle_start_state: IntTuple


def analyze_seed(spec: RecurrenceSpec, seed: Sequence[int]) -> CycleInfo:
    transition = PackedTransition(spec)
    initial = transition.encode(seed)
    mu, lam = brent_cycle(
        transition,
        initial,
        max_steps=2 * spec.state_count + 2,
    )
    state = initial
    preperiod_values: list[int] = []
    for _ in range(mu):
        preperiod_values.append(transition.first(state))
        state = transition(state)
    cycle_start = state
    period_values: list[int] = []
    for _ in range(lam):
        period_values.append(transition.first(state))
        state = transition(state)
    primitive = minimal_period(period_values)
    return CycleInfo(
        preperiod=mu,
        state_period=lam,
        output_period=len(primitive),
        preperiod_word=tuple(preperiod_values),
        period_word=primitive,
        cycle_start_state=transition.decode(cycle_start),
    )


def centered_mod3(word: Sequence[int]) -> IntTuple:
    mapping = {0: 0, 1: 1, 2: -1}
    if any(x % 3 not in mapping for x in word):
        raise ValueError("word is not a mod-3 word")
    return tuple(mapping[x % 3] for x in word)


def dft(values: Sequence[complex]) -> tuple[complex, ...]:
    n = len(values)
    if n == 0:
        return ()
    return tuple(
        sum(
            values[j] * cmath.exp(-2j * math.pi * k * j / n)
            for j in range(n)
        )
        for k in range(n)
    )


def affine_shift_symmetries(word: Sequence[int], modulus: int) -> list[dict[str, int]]:
    """Find w[n+h] = a*w[n] + b for units a modulo m."""
    s = tuple(x % modulus for x in word)
    n = len(s)
    units = [a for a in range(modulus) if math.gcd(a, modulus) == 1]
    symmetries: list[dict[str, int]] = []
    for shift in range(n):
        for a in units:
            for b in range(modulus):
                if all(s[(i + shift) % n] == (a * s[i] + b) % modulus for i in range(n)):
                    symmetries.append({"shift": shift, "a": a, "b": b})
    return symmetries


def affine_reversal_symmetries(word: Sequence[int], modulus: int) -> list[dict[str, int]]:
    """Find w[r-i] = a*w[i] + b for units a modulo m."""
    s = tuple(x % modulus for x in word)
    n = len(s)
    units = [a for a in range(modulus) if math.gcd(a, modulus) == 1]
    symmetries: list[dict[str, int]] = []
    for axis in range(n):
        for a in units:
            for b in range(modulus):
                if all(s[(axis - i) % n] == (a * s[i] + b) % modulus for i in range(n)):
                    symmetries.append({"axis": axis, "a": a, "b": b})
    return symmetries


@dataclass
class SymmetryReport:
    length: int
    primitive_word: IntTuple
    negated_word: IntTuple
    half_antiperiodic: bool
    half_length: int | None
    first_half: IntTuple | None
    second_half: IntTuple | None
    reversed_second_half: IntTuple | None
    complement_shifts: tuple[int, ...]
    globally_self_complementary_up_to_rotation: bool
    count_residues: dict[int, int]
    affine_shift_symmetries: list[dict[str, int]]
    affine_reversal_symmetries: list[dict[str, int]]
    mod3_even_frequency_max: float | None
    mod3_odd_frequency_max: float | None


def analyze_word(word: Sequence[int], modulus: int = 3) -> SymmetryReport:
    primitive = minimal_period(tuple(x % modulus for x in word))
    n = len(primitive)
    negated = negate_word(primitive, modulus)
    complement_shifts = tuple(
        h for h in range(n) if rotate_word(primitive, h) == negated
    )
    half = n // 2 if n % 2 == 0 else None
    half_antiperiodic = half is not None and half in complement_shifts
    first_half = primitive[:half] if half is not None else None
    second_half = primitive[half:] if half is not None else None
    reversed_second = tuple(reversed(second_half)) if second_half is not None else None
    counts = {r: primitive.count(r) for r in range(modulus)}

    even_max: float | None = None
    odd_max: float | None = None
    if modulus == 3 and n:
        spectrum = dft(tuple(complex(x) for x in centered_mod3(primitive)))
        even = [abs(spectrum[k]) for k in range(n) if k % 2 == 0]
        odd = [abs(spectrum[k]) for k in range(n) if k % 2 == 1]
        even_max = max(even, default=0.0)
        odd_max = max(odd, default=0.0)

    return SymmetryReport(
        length=n,
        primitive_word=primitive,
        negated_word=negated,
        half_antiperiodic=half_antiperiodic,
        half_length=half,
        first_half=first_half,
        second_half=second_half,
        reversed_second_half=reversed_second,
        complement_shifts=complement_shifts,
        globally_self_complementary_up_to_rotation=bool(complement_shifts),
        count_residues=counts,
        affine_shift_symmetries=affine_shift_symmetries(primitive, modulus),
        affine_reversal_symmetries=affine_reversal_symmetries(primitive, modulus),
        mod3_even_frequency_max=even_max,
        mod3_odd_frequency_max=odd_max,
    )


def mat_identity(n: int) -> list[list[int]]:
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def mat_mul(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]], modulus: int) -> list[list[int]]:
    rows, inner, cols = len(a), len(b), len(b[0])
    if len(a[0]) != inner:
        raise ValueError("matrix dimension mismatch")
    return [
        [sum(a[i][k] * b[k][j] for k in range(inner)) % modulus for j in range(cols)]
        for i in range(rows)
    ]


def mat_pow(matrix: Sequence[Sequence[int]], exponent: int, modulus: int) -> list[list[int]]:
    if exponent < 0:
        raise ValueError("negative matrix powers are not supported")
    result = mat_identity(len(matrix))
    base = [list(row) for row in matrix]
    e = exponent
    while e:
        if e & 1:
            result = mat_mul(result, base, modulus)
        base = mat_mul(base, base, modulus)
        e >>= 1
    return result


def companion_matrix(spec: RecurrenceSpec) -> list[list[int]]:
    k = spec.order
    matrix = [[0] * k for _ in range(k)]
    for i in range(k - 1):
        matrix[i][i + 1] = 1
    matrix[-1] = list(spec.coefficients)
    return matrix


def matrix_equal(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]]) -> bool:
    return all(tuple(x) == tuple(y) for x, y in zip(a, b))


def find_global_antiperiod(spec: RecurrenceSpec, max_h: int | None = None) -> int | None:
    """Least h>0 with M^h=-I, if found before max_h."""
    matrix = companion_matrix(spec)
    k = spec.order
    minus_identity = [
        [(-1 if i == j else 0) % spec.modulus for j in range(k)]
        for i in range(k)
    ]
    if max_h is None:
        max_h = spec.state_count
    power = mat_identity(k)
    for h in range(1, max_h + 1):
        power = mat_mul(power, matrix, spec.modulus)
        if matrix_equal(power, minus_identity):
            return h
    return None


def enumerate_cycles(spec: RecurrenceSpec, max_states: int = 2_000_000) -> list[IntTuple]:
    """Enumerate every state cycle exactly in O(m^k) time."""
    total = spec.state_count
    if total > max_states:
        raise ValueError(
            f"state space {total} exceeds max_states={max_states}; "
            "use single-seed Brent analysis or raise the limit"
        )
    transition = PackedTransition(spec)
    done = bytearray(total)
    cycles: list[IntTuple] = []
    for start in range(total):
        if done[start]:
            continue
        path: list[int] = []
        local_index: dict[int, int] = {}
        state = start
        while not done[state] and state not in local_index:
            local_index[state] = len(path)
            path.append(state)
            state = transition(state)
        if state in local_index:
            cycle_states = path[local_index[state] :]
            word = tuple(transition.first(s) for s in cycle_states)
            cycles.append(canonical_rotation(minimal_period(word)))
        for visited in path:
            done[visited] = 1
    cycles.sort(key=lambda w: (len(w), w))
    return cycles


def legacy_families(modulus: int = 3) -> tuple[RecurrenceSpec, ...]:
    families: list[RecurrenceSpec] = []
    names = {
        2: ["fibonacci"],
        3: ["padovan", "tribonacci"],
        4: ["duotetranacci", "tritetranacci", "tetranacci"],
        5: ["duopentanacci", "tripentanacci", "tetrapentanacci", "pentanacci"],
        6: ["duohexanacci", "trihexanacci", "tetrahexanacci", "pentahexanacci", "hexanacci"],
    }
    # Explicitly retain the repository's named recurrence families.
    coefficients = {
        "fibonacci": (1, 1),
        "padovan": (1, 1, 0),
        "tribonacci": (1, 1, 1),
        "duotetranacci": (1, 1, 0, 0),
        "tritetranacci": (1, 1, 1, 0),
        "tetranacci": (1, 1, 1, 1),
        "duopentanacci": (1, 1, 0, 0, 0),
        "tripentanacci": (1, 1, 1, 0, 0),
        "tetrapentanacci": (1, 1, 1, 1, 0),
        "pentanacci": (1, 1, 1, 1, 1),
        "duohexanacci": (1, 1, 0, 0, 0, 0),
        "trihexanacci": (1, 1, 1, 0, 0, 0),
        "tetrahexanacci": (1, 1, 1, 1, 0, 0),
        "pentahexanacci": (1, 1, 1, 1, 1, 0),
        "hexanacci": (1, 1, 1, 1, 1, 1),
    }
    for order in sorted(names):
        for name in names[order]:
            families.append(RecurrenceSpec(coefficients[name], modulus, name))
    return tuple(families)


def cycle_pair_classification(cycles: Sequence[IntTuple], modulus: int) -> list[dict[str, object]]:
    index = {canonical_rotation(w): i for i, w in enumerate(cycles)}
    records: list[dict[str, object]] = []
    for i, word in enumerate(cycles):
        report = analyze_word(word, modulus)
        complement = canonical_rotation(negate_word(word, modulus))
        partner = index.get(complement)
        records.append(
            {
                "cycle_id": i,
                "word": word_to_text(word),
                "length": len(word),
                "half_antiperiodic": report.half_antiperiodic,
                "complement_shifts": list(report.complement_shifts),
                "complement_partner_cycle_id": partner,
                "complement_partner_distinct": partner is not None and partner != i,
                "reversed_second_half": (
                    word_to_text(report.reversed_second_half)
                    if report.reversed_second_half is not None
                    else None
                ),
            }
        )
    return records


def analyze_family(spec: RecurrenceSpec, max_states: int) -> dict[str, object]:
    cycles = enumerate_cycles(spec, max_states=max_states)
    nonzero_cycles = [w for w in cycles if any(w)]
    return {
        "name": spec.name,
        "modulus": spec.modulus,
        "order": spec.order,
        "coefficients": list(spec.coefficients),
        "invertible": spec.invertible,
        "state_count": spec.state_count,
        "global_antiperiod_h": find_global_antiperiod(spec),
        "cycle_count": len(cycles),
        "nonzero_cycle_count": len(nonzero_cycles),
        "period_lengths": sorted({len(w) for w in cycles}),
        "cycles": cycle_pair_classification(cycles, spec.modulus),
    }


def _family_worker(payload: tuple[IntTuple, int, str, int]) -> dict[str, object]:
    coefficients, modulus, name, max_states = payload
    return analyze_family(RecurrenceSpec(coefficients, modulus, name), max_states)


def scan_binary_families(
    min_order: int,
    max_order: int,
    modulus: int,
    max_states: int,
    workers: int,
) -> list[dict[str, object]]:
    payloads: list[tuple[IntTuple, int, str, int]] = []
    for order in range(min_order, max_order + 1):
        for mask in range(1, 1 << order):
            coefficients = tuple((mask >> i) & 1 for i in range(order))
            name = f"binary-k{order}-coeff{''.join(map(str, coefficients))}"
            payloads.append((coefficients, modulus, name, max_states))
    if workers <= 1:
        return [_family_worker(payload) for payload in payloads]
    records: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {executor.submit(_family_worker, p): p for p in payloads}
        for future in as_completed(future_map):
            records.append(future.result())
    records.sort(key=lambda r: (r["order"], r["coefficients"]))
    return records


def json_ready(value: object) -> object:
    if hasattr(value, "__dataclass_fields__"):
        return {k: json_ready(v) for k, v in asdict(value).items()}
    if isinstance(value, tuple):
        return [json_ready(v) for v in value]
    if isinstance(value, list):
        return [json_ready(v) for v in value]
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    return value


def write_json(path: str | None, data: object) -> None:
    text = json.dumps(json_ready(data), ensure_ascii=False, indent=2, sort_keys=True)
    if path:
        Path(path).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


def command_word(args: argparse.Namespace) -> None:
    report = analyze_word(args.word, args.modulus)
    write_json(args.json_out, report)


def command_recurrence(args: argparse.Namespace) -> None:
    spec = RecurrenceSpec(args.coefficients, args.modulus, args.name)
    info = analyze_seed(spec, args.seed)
    result = {
        "spec": {
            "name": spec.name,
            "modulus": spec.modulus,
            "coefficients": list(spec.coefficients),
            "order": spec.order,
            "invertible": spec.invertible,
        },
        "cycle": info,
        "symmetry": analyze_word(info.period_word, spec.modulus),
        "global_antiperiod_h": find_global_antiperiod(spec, args.max_matrix_h),
    }
    write_json(args.json_out, result)


def command_cycles(args: argparse.Namespace) -> None:
    spec = RecurrenceSpec(args.coefficients, args.modulus, args.name)
    write_json(args.json_out, analyze_family(spec, args.max_states))


def command_legacy(args: argparse.Namespace) -> None:
    records = [analyze_family(spec, args.max_states) for spec in legacy_families(args.modulus)]
    write_json(args.json_out, records)
    if args.csv_out:
        with Path(args.csv_out).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "family",
                    "order",
                    "coefficients",
                    "invertible",
                    "global_antiperiod_h",
                    "cycle_id",
                    "word",
                    "length",
                    "half_antiperiodic",
                    "complement_partner_cycle_id",
                ]
            )
            for family in records:
                for cycle in family["cycles"]:
                    writer.writerow(
                        [
                            family["name"],
                            family["order"],
                            "".join(map(str, family["coefficients"])),
                            family["invertible"],
                            family["global_antiperiod_h"],
                            cycle["cycle_id"],
                            cycle["word"],
                            cycle["length"],
                            cycle["half_antiperiodic"],
                            cycle["complement_partner_cycle_id"],
                        ]
                    )


def command_scan(args: argparse.Namespace) -> None:
    records = scan_binary_families(
        args.min_order,
        args.max_order,
        args.modulus,
        args.max_states,
        args.workers,
    )
    write_json(args.json_out, records)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Exact modular linear-recurrence period and symmetry analyzer"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    word_parser = subparsers.add_parser("word", help="analyze an explicit periodic word")
    word_parser.add_argument("--word", type=parse_int_tuple, required=True)
    word_parser.add_argument("--modulus", type=int, default=3)
    word_parser.add_argument("--json-out")
    word_parser.set_defaults(func=command_word)

    recurrence_parser = subparsers.add_parser("recurrence", help="analyze one recurrence seed")
    recurrence_parser.add_argument("--coefficients", type=parse_int_tuple, required=True)
    recurrence_parser.add_argument("--seed", type=parse_int_tuple, required=True)
    recurrence_parser.add_argument("--modulus", type=int, default=3)
    recurrence_parser.add_argument("--name", default="custom")
    recurrence_parser.add_argument("--max-matrix-h", type=int)
    recurrence_parser.add_argument("--json-out")
    recurrence_parser.set_defaults(func=command_recurrence)

    cycles_parser = subparsers.add_parser("cycles", help="enumerate every state cycle of one recurrence")
    cycles_parser.add_argument("--coefficients", type=parse_int_tuple, required=True)
    cycles_parser.add_argument("--modulus", type=int, default=3)
    cycles_parser.add_argument("--name", default="custom")
    cycles_parser.add_argument("--max-states", type=int, default=2_000_000)
    cycles_parser.add_argument("--json-out")
    cycles_parser.set_defaults(func=command_cycles)

    legacy_parser = subparsers.add_parser("legacy", help="scan the fifteen recurrence families in the repository")
    legacy_parser.add_argument("--modulus", type=int, default=3)
    legacy_parser.add_argument("--max-states", type=int, default=2_000_000)
    legacy_parser.add_argument("--json-out")
    legacy_parser.add_argument("--csv-out")
    legacy_parser.set_defaults(func=command_legacy)

    scan_parser = subparsers.add_parser("scan", help="scan all nonzero binary coefficient masks")
    scan_parser.add_argument("--min-order", type=int, default=2)
    scan_parser.add_argument("--max-order", type=int, default=6)
    scan_parser.add_argument("--modulus", type=int, default=3)
    scan_parser.add_argument("--max-states", type=int, default=2_000_000)
    scan_parser.add_argument("--workers", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    scan_parser.add_argument("--json-out")
    scan_parser.set_defaults(func=command_scan)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
