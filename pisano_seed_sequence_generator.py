#!/usr/bin/env python3
"""Generate integer recurrence candidates from cyclic windows of modular periods.

This tool implements the construction used by OEIS A276275:

    modular period: 1112201210010
    recurrence:     a(n+3) = a(n) + a(n+1)
    one-based start: 4
    seed:           220
    integer lift:   2,2,0,4,2,4,6,...

Coefficient convention
----------------------
Coefficients are ordered from the oldest state coordinate to the newest:

    a[n+k] = c[0] a[n] + c[1] a[n+1] + ... + c[k-1] a[n+k-1].

Thus the Padovan recurrence a(n)=a(n-2)+a(n-3) has coefficients (1,1,0).
The corresponding OEIS/Mathematica LinearRecurrence signature is reversed:
(0,1,1).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

try:
    import modular_recurrence_analyzer as mra
except ImportError as exc:  # pragma: no cover - user-facing error
    raise SystemExit(
        "modular_recurrence_analyzer.py must be in the same directory or on PYTHONPATH"
    ) from exc


IntTuple = tuple[int, ...]


def parse_int_tuple(text: str) -> IntTuple:
    return mra.parse_int_tuple(text)


def cyclic_window(period: Sequence[int], start: int, length: int) -> IntTuple:
    """Return a cyclic window; start is zero-based."""
    if not period:
        raise ValueError("period must not be empty")
    if length <= 0:
        raise ValueError("window length must be positive")
    n = len(period)
    return tuple(period[(start + j) % n] for j in range(length))


def linear_window(period: Sequence[int], start: int, length: int) -> IntTuple:
    if start < 0 or start + length > len(period):
        raise ValueError("linear window exceeds period boundaries")
    return tuple(period[start : start + length])


def generate_integer_recurrence(
    coefficients: Sequence[int],
    seed: Sequence[int],
    terms: int,
) -> IntTuple:
    coefficients = tuple(int(c) for c in coefficients)
    values = [int(x) for x in seed]
    k = len(coefficients)
    if k == 0:
        raise ValueError("coefficients must not be empty")
    if len(values) != k:
        raise ValueError(f"seed length {len(values)} does not match order {k}")
    if terms < 0:
        raise ValueError("terms must be nonnegative")
    if terms <= k:
        return tuple(values[:terms])
    while len(values) < terms:
        start = len(values) - k
        next_value = sum(coefficients[j] * values[start + j] for j in range(k))
        values.append(next_value)
    return tuple(values)


def validate_cyclic_period(
    period: Sequence[int],
    coefficients: Sequence[int],
    modulus: int,
) -> bool:
    """Check the recurrence at every cyclic index of a claimed period."""
    p = tuple(x % modulus for x in period)
    c = tuple(x % modulus for x in coefficients)
    n, k = len(p), len(c)
    if n == 0 or k == 0:
        return False
    for i in range(n):
        expected = sum(c[j] * p[(i + j) % n] for j in range(k)) % modulus
        if p[(i + k) % n] != expected:
            return False
    return True


def recurrence_denominator(coefficients: Sequence[int]) -> IntTuple:
    """Q(x) coefficients for A(x)=P(x)/Q(x), indexed by powers of x."""
    c = tuple(int(x) for x in coefficients)
    k = len(c)
    # Q(x) = 1 - c[k-1]x - c[k-2]x^2 - ... - c[0]x^k.
    return (1,) + tuple(-c[k - j] for j in range(1, k + 1))


def recurrence_numerator(
    coefficients: Sequence[int],
    seed: Sequence[int],
) -> IntTuple:
    """Numerator coefficients of the ordinary generating function at offset 0."""
    q = recurrence_denominator(coefficients)
    a = tuple(int(x) for x in seed)
    k = len(a)
    p = []
    for n in range(k):
        value = a[n]
        for j in range(1, min(n, len(q) - 1) + 1):
            value += q[j] * a[n - j]
        p.append(value)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def polynomial_text(coefficients: Sequence[int], variable: str = "x") -> str:
    parts: list[str] = []
    for power, coefficient in enumerate(coefficients):
        if coefficient == 0:
            continue
        magnitude = abs(coefficient)
        if power == 0:
            body = str(magnitude)
        elif power == 1:
            body = variable if magnitude == 1 else f"{magnitude}*{variable}"
        else:
            body = (
                f"{variable}^{power}"
                if magnitude == 1
                else f"{magnitude}*{variable}^{power}"
            )
        if not parts:
            parts.append(body if coefficient > 0 else f"-{body}")
        else:
            parts.append((" + " if coefficient > 0 else " - ") + body)
    return "".join(parts) if parts else "0"


def generating_function_text(
    coefficients: Sequence[int],
    seed: Sequence[int],
    offset: int = 0,
) -> str:
    numerator = recurrence_numerator(coefficients, seed)
    denominator = recurrence_denominator(coefficients)
    numerator_text = polynomial_text(numerator)
    if offset:
        numerator_text = f"x^{offset}*({numerator_text})"
    return f"({numerator_text})/({polynomial_text(denominator)})"


def gcd_normalized(values: Sequence[int]) -> IntTuple:
    nonzero = [abs(x) for x in values if x]
    if not nonzero:
        return tuple(values)
    divisor = math.gcd(*nonzero)
    return tuple(x // divisor for x in values)


def candidate_query(terms: Sequence[int], count: int = 12, skip: int = 0) -> str:
    selected = tuple(terms)[skip : skip + count]
    return ",".join(str(x) for x in selected)


@dataclass
class Candidate:
    family: str
    period: str
    period_length: int
    modulus: int
    recurrence_coefficients_oldest_first: list[int]
    oeis_signature_newest_first: list[int]
    start_zero_based: int
    start_one_based: int
    seed: list[int]
    seed_text: str
    cyclic: bool
    period_valid_for_recurrence: bool
    terms: list[int]
    gcd_normalized_terms: list[int]
    generating_function_offset_0: str
    generating_function_offset_1: str
    oeis_search_query: str
    provenance_comment: str
    manual_submission_name: str
    manual_submission_formula: str


def build_candidate(
    period: Sequence[int],
    coefficients: Sequence[int],
    modulus: int,
    start: int,
    terms: int,
    family: str,
    cyclic: bool,
) -> Candidate:
    order = len(coefficients)
    seed = (
        cyclic_window(period, start, order)
        if cyclic
        else linear_window(period, start, order)
    )
    generated = generate_integer_recurrence(coefficients, seed, terms)
    period_text = mra.word_to_text(period)
    seed_text = mra.word_to_text(seed)
    signature = list(reversed(tuple(coefficients)))
    recurrence_formula = " + ".join(
        f"{coefficient}*a(n-{order-j})"
        for j, coefficient in enumerate(coefficients)
        if coefficient
    ) or "0"
    return Candidate(
        family=family,
        period=period_text,
        period_length=len(period),
        modulus=modulus,
        recurrence_coefficients_oldest_first=list(coefficients),
        oeis_signature_newest_first=signature,
        start_zero_based=start,
        start_one_based=start + 1,
        seed=list(seed),
        seed_text=seed_text,
        cyclic=cyclic,
        period_valid_for_recurrence=validate_cyclic_period(
            period, coefficients, modulus
        ),
        terms=list(generated),
        gcd_normalized_terms=list(gcd_normalized(generated)),
        generating_function_offset_0=generating_function_text(
            coefficients, seed, offset=0
        ),
        generating_function_offset_1=generating_function_text(
            coefficients, seed, offset=1
        ),
        oeis_search_query=candidate_query(generated),
        provenance_comment=(
            f"Generated from the length-{order} window beginning at term "
            f"{start + 1} of the modulo-{modulus} period {period_text}; "
            "the lifted integer sequence satisfies the same linear recurrence."
        ),
        manual_submission_name=(
            f"{family}-type sequence with seed {seed_text} extracted from "
            f"a modulo-{modulus} period"
        ),
        manual_submission_formula=f"a(n) = {recurrence_formula}",
    )


def all_candidates(
    period: Sequence[int],
    coefficients: Sequence[int],
    modulus: int,
    terms: int,
    family: str,
    cyclic: bool,
    deduplicate: str,
) -> list[Candidate]:
    limit = len(period) if cyclic else len(period) - len(coefficients) + 1
    if limit < 1:
        raise ValueError("period is shorter than recurrence order")
    candidates: list[Candidate] = []
    seen: set[tuple[int, ...]] = set()
    for start in range(limit):
        candidate = build_candidate(
            period,
            coefficients,
            modulus,
            start,
            terms,
            family,
            cyclic,
        )
        if deduplicate == "seed":
            key = tuple(candidate.seed)
        elif deduplicate == "terms":
            key = tuple(candidate.terms)
        elif deduplicate == "normalized":
            key = tuple(candidate.gcd_normalized_terms)
        else:
            key = (start,)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(candidate)
    return candidates


def write_json(path: str | None, data: object) -> None:
    text = json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True)
    if path:
        Path(path).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


def write_csv(path: str, candidates: Sequence[Candidate]) -> None:
    fields = [
        "family",
        "period",
        "modulus",
        "start_one_based",
        "seed_text",
        "recurrence_coefficients_oldest_first",
        "oeis_signature_newest_first",
        "period_valid_for_recurrence",
        "terms",
        "gcd_normalized_terms",
        "generating_function_offset_0",
        "oeis_search_query",
    ]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for candidate in candidates:
            row = asdict(candidate)
            row["recurrence_coefficients_oldest_first"] = ",".join(
                map(str, candidate.recurrence_coefficients_oldest_first)
            )
            row["oeis_signature_newest_first"] = ",".join(
                map(str, candidate.oeis_signature_newest_first)
            )
            row["terms"] = ",".join(map(str, candidate.terms))
            row["gcd_normalized_terms"] = ",".join(
                map(str, candidate.gcd_normalized_terms)
            )
            writer.writerow({field: row[field] for field in fields})


def command_single(args: argparse.Namespace) -> None:
    period = args.period
    start = args.start - 1 if args.one_based else args.start
    candidate = build_candidate(
        period,
        args.coefficients,
        args.modulus,
        start,
        args.terms,
        args.family,
        not args.linear_window,
    )
    write_json(args.json_out, asdict(candidate))


def command_all(args: argparse.Namespace) -> None:
    candidates = all_candidates(
        args.period,
        args.coefficients,
        args.modulus,
        args.terms,
        args.family,
        not args.linear_window,
        args.deduplicate,
    )
    data = [asdict(candidate) for candidate in candidates]
    write_json(args.json_out, data)
    if args.csv_out:
        write_csv(args.csv_out, candidates)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Lift recurrence-order windows of modular periods to integer "
            "linear-recurrence candidate sequences"
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--period", type=parse_int_tuple, required=True)
        p.add_argument("--coefficients", type=parse_int_tuple, required=True)
        p.add_argument("--modulus", type=int, default=3)
        p.add_argument("--terms", type=int, default=40)
        p.add_argument("--family", default="custom")
        p.add_argument(
            "--linear-window",
            action="store_true",
            help="do not wrap windows around the end of the period",
        )
        p.add_argument("--json-out")

    single = subparsers.add_parser(
        "single", help="generate one candidate from one period window"
    )
    add_common(single)
    single.add_argument("--start", type=int, required=True)
    single.add_argument(
        "--one-based",
        action="store_true",
        help="interpret --start as one-based (OEIS prose convention)",
    )
    single.set_defaults(func=command_single)

    all_parser = subparsers.add_parser(
        "all", help="generate candidates from every recurrence-order window"
    )
    add_common(all_parser)
    all_parser.add_argument(
        "--deduplicate",
        choices=("none", "seed", "terms", "normalized"),
        default="seed",
    )
    all_parser.add_argument("--csv-out")
    all_parser.set_defaults(func=command_all)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if len(args.coefficients) == 0:
        parser.error("--coefficients must not be empty")
    if args.modulus < 2:
        parser.error("--modulus must be at least 2")
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
