#!/usr/bin/env python3
"""Search and compare generated candidate sequences against the OEIS.

The OEIS has no automatic publication endpoint.  It does, however, expose
machine-readable JSON search results via:

    https://oeis.org/search?q=...&fmt=json

This script is deliberately read-only.  It can:

1. query the OEIS JSON search endpoint with conservative throttling and cache;
2. compare candidates against a local stripped or stripped.gz OEIS data file;
3. rank exact, shifted-prefix, and normalized matches;
4. prepare a review report for manual submission.

It does not log in, edit, or submit entries.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


IntTuple = tuple[int, ...]


def parse_terms(value: str | Sequence[int]) -> IntTuple:
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return ()
        return tuple(int(x.strip()) for x in text.replace(" ", ",").split(",") if x.strip())
    return tuple(int(x) for x in value)


def gcd_normalized(values: Sequence[int]) -> IntTuple:
    nonzero = [abs(x) for x in values if x]
    if not nonzero:
        return tuple(values)
    divisor = math.gcd(*nonzero)
    result = tuple(x // divisor for x in values)
    # Canonical sign: first nonzero term positive.
    first = next((x for x in result if x), 0)
    return tuple(-x for x in result) if first < 0 else result


def longest_common_prefix(a: Sequence[int], b: Sequence[int]) -> int:
    length = 0
    for x, y in zip(a, b):
        if x != y:
            break
        length += 1
    return length


def best_shifted_prefix(
    candidate: Sequence[int],
    target: Sequence[int],
    max_shift: int,
) -> tuple[int, int, str]:
    """Return (matched length, shift, orientation).

    shift >= 0 means target[shift:] is compared with candidate.
    orientation may be exact, negated, or gcd_normalized.
    """
    best_result = (0, 0, "exact")
    best_rank = (-1, -1, 0)
    orientation_priority = {"exact": 2, "gcd_normalized": 1, "negated": 0}
    variants = [
        ("exact", tuple(candidate)),
        ("gcd_normalized", gcd_normalized(candidate)),
        ("negated", tuple(-x for x in candidate)),
    ]
    target_variants = {
        "exact": tuple(target),
        "negated": tuple(target),
        "gcd_normalized": gcd_normalized(target),
    }
    for orientation, source in variants:
        comparison_target = target_variants[orientation]
        for shift in range(min(max_shift, len(comparison_target)) + 1):
            length = longest_common_prefix(source, comparison_target[shift:])
            rank = (length, orientation_priority[orientation], -shift)
            if rank > best_rank:
                best_rank = rank
                best_result = (length, shift, orientation)
    return best_result


def parse_oeis_data_field(data: object) -> IntTuple:
    if data is None:
        return ()
    if isinstance(data, list):
        return tuple(int(x) for x in data)
    text = str(data).strip().strip(",")
    if not text:
        return ()
    return tuple(int(x.strip()) for x in text.split(",") if x.strip())


def normalize_json_results(payload: object) -> list[dict[str, object]]:
    if isinstance(payload, dict):
        results = payload.get("results", [])
        if isinstance(results, list):
            return [r for r in results if isinstance(r, dict)]
        return []
    if isinstance(payload, list):
        return [r for r in payload if isinstance(r, dict)]
    return []


class OEISClient:
    def __init__(
        self,
        cache_dir: Path,
        delay: float = 2.0,
        timeout: float = 30.0,
        user_agent: str = (
            "modular-recurrence-candidate-review/1.0 "
            "(read-only OEIS comparison tool)"
        ),
    ) -> None:
        self.cache_dir = cache_dir
        self.delay = max(0.0, delay)
        self.timeout = timeout
        self.user_agent = user_agent
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._last_request = 0.0

    def search(self, query: str, refresh: bool = False) -> list[dict[str, object]]:
        cache_key = hashlib.sha256(query.encode("utf-8")).hexdigest()
        cache_path = self.cache_dir / f"{cache_key}.json"
        if cache_path.exists() and not refresh:
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            return normalize_json_results(payload)

        elapsed = time.monotonic() - self._last_request
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)

        url = "https://oeis.org/search?" + urllib.parse.urlencode(
            {"q": query, "fmt": "json"}
        )
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": self.user_agent,
                "Accept": "application/json",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                raw = response.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            raise RuntimeError(f"OEIS HTTP error {exc.code} for query {query!r}") from exc
        except urllib.error.URLError as exc:
            raise RuntimeError(f"OEIS network error for query {query!r}: {exc}") from exc
        self._last_request = time.monotonic()
        payload = json.loads(raw)
        cache_path.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        return normalize_json_results(payload)


@dataclass
class Match:
    candidate_id: str
    candidate_terms: list[int]
    query: str
    source: str
    oeis_number: str | None
    oeis_name: str | None
    oeis_terms: list[int]
    matched_prefix_length: int
    target_shift: int
    orientation: str
    exact_prefix: bool
    review_url: str | None


def oeis_number(record: dict[str, object]) -> str | None:
    number = record.get("number")
    if number is None:
        return None
    text = str(number)
    if text.upper().startswith("A"):
        return text.upper()
    try:
        return f"A{int(text):06d}"
    except ValueError:
        return text


def score_record(
    candidate_id: str,
    candidate_terms: Sequence[int],
    query: str,
    record: dict[str, object],
    source: str,
    max_shift: int,
) -> Match:
    terms = parse_oeis_data_field(record.get("data"))
    matched, shift, orientation = best_shifted_prefix(
        candidate_terms, terms, max_shift
    )
    number = oeis_number(record)
    return Match(
        candidate_id=candidate_id,
        candidate_terms=list(candidate_terms),
        query=query,
        source=source,
        oeis_number=number,
        oeis_name=str(record.get("name", "")) or None,
        oeis_terms=list(terms),
        matched_prefix_length=matched,
        target_shift=shift,
        orientation=orientation,
        exact_prefix=(orientation == "exact" and shift == 0 and matched == len(candidate_terms)),
        review_url=f"https://oeis.org/{number}" if number else None,
    )


def candidate_query(
    terms: Sequence[int],
    query_terms: int,
    skip: int,
) -> str:
    selected = tuple(terms)[skip : skip + query_terms]
    return ",".join(str(x) for x in selected)


def load_candidates(path: Path) -> list[tuple[str, IntTuple, dict[str, object]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload if isinstance(payload, list) else [payload]
    result = []
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            continue
        terms = parse_terms(record.get("terms", []))
        candidate_id = str(
            record.get(
                "candidate_id",
                f"{record.get('family', 'candidate')}:{record.get('start_one_based', index + 1)}",
            )
        )
        result.append((candidate_id, terms, record))
    return result


def open_stripped(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def iter_stripped(path: Path) -> Iterator[tuple[str, IntTuple]]:
    with open_stripped(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            try:
                number, data = line.split(" ", 1)
            except ValueError:
                continue
            number = number.strip()
            if not number.startswith("A"):
                continue
            try:
                terms = parse_oeis_data_field(data)
            except ValueError:
                continue
            yield number, terms


def offline_matches(
    candidate_id: str,
    candidate_terms: Sequence[int],
    stripped_path: Path,
    max_shift: int,
    min_match: int,
    limit: int,
) -> list[Match]:
    matches: list[Match] = []
    query = candidate_query(candidate_terms, min(12, len(candidate_terms)), 0)
    for number, terms in iter_stripped(stripped_path):
        matched, shift, orientation = best_shifted_prefix(
            candidate_terms, terms, max_shift
        )
        if matched < min_match:
            continue
        matches.append(
            Match(
                candidate_id=candidate_id,
                candidate_terms=list(candidate_terms),
                query=query,
                source="local_stripped",
                oeis_number=number,
                oeis_name=None,
                oeis_terms=list(terms),
                matched_prefix_length=matched,
                target_shift=shift,
                orientation=orientation,
                exact_prefix=(
                    orientation == "exact"
                    and shift == 0
                    and matched == len(candidate_terms)
                ),
                review_url=f"https://oeis.org/{number}",
            )
        )
    matches.sort(
        key=lambda m: (
            m.matched_prefix_length,
            m.exact_prefix,
            -m.target_shift,
        ),
        reverse=True,
    )
    return matches[:limit]


def online_matches(
    client: OEISClient,
    candidate_id: str,
    candidate_terms: Sequence[int],
    query_terms: int,
    skip: int,
    max_shift: int,
    limit: int,
    refresh: bool,
) -> list[Match]:
    query = candidate_query(candidate_terms, query_terms, skip)
    records = client.search(query, refresh=refresh)
    matches = [
        score_record(
            candidate_id,
            candidate_terms,
            query,
            record,
            "oeis_json",
            max_shift,
        )
        for record in records
    ]
    matches.sort(
        key=lambda m: (
            m.matched_prefix_length,
            m.exact_prefix,
            -m.target_shift,
        ),
        reverse=True,
    )
    return matches[:limit]


def manual_review_draft(
    metadata: dict[str, object],
    terms: Sequence[int],
    matches: Sequence[Match],
) -> str:
    coefficients = metadata.get("recurrence_coefficients_oldest_first", [])
    signature = metadata.get("oeis_signature_newest_first", [])
    seed = metadata.get("seed", [])
    period = metadata.get("period", "")
    start = metadata.get("start_one_based", "")
    family = metadata.get("family", "custom")
    lines = [
        f"Candidate: {family}, seed {seed}",
        f"Data: {','.join(map(str, terms))}",
        f"Recurrence coefficients (oldest first): {coefficients}",
        f"OEIS/LinearRecurrence signature (newest first): {signature}",
        f"Source period: {period}",
        f"Window start (one-based): {start}",
        f"Generating function: {metadata.get('generating_function_offset_1', '')}",
        "",
        "Closest OEIS matches:",
    ]
    if matches:
        for match in matches:
            lines.append(
                f"- {match.oeis_number or '?'}: match={match.matched_prefix_length}, "
                f"shift={match.target_shift}, orientation={match.orientation}, "
                f"name={match.oeis_name or ''}"
            )
    else:
        lines.append("- none returned by the selected comparison method")
    lines.extend(
        [
            "",
            "This draft is for human review only.",
            "Submission must be made manually through the OEIS contribution interface.",
        ]
    )
    return "\n".join(lines)


def write_csv(path: Path, matches: Sequence[Match]) -> None:
    fields = list(asdict(matches[0]).keys()) if matches else [
        "candidate_id",
        "candidate_terms",
        "query",
        "source",
        "oeis_number",
        "oeis_name",
        "oeis_terms",
        "matched_prefix_length",
        "target_shift",
        "orientation",
        "exact_prefix",
        "review_url",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for match in matches:
            row = asdict(match)
            row["candidate_terms"] = ",".join(map(str, match.candidate_terms))
            row["oeis_terms"] = ",".join(map(str, match.oeis_terms))
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Read-only OEIS search and comparison for generated candidates"
    )
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--mode", choices=("online", "offline"), default="online")
    parser.add_argument("--stripped", type=Path)
    parser.add_argument("--query-terms", type=int, default=10)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--max-shift", type=int, default=3)
    parser.add_argument("--min-match", type=int, default=6)
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--delay", type=float, default=2.0)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--cache-dir", type=Path, default=Path(".oeis_cache"))
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--csv-out", type=Path)
    parser.add_argument("--draft-dir", type=Path)
    args = parser.parse_args()

    if args.mode == "offline" and args.stripped is None:
        parser.error("--stripped is required in offline mode")
    candidates = load_candidates(args.candidates)
    client = (
        OEISClient(args.cache_dir, args.delay, args.timeout)
        if args.mode == "online"
        else None
    )
    all_matches: list[Match] = []
    reports: list[dict[str, object]] = []
    if args.draft_dir:
        args.draft_dir.mkdir(parents=True, exist_ok=True)

    for candidate_id, terms, metadata in candidates:
        if args.mode == "online":
            assert client is not None
            matches = online_matches(
                client,
                candidate_id,
                terms,
                args.query_terms,
                args.skip,
                args.max_shift,
                args.limit,
                args.refresh,
            )
        else:
            matches = offline_matches(
                candidate_id,
                terms,
                args.stripped,
                args.max_shift,
                args.min_match,
                args.limit,
            )
        all_matches.extend(matches)
        report = {
            "candidate_id": candidate_id,
            "metadata": metadata,
            "matches": [asdict(match) for match in matches],
        }
        reports.append(report)
        if args.draft_dir:
            safe_name = "".join(
                ch if ch.isalnum() or ch in "-_" else "_"
                for ch in candidate_id
            )
            (args.draft_dir / f"{safe_name}.txt").write_text(
                manual_review_draft(metadata, terms, matches) + "\n",
                encoding="utf-8",
            )

    output = json.dumps(reports, ensure_ascii=False, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(output + "\n", encoding="utf-8")
    else:
        print(output)
    if args.csv_out:
        write_csv(args.csv_out, all_matches)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
