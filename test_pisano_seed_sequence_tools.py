import json
import tempfile
import unittest
from pathlib import Path

import pisano_seed_sequence_generator as gen
import oeis_sequence_compare as osc


class SeedGeneratorTests(unittest.TestCase):
    def test_a276275_seed_and_terms(self):
        period = tuple(int(x) for x in "1112201210010")
        candidate = gen.build_candidate(
            period=period,
            coefficients=(1, 1, 0),
            modulus=3,
            start=3,  # fourth term, zero-based
            terms=15,
            family="padovan",
            cyclic=True,
        )
        self.assertEqual(candidate.seed, [2, 2, 0])
        self.assertEqual(
            candidate.terms,
            [2, 2, 0, 4, 2, 4, 6, 6, 10, 12, 16, 22, 28, 38, 50],
        )
        self.assertTrue(candidate.period_valid_for_recurrence)
        self.assertEqual(candidate.oeis_signature_newest_first, [0, 1, 1])
        self.assertEqual(
            candidate.generating_function_offset_0,
            "(2 + 2*x - 2*x^2)/(1 - x^2 - x^3)",
        )

    def test_all_windows(self):
        candidates = gen.all_candidates(
            tuple(int(x) for x in "1112201210010"),
            (1, 1, 0),
            3,
            20,
            "padovan",
            True,
            "seed",
        )
        seeds = {tuple(candidate.seed) for candidate in candidates}
        self.assertIn((2, 2, 0), seeds)
        self.assertIn((2, 1, 0), seeds)


class OEISCompareTests(unittest.TestCase):
    def test_json_result_normalization_and_match(self):
        payload = {
            "count": 1,
            "results": [
                {
                    "number": 276275,
                    "data": "2,2,0,4,2,4,6,6,10,12,16,22",
                    "name": "Padovan like sequence",
                }
            ],
        }
        records = osc.normalize_json_results(payload)
        match = osc.score_record(
            "padovan:4",
            (2, 2, 0, 4, 2, 4, 6, 6, 10, 12),
            "2,2,0,4,2,4,6,6,10,12",
            records[0],
            "test",
            3,
        )
        self.assertEqual(match.oeis_number, "A276275")
        self.assertEqual(match.matched_prefix_length, 10)
        self.assertTrue(match.exact_prefix)

    def test_offline_stripped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stripped"
            path.write_text(
                "# test\n"
                "A276275 ,2,2,0,4,2,4,6,6,10,12,16,22,\n"
                "A000045 ,0,1,1,2,3,5,8,13,\n",
                encoding="utf-8",
            )
            matches = osc.offline_matches(
                "candidate",
                (2, 2, 0, 4, 2, 4, 6, 6, 10, 12),
                path,
                max_shift=2,
                min_match=6,
                limit=5,
            )
            self.assertEqual(matches[0].oeis_number, "A276275")
            self.assertEqual(matches[0].matched_prefix_length, 10)


if __name__ == "__main__":
    unittest.main()
