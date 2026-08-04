import unittest
import modular_recurrence_analyzer as mra


class AnalyzerTests(unittest.TestCase):
    def test_fibonacci_mod3(self):
        recurrence = mra.RecurrenceSpec((1, 1), 3, 'fibonacci')
        info = mra.analyze_seed(recurrence, (0, 1))
        self.assertEqual(info.preperiod, 0)
        self.assertEqual(info.period_word, (0, 1, 1, 2, 0, 2, 2, 1))
        self.assertEqual(info.output_period, 8)
        self.assertEqual(mra.find_global_antiperiod(recurrence), 4)

    def test_fibonacci_half_antiperiod(self):
        report = mra.analyze_word((0, 1, 1, 2, 0, 2, 2, 1), 3)
        self.assertTrue(report.half_antiperiodic)
        self.assertEqual(report.first_half, (0, 1, 1, 2))
        self.assertEqual(report.second_half, (0, 2, 2, 1))
        self.assertEqual(report.reversed_second_half, (1, 2, 2, 0))
        self.assertLess(report.mod3_even_frequency_max, 1e-9)

    def test_tribonacci_repository_example(self):
        recurrence = mra.RecurrenceSpec((1, 1, 1), 3, 'tribonacci')
        info = mra.analyze_seed(recurrence, (2, 1, 2))
        self.assertEqual(
            info.period_word,
            (2, 1, 2, 2, 2, 0, 1, 0, 1, 2, 0, 0, 2),
        )

    def test_tritetranacci_cycle_classes(self):
        recurrence = mra.RecurrenceSpec((1, 1, 1, 0), 3, 'tritetranacci')
        cycles = mra.enumerate_cycles(recurrence)
        words = {mra.word_to_text(word) for word in cycles}
        expected = {
            '0',
            '12',
            '011022',
            '00111201',
            '00222102',
            '01220211',
            '000101122122220012101021',
            '000202211211110021202012',
        }
        self.assertEqual(words, expected)

    def test_complement_pair_vs_self_antiperiod(self):
        self.assertTrue(mra.analyze_word((0, 1, 1, 0, 2, 2), 3).half_antiperiodic)
        a = (0, 0, 1, 1, 1, 2, 0, 1)
        b = mra.negate_word(a, 3)
        self.assertFalse(mra.analyze_word(a, 3).half_antiperiodic)
        self.assertEqual(b, (0, 0, 2, 2, 2, 1, 0, 2))


if __name__ == '__main__':
    unittest.main()
