"""Tests for the [semester] mode switch and round sequencing."""

import unittest
from configparser import ConfigParser

from astroq.splan import MODE_SEQUENCES, _ROUND_SPECS, SemesterPlanner


def _planner_with_config(text):
    """Minimal planner stub: only ``config`` is needed by _resolve_mode."""
    sp = SemesterPlanner.__new__(SemesterPlanner)
    sp.config = ConfigParser()
    sp.config.read_string(text)
    return sp


class TestResolveMode(unittest.TestCase):
    def test_default_is_round1(self):
        sp = _planner_with_config("[semester]\n")
        mode, sequence = sp._resolve_mode()
        self.assertEqual(mode, "round1")
        self.assertEqual(sequence, ("Round1",))

    def test_full_runs_five_rounds(self):
        sp = _planner_with_config("[semester]\nmode = full\n")
        mode, sequence = sp._resolve_mode()
        self.assertEqual(mode, "full")
        self.assertEqual(
            sequence, ("Round1", "Round2", "Round3", "Round4", "UpcomingNight")
        )

    def test_unknown_mode_raises(self):
        sp = _planner_with_config("[semester]\nmode = bonus\n")
        with self.assertRaisesRegex(ValueError, "mode='bonus' invalid"):
            sp._resolve_mode()

    def test_legacy_keys_raise_with_migration_hint(self):
        for legacy in ("run_bonus_round", "run_upcoming_night_round"):
            sp = _planner_with_config(f"[semester]\n{legacy} = True\n")
            with self.assertRaisesRegex(ValueError, "mode = round1 | full"):
                sp._resolve_mode()

    def test_every_round_has_a_build_method(self):
        for sequence in MODE_SEQUENCES.values():
            for label in sequence:
                _, build_method = _ROUND_SPECS[label]
                self.assertTrue(
                    callable(getattr(SemesterPlanner, build_method)),
                    f"{label} -> {build_method} is not a SemesterPlanner method",
                )


if __name__ == "__main__":
    unittest.main()
