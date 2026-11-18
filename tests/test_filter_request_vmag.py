import pandas as pd
import pytest

from astroq.scripts.filter_request_vmag import apply_vmag_max


@pytest.mark.parametrize(
    "vmag_max,expected_inactive",
    [
        (8.0, [False, True, False]),
        (10.0, [False, False, False]),
        (5.0, [True, True, False]),
    ],
)
def test_apply_vmag_max(vmag_max, expected_inactive):
    df = pd.DataFrame(
        {
            "unique_id": ["bright", "faint", "none"],
            "Vmag": [6.0, 9.0, "none"],
            "inactive": [False, False, False],
        }
    )
    out, n_faint = apply_vmag_max(df, vmag_max)
    assert n_faint == sum(1 for v in [6.0, 9.0] if v > vmag_max)
    assert out["inactive"].tolist() == expected_inactive
