import json
import numpy as np
from astropy.table import Table
from click.testing import CliRunner

from gwassociation.cli import main


def _write_skymap(path, prob):
    t = Table()
    t["PROB"] = np.asarray(prob, dtype=float)
    t.write(path, overwrite=True)


def test_pairwise_lensing_associations_cli(tmp_path):
    skydir = tmp_path / "maps"
    skydir.mkdir()
    _write_skymap(skydir / "A.fits", np.ones(12) / 12)
    _write_skymap(skydir / "B.fits", np.ones(12) / 12)
    _write_skymap(skydir / "C.fits", np.array([1] + [0] * 11, dtype=float))

    out = tmp_path / "overlaps.json"
    runner = CliRunner()
    result = runner.invoke(main, [
        "pairwise-lensing-associations",
        "--skymap-dir", str(skydir),
        "--pattern", "*.fits",
        "--out", str(out),
        "--top-k", "2",
    ])

    assert result.exit_code == 0
    payload = json.loads(out.read_text())
    assert len(payload) == 3
    assert "A,B" in payload
    assert np.isclose(payload["A,B"], 1.0)
