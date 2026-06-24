# gwassociation: GW Association & Overlap Screening

`gwassociation` is a lightweight Python package with two complementary tools:

- **GW–EM association** (`gwassociation odds`): evaluate whether an
  electromagnetic transient is the counterpart to a gravitational-wave (GW)
  event by combining sky-position, distance, and timing overlap terms into
  posterior odds.
- **GW–GW overlap screening** (`gwassociation screen`): screen *pairs* of GW
  sky maps for unusually high spatial overlap — a first-pass filter for
  strongly lensed event pairs (two detections of one source) and for
  prioritising multimessenger follow-up. This reproduces the analysis in
  *"Screening for Gravitationally Lensed Gravitational-Wave Events via Skymap
  Overlap"* (Pak & Magaña Hernández).

Both share the Singer et al. (2014) sky-map overlap integral.

## Install

```console
pip install -e .
```

Optional integrations are available as extras:

```console
pip install -e ".[healpy,ligo,hdf5,dev]"
```

The overlap-screening pipeline needs a heavier dependency set (GraceDB client,
`ligo.skymap`, `healpy`). Install it with the `screen` extra:

```console
pip install -e ".[screen]"
```

## Minimal Python API

```python
from gwassociation import Association

assoc = Association(
    "fits_files/S190425z_bayestar.fits.gz,0",
    {
        "name": "AT2024abc",
        "ra": 120.5,
        "dec": -30.2,
        "z": 0.05,
        "time": 1234567890.0,
        "gw_time": 1234567880.0,
    },
)
results = assoc.compute_odds(em_model="kilonova")
print(results["confidence"], results["posterior_odds"])
```

## Command line

The console script is a command group with two sub-commands.

### `gwassociation odds` — GW–EM association

```console
gwassociation odds \
  --gw-file fits_files/S190425z_bayestar.fits.gz,0 \
  --ra 120.5 \
  --dec -30.2 \
  --z 0.05 \
  --time 1234567890 \
  --gw-time 1234567880 \
  --out results
```

### `gwassociation screen` — GW–GW lensing-overlap screening

Query GraceDB, download the best sky map per event, compute the overlap
integral for every unique pair, and rank the highest-overlap pairs:

```console
gwassociation screen \
  --time-window "32 months ago .. now" \
  --far-threshold 2.3e-5 \
  --bbh-threshold 0.1 \
  --workers 6 \
  --top-n 30 \
  --min-days-apart 1 \
  --out screen_out
```

This writes to `--out`:

- `overlaps.json` — overlap integral for every pair (reusable; see below)
- `top_pairs.csv` — ranked top-N pairs with empirical p-values
- `overlap_histogram.png` — distribution of overlaps (full + high-overlap tail)
- `survival_function.png` — empirical survival function / p-value curve
- `joint_<e1>_<e2>.png` — joint `P1·P2` sky map for the top pair(s)
- `summary.json` — parameters, distribution stats, corrections, and rankings

`--min-days-apart 1` excludes same-day duplicate detections of a single
trigger. The all-pairs step is the expensive one; re-rank or re-plot without
recomputing by passing a saved overlaps file:

```console
gwassociation screen --reuse-overlaps screen_out/overlaps.json \
  --skymap-dir screen_out/skymaps --out screen_out
```

## Examples

The maintained example is intentionally small and uses the public API:

```console
python examples/minimal_script.py fits_files/S190425z_bayestar.fits.gz,0
```

## Development

```console
pip install -e ".[dev]"
pytest -q tests
```

The repository intentionally excludes generated plots, ad-hoc analysis outputs,
and notebooks so the package stays compact. The compact FITS fixture under `fits_files/` is retained because it is required for regression coverage.

## Core modules

- `gwassociation.association.Association`: high-level GW–EM API.
- `gwassociation.io`: sky-map and transient containers/loaders.
- `gwassociation.analysis`: spatial, line-of-sight distance, temporal, and odds calculations.
- `gwassociation.stats`: lower-level overlap and prior utilities.
- `gwassociation.plots` / `gwassociation.plotting`: optional plotting helpers.
- `gwassociation.screening`: GW–GW overlap screening pipeline
  (`gracedb`, `download`, `overlap`, `pipeline`, `statistics`, `plots`, `run`).
  `run.run_screening(...)` is the programmatic equivalent of `gwassociation screen`.

## License

See the project license file if one is added by the maintainers.
