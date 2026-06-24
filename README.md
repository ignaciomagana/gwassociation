# gwassociation: GW Association & Overlap Screening

`gwassociation` provides two command-line tools that share the Singer et al.
(2014) sky-map overlap integral:

- **`gwassociation odds`** — GW–EM association: posterior odds that an EM
  transient is the counterpart to a gravitational-wave event (sky, distance,
  and time overlap).
- **`gwassociation screen`** — GW–GW overlap screening: rank pairs of GW sky
  maps by spatial overlap to flag candidate gravitationally lensed pairs.

## Install

```console
pip install -e .            # core (odds)
pip install -e ".[screen]"  # adds screening (GraceDB client, ligo.skymap, healpy)
```

## `gwassociation odds` — GW–EM association

```console
gwassociation odds --gw-file SKYMAP.fits --ra RA --dec DEC \
  --z Z --time T --gw-time T_GW --out OUTDIR
```

Or from Python:

```python
from gwassociation import Association
results = Association("SKYMAP.fits", transient).compute_odds(em_model="kilonova")
```

## `gwassociation screen` — GW–GW overlap screening

**Run on new data** (query GraceDB, download sky maps, compute all pairs):

```console
gwassociation screen --time-window "START .. END" --workers N --out OUTDIR
```

**Re-use existing O4 overlaps** (skip query/download/compute; re-rank and
re-plot instantly from a previous run):

```console
gwassociation screen --reuse-overlaps OVERLAPS.json --skymap-dir SKYMAP_DIR --out OUTDIR
```

Outputs in `--out`: `overlaps.json`, `top_pairs.csv`, `overlap_histogram.png`,
`survival_function.png`, joint sky maps for the top pairs, and `summary.json`.
Add `--min-days-apart 1` to drop same-day duplicate detections of one trigger.
The Python entry point is `gwassociation.screening.run.run_screening`.

## Development

```console
pip install -e ".[dev]"
pytest -q tests
```

## Modules

- `gwassociation.association` — GW–EM odds API.
- `gwassociation.screening` — overlap-screening pipeline (`gracedb`, `download`,
  `overlap`, `pipeline`, `statistics`, `plots`, `run`).
- `gwassociation.io` / `analysis` / `stats` / `plotting` — loaders, overlap
  terms, priors, and plotting helpers.
