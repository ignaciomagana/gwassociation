Utilities to ingest GW volume localizations and evaluate association with arbitrary transients.
```bash
conda create -n gwPackage python=3.11 -y
conda activate gwPackage
pip install -e .
pytest -q


rm -rf build dist ./*.egg-info(N) src/*.egg-info(N)
python -m pip install -U pip setuptools wheel
python -m pip install -e .
python -m examples.minimal_script
