### Install

pip install gw-assoc
from gw_assoc import Association### 4. (Optional) Tag the release in git

From `gwPackage`:

git tag v0.1.0
git push origin v0.1.0### 5. For future changes

- Bump version (e.g. `0.1.1`) in your config
- Rebuild and re-upload:

python -m build
python -m twine upload dist/*If you tell me the exact output of `twine upload` (success message), I can confirm everything looks consistent with PyPI’s expectations.
