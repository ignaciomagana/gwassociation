# GW-EM Association Framework

A comprehensive Python framework for evaluating associations between gravitational wave (GW) events and electromagnetic (EM) transients using Bayesian statistics.

## Quick Info

The GW-EM Association Framework computes the probability that an electromagnetic transient is associated with a gravitational wave event using Bayesian statistics. It evaluates spatial, distance, and temporal overlap between GW skymaps and EM observations.

### Quick Start

**Command Line Interface:**
```bash
gw-assoc --gw-file skymap.fits --ra 120.5 --dec -30.0 --z 0.05 --time 1234567890
```

**Python API:**
```python
from gw_assoc import Association

# Create association analysis
assoc = Association("skymap.fits", {
    "ra": 120.5,      # Right ascension in degrees
    "dec": -30.0,     # Declination in degrees
    "z": 0.05,        # Redshift (optional)
    "z_err": 0.003,   # Redshift uncertainty (optional)
    "time": 1234567890,  # Detection time (GPS)
    "gw_time": 1234567889  # GW event time (GPS, optional)
})

# Compute association odds
results = assoc.compute_odds(
    em_model='kilonova',  # or 'grb', 'afterglow'
    prior_odds=1.0,
    chance_coincidence_rate=1e-4
)

# Get results
print(f"P(Associated) = {results['confidence']:.1%}")
print(f"Posterior Odds = {results['posterior_odds']:.3e}")
print(f"Bayes Factor = {results['bayes_factor']:.3e}")

# Rank multiple candidates
candidates = [
    {"ra": 120.5, "dec": -30.0, "z": 0.05, "time": 1234567890},
    {"ra": 121.0, "dec": -29.5, "z": 0.051, "time": 1234567900}
]
rankings = assoc.rank_candidates(candidates)
```

### Key Features

- **Bayesian Association Analysis**: Compute posterior odds for GW-EM associations
- **Multiple Overlap Integrals**: Spatial (I_Ω), distance (I_DL), and temporal (I_t)
- **EM Models**: Support for kilonova, GRB, and afterglow light curve models
- **Candidate Ranking**: Rank multiple EM candidates by association probability
- **Publication-Quality Plots**: Generate figures for papers and presentations
- **Command-Line Interface**: Easy-to-use CLI for quick analysis
- **3D Skymap Support**: Handles both 2D and 3D GW skymaps with distance information
- **Dual-Skymap Coincidence**: Compare two skymaps (GW vs EM) and compute radial overlap following *Coincident Detection Significance in Multimessenger Astronomy*

## Installation

### Prerequisites

- Python ≥ 3.7
- pip (Python package manager)

### Basic Installation

1. **Clone or download the repository:**
   ```bash
   cd gwPackage
   ```

2. **Install required dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Install the package in development mode:**
   ```bash
   pip install -e .
   ```

### Required Dependencies

The following packages are required and will be installed automatically:

- `numpy>=1.19.0` - Numerical computations
- `scipy>=1.5.0` - Scientific computing and statistical functions
- `astropy>=4.0` - Astronomy and astrophysics utilities
- `matplotlib>=3.3.0` - Plotting and visualization
- `click>=7.0` - Command-line interface framework
- `healpy>=1.14.0` - HEALPix sky map handling

### Recommended Dependencies

For enhanced functionality:

- `pandas>=1.1.0` - Data handling and manipulation
- `seaborn>=0.11.0` - Enhanced statistical plots

### Optional Dependencies

For advanced features (requires LIGO-Virgo-KAGRA access):

- `ligo.skymap>=1.0.0` - Advanced GW sky map handling
- `ligo-gracedb>=2.0.0` - GraceDB access for GW event data

### Verification

After installation, verify the package works correctly:

```bash
python test_gw_assoc.py
```

For minimal testing (without optional dependencies):

```bash
python test_gw_assoc.py --minimal
```

### Development Installation

For development with additional tools:

```bash
pip install -e ".[dev]"
```

This includes:
- `pytest>=6.0.0` - Testing framework
- `pytest-cov>=2.10.0` - Test coverage
- `black>=20.8b1` - Code formatting
- `flake8>=3.8.0` - Code linting

## How to Use with Real Inputs

This section provides step-by-step instructions for using the framework with real GW events and EM transients from actual observations.

### Getting Real GW Skymaps

#### From GraceDB (Public Events)

1. **Visit GraceDB:**
   - Public events: https://gracedb.ligo.org/superevents/public/O4/
   - Browse events or search for a specific event (e.g., `S250818k`)

2. **Download the Skymap:**
   - Click on an event to view details
   - Download the `bayestar.fits.gz` or `bilby.fits.gz` file
   - These are 3D skymaps with distance information (preferred)
   - Alternatively, 2D skymaps (`*.fits`) are also supported

3. **Example Event:**
   ```bash
   # Download S250818k skymap
   wget https://gracedb.ligo.org/api/superevents/S250818k/files/bayestar.fits.gz
   ```

#### From LIGO-Virgo-KAGRA Alerts

- **GCN Circulars**: Check GCN notices for public alerts
- **Public Data Releases**: Download from official data releases
- **API Access**: Use `ligo-gracedb` Python package (requires authentication for private events)

### Preparing EM Transient Data

#### Option 1: Python Dictionary

The simplest way is to create a dictionary with transient information:

```python
transient_info = {
    "name": "AT2024abc",           # Transient name (optional)
    "ra": 192.42625,                # Right ascension in degrees
    "dec": 34.82472,                # Declination in degrees
    "z": 0.438,                     # Redshift (optional but recommended)
    "z_err": 0.005,                 # Redshift uncertainty (optional)
    "time": 1242442967.447,         # Detection time (GPS seconds)
    "gw_time": 1242442965.0,        # GW event time (GPS, optional)
    "magnitude": 18.5,              # Apparent magnitude (optional)
    "filter_band": "r"              # Filter band (optional)
}
```

**Time Formats:**
- **GPS time**: Seconds since GPS epoch (recommended)
- **MJD**: Modified Julian Date (will be converted automatically)
- If `gw_time` is not provided, it will be estimated from the transient time

#### Option 2: JSON File

Create a JSON file with transient data:

```json
{
  "name": "AT2024abc",
  "ra": 192.42625,
  "dec": 34.82472,
  "z": 0.438,
  "z_err": 0.005,
  "time": 1242442967.447,
  "gw_time": 1242442965.0,
  "magnitude": 18.5,
  "filter_band": "r"
}
```

Load and use:
```python
import json
from gw_assoc import Association

# Load transient data
with open('transient.json', 'r') as f:
    transient_info = json.load(f)

# Create association
assoc = Association("S250818k_bayestar.fits.gz", transient_info)
results = assoc.compute_odds()
```

#### Option 3: CSV File

Create a CSV file with multiple transients:

```csv
name,ra,dec,z,z_err,time,magnitude,filter_band
AT2024abc,192.42625,34.82472,0.438,0.005,1242442967.447,18.5,r
AT2024def,192.43000,34.82000,0.440,0.006,1242442970.0,19.2,g
AT2024ghi,192.42000,34.83000,,,1242442968.0,17.8,i
```

Load and process:
```python
from gw_assoc import Association
from gw_assoc.ingest import ingest_transient_list

# Load transient list
candidates = ingest_transient_list('transients.csv')

# Create association (GW time from skymap or provided)
assoc = Association("S250818k_bayestar.fits.gz", {"gw_time": 1242442965.0})

# Rank all candidates
rankings = assoc.rank_candidates(candidates)

# Display results
for i, ranking in enumerate(rankings):
    print(f"{i+1}. {ranking['candidate'].name}: "
          f"P(Associated) = {ranking['probability']:.1%}")
```

### Complete Workflow Examples

#### Example 1: Single Transient Analysis

**Scenario:** You have a single EM transient candidate and want to check if it's associated with a GW event.

```python
from gw_assoc import Association

# 1. Prepare GW skymap and transient data
skymap_file = "S250818k_bayestar.fits.gz"
transient_info = {
    "name": "AT2024abc",
    "ra": 192.42625,      # From your observations
    "dec": 34.82472,      # From your observations
    "z": 0.438,           # From spectroscopy
    "z_err": 0.005,       # Measurement uncertainty
    "time": 1242442967.447,  # Detection time (GPS)
    "gw_time": 1242442965.0  # From GraceDB event page
}

# 2. Create association
assoc = Association(skymap_file, transient_info)

# 3. Compute odds (adjust parameters as needed)
results = assoc.compute_odds(
    em_model='kilonova',          # or 'grb', 'afterglow'
    prior_odds=1.0,               # BNS = 1.0, NSBH = 0.1, BBH = 0.01
    chance_coincidence_rate=1e-4, # Expected false alarm rate
    H0_uncertainty=7.0            # km/s/Mpc
)

# 4. Display results
print(f"Transient: {transient_info['name']}")
print(f"Position: RA={transient_info['ra']}°, Dec={transient_info['dec']}°")
print(f"Redshift: z={transient_info['z']} ± {transient_info['z_err']}")
print(f"\nOverlap Integrals:")
print(f"  Spatial (I_Ω):  {results['I_omega']:.3e}")
print(f"  Distance (I_DL): {results['I_dl']:.3e}")
print(f"  Temporal (I_t):  {results['I_t']:.3e}")
print(f"\nStatistics:")
print(f"  Bayes Factor:    {results['bayes_factor']:.3e}")
print(f"  Posterior Odds:  {results['posterior_odds']:.3e}")
print(f"  Log₁₀ Odds:      {results['log_posterior_odds']:.2f}")
print(f"  P(Associated):   {results['confidence']:.1%}")
print(f"\nDecision: {'✓ ASSOCIATED' if results['associated'] else '✗ NOT ASSOCIATED'}")

# 5. Generate plots
assoc.plot_skymap("association_skymap.png")
print(f"\nSaved skymap plot: association_skymap.png")
```

#### Example 2: Multiple Candidates Ranking

**Scenario:** You have multiple EM transient candidates and want to rank them by association probability.

```python
from gw_assoc import Association

# 1. Prepare GW skymap
skymap_file = "S250818k_bayestar.fits.gz"
gw_time = 1242442965.0  # From GraceDB

# 2. Prepare candidate list
candidates = [
    {
        "name": "AT2024abc",
        "ra": 192.42625,
        "dec": 34.82472,
        "z": 0.438,
        "z_err": 0.005,
        "time": 1242442967.447,
        "magnitude": 18.5
    },
    {
        "name": "AT2024def",
        "ra": 192.43000,
        "dec": 34.82000,
        "z": 0.440,
        "z_err": 0.006,
        "time": 1242442970.0,
        "magnitude": 19.2
    },
    {
        "name": "AT2024ghi",
        "ra": 192.42000,
        "dec": 34.83000,
        "z": None,  # No redshift available
        "time": 1242442968.0,
        "magnitude": 17.8
    }
]

# 3. Create association
assoc = Association(skymap_file, {"gw_time": gw_time})

# 4. Rank candidates
rankings = assoc.rank_candidates(candidates)

# 5. Display rankings
print("Candidate Rankings:")
print("-" * 80)
print(f"{'Rank':<6} {'Name':<12} {'RA':<10} {'Dec':<10} {'z':<8} {'P(Assoc)':<12} {'Decision'}")
print("-" * 80)

for i, ranking in enumerate(rankings):
    cand = ranking['candidate']
    prob = ranking['probability']
    decision = "✓ ASSOC" if prob > 0.5 else "✗ NOT ASSOC"
    z_str = f"{cand.z:.3f}" if cand.z else "N/A"
    
    print(f"{i+1:<6} {cand.name:<12} {cand.ra:<10.2f} {cand.dec:<10.2f} "
          f"{z_str:<8} {prob:<12.1%} {decision}")

# 6. Get top candidate for follow-up
top_candidate = rankings[0]
print(f"\nTop candidate: {top_candidate['candidate'].name}")
print(f"  P(Associated) = {top_candidate['probability']:.1%}")
print(f"  Posterior Odds = {top_candidate['odds']:.3e}")
```

#### Example 3: Command-Line Interface

**Scenario:** Quick analysis from the command line.

```bash
# Basic usage
gw-assoc \
  --gw-file S250818k_bayestar.fits.gz \
  --ra 192.42625 \
  --dec 34.82472 \
  --z 0.438 \
  --z-err 0.005 \
  --time 1242442967.447 \
  --gw-time 1242442965.0 \
  --model kilonova \
  --out results/ \
  --verbose

# With different EM model
gw-assoc \
  --gw-file S250818k_bayestar.fits.gz \
  --ra 192.42625 \
  --dec 34.82472 \
  --z 0.438 \
  --time 1242442967.447 \
  --model grb \
  --out results/

# Without redshift (spatial and temporal only)
gw-assoc \
  --gw-file S250818k_bayestar.fits.gz \
  --ra 192.42625 \
  --dec 34.82472 \
  --time 1242442967.447 \
  --out results/
```

#### Example 4: Batch Processing from File

**Scenario:** Process multiple candidates from a CSV or JSON file.

```python
from gw_assoc import Association
from gw_assoc.ingest import ingest_transient_list
import json

# 1. Load candidates from file
candidates = ingest_transient_list('my_transients.csv')  # or .json

# 2. Load GW event
skymap_file = "S250818k_bayestar.fits.gz"
gw_time = 1242442965.0

# 3. Create association
assoc = Association(skymap_file, {"gw_time": gw_time})

# 4. Process all candidates
rankings = assoc.rank_candidates(candidates)

# 5. Save results
results = []
for ranking in rankings:
    cand = ranking['candidate']
    results.append({
        "name": cand.name,
        "ra": cand.ra,
        "dec": cand.dec,
        "z": cand.z,
        "probability": ranking['probability'],
        "posterior_odds": ranking['odds'],
        "log_odds": ranking['log_odds'],
        "I_omega": ranking['results']['I_omega'],
        "I_dl": ranking['results']['I_dl'],
        "I_t": ranking['results']['I_t']
    })

# Save to JSON
with open('association_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"Processed {len(results)} candidates")
print(f"Results saved to: association_results.json")
```

#### Example 5: Skymap vs Skymap Coincidence (Radial Distance)

**Scenario:** Compare two full skymaps (e.g., GW and EM localization) and evaluate their coincident detection significance following *Coincident Detection Significance in Multimessenger Astronomy*.

```bash
gw-assoc \
  --gw-file S250818k_bayestar.fits.gz \
  --secondary-skymap em_localization.fits.gz \
  --secondary-time 1242443000.0 \
  --out results/ \
  --verbose
```

```python
from gw_assoc import Association

assoc = Association(
    "S250818k_bayestar.fits.gz",
    transient_info=None,
    secondary_skymap="em_localization.fits.gz",
    secondary_event_time=1242443000.0
)

results = assoc.compute_odds()
print(f"Spatial overlap (I_Ω) = {results['I_omega']:.3e}")
print(f"Radial overlap (I_DL) = {results['I_dl']:.3e}")
print(f"P(Associated) = {results['confidence']:.1%}")
```

In this mode:
- The framework loads both skymaps, validates that they share the same NSIDE, and computes angular overlap.
- Radial (distance) overlap is computed using the joint line-of-sight integral described in the paper, combining per-pixel distance posteriors.
- Temporal overlap defaults to 1, but you can pass `I_t` manually if the secondary skymap has an associated time window.

### Working with 3D Skymaps

3D skymaps (with distance information) provide more accurate distance overlap calculations:

```python
from gw_assoc import Association

# 3D skymaps automatically provide distance information
assoc = Association("S250818k_bayestar.fits.gz", {
    "ra": 192.42625,
    "dec": 34.82472,
    "z": 0.438,
    "z_err": 0.005,
    "time": 1242442967.447,
    "gw_time": 1242442965.0
})

# The framework automatically detects 3D skymaps
results = assoc.compute_odds()
print(f"Distance overlap (I_DL): {results['I_dl']:.3e}")

# For 2D skymaps, distance overlap will be 1.0 (no distance constraint)
```

### Time Format Conversions

If you have times in different formats, convert them:

```python
from astropy.time import Time

# Convert MJD to GPS
mjd_time = 58630.5
t = Time(mjd_time, format='mjd')
gps_time = t.gps  # GPS seconds

# Convert ISO format to GPS
iso_time = "2024-01-15T12:34:56"
t = Time(iso_time, format='iso')
gps_time = t.gps

# Convert GPS to MJD
gps_time = 1242442967.447
t = Time(gps_time, format='gps')
mjd_time = t.mjd

# Use in transient info
transient_info = {
    "ra": 192.42625,
    "dec": 34.82472,
    "time": gps_time,  # Use GPS time
    "gw_time": gw_gps_time
}
```

### Common Real-World Scenarios

#### Scenario 1: GW Follow-up Campaign

```python
# After receiving a GW alert and conducting observations
from gw_assoc import Association
from gw_assoc.ingest import ingest_transient_list

# 1. Download skymap from GraceDB (manual or automated)
skymap_file = "S250818k_bayestar.fits.gz"

# 2. Load candidates from your observation pipeline
candidates = ingest_transient_list('observed_candidates.csv')

# 3. Analyze associations
assoc = Association(skymap_file, {"gw_time": gw_time})
rankings = assoc.rank_candidates(candidates)

# 4. Select top candidates for spectroscopy
top_3 = rankings[:3]
for candidate in top_3:
    print(f"Priority target: {candidate['candidate'].name}")
    print(f"  RA: {candidate['candidate'].ra}°, Dec: {candidate['candidate'].dec}°")
    print(f"  P(Associated) = {candidate['probability']:.1%}")
```

#### Scenario 2: Retrospective Analysis

```python
# Analyzing historical events with known associations
from gw_assoc import Association

# GW170817-like analysis
assoc = Association("GW170817_skymap.fits.gz", {
    "name": "AT2017gfo",
    "ra": 197.45,
    "dec": -23.38,
    "z": 0.0098,
    "z_err": 0.0001,
    "time": 1187008882.4,  # Kilonova detection time
    "gw_time": 1187008882.43  # GW merger time
})

results = assoc.compute_odds(em_model='kilonova', prior_odds=1.0)
print(f"GW170817-AT2017gfo association:")
print(f"  P(Associated) = {results['confidence']:.1%}")
print(f"  Posterior Odds = {results['posterior_odds']:.3e}")
```

#### Scenario 3: Missing Data Handling

```python
# Handle cases where some data is missing
from gw_assoc import Association

# No redshift available
assoc1 = Association("skymap.fits.gz", {
    "ra": 192.42625,
    "dec": 34.82472,
    "time": 1242442967.447
    # No z: distance overlap will be 1.0
})
results1 = assoc1.compute_odds()
print(f"Spatial + temporal only: P = {results1['confidence']:.1%}")

# No time available
assoc2 = Association("skymap.fits.gz", {
    "ra": 192.42625,
    "dec": 34.82472,
    "z": 0.438
    # No time: temporal overlap will be 1.0
})
results2 = assoc2.compute_odds()
print(f"Spatial + distance only: P = {results2['confidence']:.1%}")

# Position only
assoc3 = Association("skymap.fits.gz", {
    "ra": 192.42625,
    "dec": 34.82472
    # Spatial overlap only
})
results3 = assoc3.compute_odds()
print(f"Spatial only: P = {results3['confidence']:.1%}")
```

### Tips for Real Data

1. **GW Event Time**: Always get the GW event time from GraceDB for accurate temporal calculations
2. **Redshift Quality**: Higher quality redshifts (smaller uncertainties) improve distance overlap calculations
3. **3D Skymaps**: Prefer 3D skymaps (bayestar.fits.gz, bilby.fits.gz) over 2D for better distance constraints
4. **EM Model Selection**: Choose the appropriate model:
   - `kilonova`: Optical/NIR transients (hours to days)
   - `grb`: Gamma-ray bursts (seconds)
   - `afterglow`: GRB afterglows (days to weeks)
5. **Prior Odds**: Adjust based on GW source type (BNS, NSBH, BBH)
6. **Multiple Candidates**: Always rank multiple candidates to identify the most likely association

### Troubleshooting

**Issue: Skymap file not found**
```python
# Check if file exists
from pathlib import Path
if not Path("skymap.fits.gz").exists():
    print("Download skymap from GraceDB first")
```

**Issue: Invalid time format**
```python
# Convert to GPS time
from astropy.time import Time
t = Time(your_time, format='your_format')
gps_time = t.gps
```

**Issue: Missing distance information**
```python
# Check if skymap is 3D
from gw_assoc.io.skymap import load_gw_skymap
skymap_data = load_gw_skymap("skymap.fits.gz")
if skymap_data.get('is_3d'):
    print("3D skymap detected - distance overlap will be calculated")
else:
    print("2D skymap - distance overlap = 1.0")
```

## Math

This framework implements a Bayesian statistical framework for evaluating GW-EM associations based on the formalism developed by Ashton et al. (2018, 2021). The core calculation computes the posterior odds that an EM transient is associated with a GW event.

### Bayesian Framework

The posterior odds for association are calculated as:

```
O_posterior = O_prior × BF
```

where `O_prior` is the prior odds ratio and `BF` is the Bayes factor.

### Bayes Factor

The Bayes factor compares the probability of the data under the hypothesis that the transient is associated with the GW event versus the hypothesis that it is not:

```
BF = P(data | associated) / P(data | not associated)
   = (I_Ω × I_DL × I_t) / P_chance
```

where:
- `I_Ω`: Spatial overlap integral
- `I_DL`: Distance (luminosity distance) overlap integral
- `I_t`: Temporal overlap integral
- `P_chance`: Chance coincidence probability

### Spatial Overlap Integral (I_Ω)

The spatial overlap integral measures the agreement between the GW sky localization and the EM transient position:

```
I_Ω = ∫ p_GW(Ω) × p_EM(Ω) / π_sky(Ω) dΩ
```

where:
- `p_GW(Ω)` is the GW sky localization probability density
- `p_EM(Ω)` is the EM transient position probability density (typically a point source or Gaussian)
- `π_sky(Ω)` is the prior sky probability (uniform: 1/(4π) steradians)

For a point source EM transient at position (RA, Dec), this simplifies to:

```
I_Ω = p_GW(RA, Dec) / (1/(4π))
```

where `p_GW(RA, Dec)` is the GW probability density at the transient position, normalized by the pixel area.

### Distance Overlap Integral (I_DL)

The distance overlap integral measures the agreement between the GW distance posterior and the EM transient distance (derived from redshift):

```
I_DL = ∫ p_GW(DL | Ω) × p_EM(DL) / π_DL(DL) dDL
```

where:
- `p_GW(DL | Ω)` is the GW distance posterior at the line of sight (for 3D skymaps)
- `p_EM(DL)` is the EM distance probability distribution (derived from redshift with uncertainties)
- `π_DL(DL)` is the distance prior (typically uniform in comoving volume: ∝ DL²)

For 3D skymaps, the line-of-sight distance density is:

```
p_LOS(DL | Ω) = DL² × Normal(DL; μ(Ω), σ(Ω)) × distnorm(Ω)
```

where:
- `μ(Ω)` and `σ(Ω)` are the per-pixel distance mean and standard deviation from the 3D skymap
- `distnorm(Ω)` is the per-pixel normalization factor
- The `DL²` factor accounts for the comoving volume prior

For EM transients, the distance is computed from redshift:

```
DL(z) = (c/H₀) × ∫₀^z dz' / E(z')
```

with uncertainties from:
- Redshift measurement error
- Peculiar velocity (~300 km/s)
- Hubble constant uncertainty (H₀ ≈ 73 ± 7 km/s/Mpc)

### Temporal Overlap Integral (I_t)

The temporal overlap integral accounts for the expected time delay between the GW merger and EM emission:

```
I_t = p(t_EM | t_GW, model)
```

where the probability depends on the EM counterpart model:

**Kilonova model:**
- Peak emission: ~1 day after merger
- Light curve: Log-normal rise with exponential decay
- Typical width: ~2 days

**GRB model:**
- Prompt emission: ~seconds after merger
- Gaussian temporal profile with σ ≈ 5 seconds

**Afterglow model:**
- Peak: ~1 day after merger
- Power-law decay: t^(-0.7) for t > t_peak

### Chance Coincidence Probability (P_chance)

The chance coincidence probability accounts for the expected rate of unrelated transients:

```
P_chance = R_EM × Δt × ΔΩ
```

where:
- `R_EM` is the all-sky transient rate (per day per square degree)
- `Δt` is the time window (typically days)
- `ΔΩ` is the searched sky area (square degrees)

### Posterior Probability

The posterior probability of association is:

```
P(associated | data) = O_posterior / (1 + O_posterior)
                    = 1 - 1/(1 + O_posterior)
```

A transient is considered associated if `O_posterior > 1` (or `P(associated) > 0.5`).

### Prior Odds

The prior odds depend on the GW source type:

- **BNS (Binary Neutron Star)**: `O_prior ≈ 1.0` (EM emission expected)
- **NSBH (Neutron Star-Black Hole)**: `O_prior ≈ 0.1` (EM emission possible)
- **BBH (Binary Black Hole)**: `O_prior ≈ 0.01` (EM emission unlikely)

## Citations

This framework implements statistical methods from the following papers and articles. Please cite these works when using this package in your research.

### Primary Methods

1. **Ashton, G., et al. (2018, 2021)** - Bayesian framework for GW-EM associations
   - Original formulation of the Bayesian association framework
   - Development of overlap integral formalism
   - Implementation of spatial, distance, and temporal overlap calculations

2. **Singer, L. P., & Price, L. R. (2016)** - Rapid sky localization with BAYESTAR
   - Rapid Bayesian sky localization algorithm
   - HEALPix skymap generation and handling
   - Reference: *Physical Review D*, 93, 024013
   - DOI: 10.1103/PhysRevD.93.024013

### Gravitational Wave Sky Localization

3. **LIGO-Virgo-KAGRA Collaboration (O4)** - Recent observational runs
   - O4 observing run papers on GW follow-up strategies
   - Sky localization improvements and 3D skymaps
   - Distance estimation methods

### Distance and Cosmology

4. **Planck Collaboration (2015)** - Cosmological parameters
   - Hubble constant and cosmological parameter measurements
   - Used for redshift-distance conversions
   - Reference: *Astronomy & Astrophysics*, 594, A13
   - DOI: 10.1051/0004-6361/201525830

### Electromagnetic Counterparts

5. **Kilonova models** - Various authors
   - Kilonova light curve models and temporal profiles
   - Expected time delays and light curve evolution
   - References in multi-messenger astronomy literature

6. **GRB afterglow models** - Various authors
   - Gamma-ray burst prompt and afterglow emission
   - Temporal profiles and light curve modeling
   - References in GRB and multi-messenger literature

### Software and Tools

7. **HEALPix** - Hierarchical Equal Area isoLatitude Pixelization
   - Sky map pixelization scheme
   - Górski, K. M., et al. (2005)
   - Reference: *The Astrophysical Journal*, 622, 759
   - DOI: 10.1086/427976

8. **Astropy Collaboration** - Astropy Project
   - Astronomy and astrophysics Python package
   - Cosmology calculations and coordinate transformations
   - Reference: *Astronomy & Astrophysics*, 558, A33
   - DOI: 10.1051/0004-6361/201322068

### Additional References

9. **Multi-messenger astronomy reviews** - Various authors
   - GW170817 and subsequent multi-messenger events
   - Association analysis methodologies
   - Follow-up strategies and best practices

10. **LIGO-Virgo-KAGRA Collaboration papers** - Various
    - Gravitational wave detection papers
    - Sky localization and distance estimation methods
    - Multi-messenger follow-up campaigns

### Citation Format

If you use this framework in your research, please cite:

```bibtex
@software{gw_assoc,
  author = {Iganacio Magana, Kaitlyn Pak},
  title = {GW-EM Association Framework},
  year = {2025},
  url = {https://github.com/ignaciomagana/gw-assoc},
  version = {0.2.0}
}
```

And please also cite the primary methodology papers (Ashton et al. 2018, 2021) and other relevant references from the list above.

## Additional Resources

### Documentation

See the `examples.py` script for detailed usage examples:
```bash
python examples.py
```

### Testing

Run the test suite to verify installation:
```bash
python test_gw_assoc.py
```

### Examples

Example scripts are available in the `examples/` directory:
- `minimal_script.py` - Basic usage example
- See `examples.py` for comprehensive examples

### Getting GW Skymaps

GW skymaps can be downloaded from:
- GraceDB: https://gracedb.ligo.org/superevents/
- LIGO-Virgo-KAGRA public alerts


## License

MIT License - see LICENSE file for details.

## Acknowledgments

This framework builds upon methods developed by the LIGO-Virgo-KAGRA collaboration and the broader multi-messenger astronomy community. We thank the developers of the underlying software packages (Astropy, HEALPix, ligo.skymap, etc.) for their invaluable tools.
