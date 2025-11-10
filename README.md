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
  author = {Your Name},
  title = {GW-EM Association Framework},
  year = {2024},
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

## Contributing

Contributions are welcome! Please submit pull requests or open issues on GitHub.

## License

MIT License - see LICENSE file for details.

## Acknowledgments

This framework builds upon methods developed by the LIGO-Virgo-KAGRA collaboration and the broader multi-messenger astronomy community. We thank the developers of the underlying software packages (Astropy, HEALPix, ligo.skymap, etc.) for their invaluable tools.
