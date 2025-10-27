'''Galaxy density and rate calculations'''

import numpy as np
from astropy.cosmology import Planck15

def galaxy_density(z, galaxy_type='all'):
    '''
    Number density of galaxies at redshift z
    Returns: galaxies per Mpc^3
    '''
    # Simplified Schechter function parameters
    if galaxy_type == 'all':
        phi_star = 1e-2  # Mpc^-3
        alpha = -1.0
    elif galaxy_type == 'star_forming':
        phi_star = 5e-3
        alpha = -1.2
    else:  # elliptical
        phi_star = 3e-3
        alpha = -0.8
    
    # Evolution with redshift (simplified)
    evolution = (1 + z)**(-1)
    
    return phi_star * evolution

def expected_transients(volume, rate=1e-4):
    '''
    Expected number of transients in a volume
    
    Parameters:
    volume: Comoving volume in Mpc^3
    rate: Transient rate per Mpc^3 per year
    '''
    observation_time = 1.0  # year
    return volume * rate * observation_time

def false_alarm_probability(n_candidates, search_area, search_time):
    '''
    Probability of chance coincidence
    
    Parameters:
    n_candidates: Number of candidates found
    search_area: Sky area searched (square degrees)
    search_time: Time window (days)
    '''
    # All-sky transient rate (rough estimate)
    all_sky_rate = 100  # per day for magnitude < 20
    
    # Scale by search area (all sky = 41253 sq deg)
    area_fraction = search_area / 41253
    
    # Expected false alarms
    expected_false = all_sky_rate * area_fraction * search_time
    
    # Poisson probability
    from scipy.stats import poisson
    return 1 - poisson.cdf(n_candidates - 1, expected_false)