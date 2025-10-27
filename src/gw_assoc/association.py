from .io import load_gw_skymap, GWEvent
from .io.transient import Transient
from .analysis import compute_posterior_odds
from .plotting.skymap import plot_skymap

class Association:
    """High-level wrapper for evaluating GW-EM associations"""

    def __init__(self, gw_file: str, transient_info: dict):
        """
        Initialize Association
        
        Parameters
        ----------
        gw_file : str
            Path to GW skymap FITS file
        transient_info : dict
            Dictionary with transient information
            Should contain: ra, dec, z (optional), time (optional)
            Can also contain: gw_time for the GW event time
        """
        # Extract GW event time if provided, otherwise use a default
        gw_time = transient_info.pop('gw_time', None)
        if gw_time is None:
            # If no GW time provided, estimate from transient time
            if 'time' in transient_info and transient_info['time'] is not None:
                gw_time = transient_info['time'] - 86400  # Default: 1 day before transient
            else:
                gw_time = 0.0  # Fallback default
        
        # Create GW event
        self.gw = GWEvent(
            skymap_path=gw_file,
            event_time=gw_time
        )
        
        # Create transient (without gw_time which isn't a Transient parameter)
        self.transient = Transient(**transient_info)
        
        # Store GW time separately if needed
        self.gw_time = gw_time

    def compute_odds(self, **kwargs):
        """
        Compute association odds
        
        Parameters
        ----------
        **kwargs : dict
            Additional parameters for odds calculation:
            - em_model: 'kilonova', 'grb', or 'afterglow'
            - prior_odds: Prior odds ratio (default 1.0)
            - chance_coincidence_rate: Rate of chance coincidences
            - H0_uncertainty: Hubble constant uncertainty
        
        Returns
        -------
        dict
            Results dictionary with odds and probabilities
        """
        # Make sure skymap is loaded
        if self.gw.skymap is None:
            self.gw.load_skymap()
        
        return compute_posterior_odds(self.gw, self.transient, **kwargs)

    def plot_skymap(self, out_file: str = "skymap.png"):
        """Plot skymap with transient"""
        # Make sure skymap is loaded
        if self.gw.skymap is None:
            self.gw.load_skymap()
            
        plot_skymap(self.gw, self.transient, out_file)
        
    def rank_candidates(self, candidates_list):
        """
        Rank multiple candidates by association probability
        
        Parameters
        ----------
        candidates_list : list
            List of dictionaries with candidate information
        
        Returns
        -------
        list
            Sorted list of candidates with odds and probabilities
        """
        rankings = []
        
        for candidate_info in candidates_list:
            # Remove gw_time if present (not a Transient parameter)
            candidate_info_clean = candidate_info.copy()
            candidate_info_clean.pop('gw_time', None)
            
            # Create candidate transient
            candidate = Transient(**candidate_info_clean)
            
            # Temporarily set as the transient for this association
            original_transient = self.transient
            self.transient = candidate
            
            # Compute odds
            result = self.compute_odds()
            
            # Add to rankings
            rankings.append({
                'candidate': candidate,
                'odds': result['posterior_odds'],
                'probability': result['confidence'],
                'log_odds': result['log_posterior_odds'],
                'results': result  # Full results if needed
            })
            
            # Restore original transient
            self.transient = original_transient
        
        # Sort by odds (highest first)
        rankings.sort(key=lambda x: x['odds'], reverse=True)
        
        return rankings