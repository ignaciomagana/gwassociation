from .io import load_gw_skymap, GWEvent
from .io.transient import Transient
from .analysis import compute_posterior_odds
from .plotting.skymap import plot_skymap

# --- Add this near the top of association.py ---
from typing import Dict, Any

def _normalize_candidate_info(info: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize a candidate dict to the Transient(...) signature.
    Accepts keys like:
      - 'ra','dec' (deg)  ✅ target
      - 'ra_deg','dec_deg' (deg)
      - 'ra_rad','dec_rad' (radians)
      - 'coord' = astropy.coordinates.SkyCoord
    Also passes through z, z_err, time if present.
    """
    out = dict(info)
    out.pop("gw_time", None)  # not a Transient arg

    # Already good?
    if "ra" in out and "dec" in out:
        return out

    # ra/dec in degrees with different names
    if "ra_deg" in out and "dec_deg" in out:
        out["ra"] = float(out.pop("ra_deg"))
        out["dec"] = float(out.pop("dec_deg"))
        return out

    # radians -> degrees
    if "ra_rad" in out and "dec_rad" in out:
        import numpy as np
        out["ra"] = float(np.rad2deg(out.pop("ra_rad")))
        out["dec"] = float(np.rad2deg(out.pop("dec_rad")))
        return out

    # SkyCoord
    if "coord" in out:
        try:
            from astropy.coordinates import SkyCoord
            c = out.pop("coord")
            if isinstance(c, SkyCoord):
                out["ra"] = float(c.ra.deg)
                out["dec"] = float(c.dec.deg)
                return out
        except Exception:
            pass

    raise ValueError(
        "Candidate is missing required sky coordinates. "
        "Provide ('ra','dec') in degrees, or ('ra_deg','dec_deg'), "
        "or ('ra_rad','dec_rad') in radians, or 'coord'=SkyCoord."
    )

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
            candidate_info_clean = _normalize_candidate_info(candidate_info)

            candidate = Transient(**candidate_info_clean)

            original_transient = self.transient
            self.transient = candidate

            result = self.compute_odds()

            rankings.append({
                'candidate': candidate,
                'odds': result['posterior_odds'],
                'probability': result['confidence'],
                'log_odds': result['log_posterior_odds'],
                'results': result
            })

            self.transient = original_transient
        
        # Sort by odds (highest first)
        rankings.sort(key=lambda x: x['odds'], reverse=True)
        
        return rankings