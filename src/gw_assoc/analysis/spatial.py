import numpy as np
import healpy as hp
from scipy import integrate
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy_healpix as ah

def skymap_overlap_integral(gw_skymap, ext_skymap=None,
                            ra=None, dec=None,
                            gw_nested=True, ext_nested=True):
    '''
    Sky map overlap integral between two sky maps.
    [Your existing implementation from spatial.py]
    '''
    # [Keep your existing implementation here]
    # Just showing the structure
    
    if ra is not None and dec is not None:
        # Point source case
        theta = np.radians(90 - dec)
        phi = np.radians(ra)
        
        if isinstance(gw_skymap, dict):
            skymap_data = gw_skymap.get('data')
            nside = gw_skymap.get('nside')
        else:
            skymap_data = gw_skymap
            nside = hp.npix2nside(len(skymap_data))
        
        if skymap_data is not None:
            ipix = hp.ang2pix(nside, theta, phi, nest=gw_nested)
            norm = np.sum(skymap_data)
            if norm > 0:
                return skymap_data[ipix] * len(skymap_data) / norm
    
    return 0.0

class SpatialOverlap:
    '''Calculate spatial overlap integral I_Ω between GW skymap and EM position'''
    
    @staticmethod
    def compute(gw_event, em_transient, search_radius: float = 0.0) -> float:
        '''
        Calculate spatial overlap integral I_Ω
        
        Parameters:
        -----------
        gw_event: GWEvent object with loaded skymap
        em_transient: Transient object with position
        search_radius: Error radius in degrees (default 0 for point source)
        
        Returns:
        --------
        I_omega: Spatial overlap probability
        '''
        if gw_event.skymap is None:
            gw_event.load_skymap()
        
        # Use the existing skymap_overlap_integral function
        return skymap_overlap_integral(
            {'data': gw_event.skymap, 'nside': gw_event.nside},
            ra=em_transient.ra,
            dec=em_transient.dec,
            gw_nested=True
        )