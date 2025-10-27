import numpy as np
import healpy as hp
from typing import Dict, Optional, Tuple
try:
    import ligo.skymap.io
    import ligo.skymap.distance
    LIGO_SKYMAP_AVAILABLE = True
except ImportError:
    LIGO_SKYMAP_AVAILABLE = False

def load_gw_skymap(path: str) -> Dict:
    '''
    Load GW skymap from FITS file
    
    Returns dict with skymap data and metadata
    '''
    if LIGO_SKYMAP_AVAILABLE:
        try:
            skymap_data = ligo.skymap.io.read_sky_map(path)
            
            if len(skymap_data) == 4:  # 2D skymap
                return {
                    'file': path,
                    'kind': 'gw_skymap_2d',
                    'data': skymap_data[0],
                    'header': skymap_data[1],
                    'nside': hp.npix2nside(len(skymap_data[0])),
                    'is_3d': False
                }
            else:  # 3D skymap
                distances = ligo.skymap.distance.parameters_to_marginal_moments(
                    skymap_data[0], skymap_data[1]
                )
                return {
                    'file': path,
                    'kind': 'gw_skymap_3d',
                    'data': skymap_data,
                    'distances': distances,
                    'is_3d': True
                }
        except Exception as e:
            print(f"Error loading with ligo.skymap: {e}")
    
    # Fallback to simple healpy loading
    try:
        skymap, header = hp.read_map(path, h=True, verbose=False)
        return {
            'file': path,
            'kind': 'gw_skymap_2d',
            'data': skymap,
            'header': dict(header),
            'nside': hp.npix2nside(len(skymap)),
            'is_3d': False
        }
    except Exception as e:
        print(f"Error loading skymap: {e}")
        # Return minimal placeholder
        return {'file': path, 'kind': 'gw_skymap', 'data': None}

class GWEvent:
    '''Container for GW event data'''
    def __init__(self, skymap_path: str, event_time: float, event_name: str = None):
        self.skymap_path = skymap_path
        self.event_time = event_time  # GPS time
        self.event_name = event_name
        self.skymap = None
        self.distances = None
        self.nside = None
        self.is_3d = False
        
    def load_skymap(self):
        '''Load skymap data'''
        skymap_dict = load_gw_skymap(self.skymap_path)
        self.skymap = skymap_dict.get('data')
        self.distances = skymap_dict.get('distances')
        self.nside = skymap_dict.get('nside')
        self.is_3d = skymap_dict.get('is_3d', False)
        return skymap_dict