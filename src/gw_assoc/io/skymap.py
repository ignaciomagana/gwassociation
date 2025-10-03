def load_gw_skymap(path: str):
    """
    Potato MVP loader: return a tiny dict.
    Swap with ligo.skymap/healpy reader later.
    """
    return {"file": path, "kind": "gw_skymap"}
