# src/gw_assoc/analysis/odds.py

def compute_posterior_odds(gw_data, transient, **kwargs):
    """
    Minimal stub for posterior odds calculation.
    Replace with real math later.
    """
    return {
        "gw_file": getattr(gw_data, "get", lambda *_: None)("file"),
        "transient": {
            "ra": getattr(transient, "ra", None),
            "dec": getattr(transient, "dec", None),
            "z": getattr(transient, "z", None),
            "time": getattr(transient, "time", None),
        },
        "posterior_odds": None,
        "note": "compute_posterior_odds not implemented yet",
    }
