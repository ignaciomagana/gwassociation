def compute_posterior_odds(gw, transient, **kwargs):
    """
    Potato MVP:
      I_omega : placeholder sky overlap (constant)
      I_dl    : placeholder distance overlap (constant)
      I_t     : placeholder temporal overlap (constant)
      odds    : product of the three
    """
    I_omega = kwargs.get("I_omega", 0.01)  # pretend sky prob
    I_dl    = kwargs.get("I_dl",    0.50)  # pretend distance overlap
    I_t     = kwargs.get("I_t",     0.30)  # pretend time overlap
    odds    = I_omega * I_dl * I_t

    return {
        "gw_file": gw.get("file") if isinstance(gw, dict) else getattr(gw, "file", None),
        "transient": {
            "ra": getattr(transient, "ra", None),
            "dec": getattr(transient, "dec", None),
            "z": getattr(transient, "z", None),
            "time": getattr(transient, "time", None),
        },
        "I_omega": I_omega,
        "I_dl": I_dl,
        "I_t": I_t,
        "odds": odds,
        "note": "Potato calculation (stub)",
    }
