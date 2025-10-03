import matplotlib.pyplot as plt

def plot_skymap(gw, transient, out_file="skymap.png"):
    """
    Potato MVP: blank canvas with a point at (RA, Dec).
    Replace with healpy mollview later.
    """
    fig, ax = plt.subplots()
    ax.set_title("Skymap (potato MVP)")
    if getattr(transient, "ra", None) is not None and getattr(transient, "dec", None) is not None:
        ax.plot(transient.ra, transient.dec, "o", label="Transient")
        ax.set_xlabel("RA [deg]"); ax.set_ylabel("Dec [deg]"); ax.legend()
    fig.savefig(out_file, bbox_inches="tight")
    plt.close(fig)
