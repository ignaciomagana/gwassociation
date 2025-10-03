import os, json, pathlib, click
print("[gw-assoc] cli module imported")

import json
import pathlib
import click
from gw_assoc import Association

@click.command()
@click.option("--gw-file", type=click.Path(exists=True, dir_okay=False), required=True,
              help="Path to GW skymap file (FITS or placeholder).")
@click.option("--ra", type=float, required=True, help="Transient RA [deg].")
@click.option("--dec", type=float, required=True, help="Transient Dec [deg].")
@click.option("--z", type=float, default=None, help="Transient redshift.")
@click.option("--time", "ttime", type=float, required=True, help="Transient time (GPS or MJD).")
@click.option("--out", "outdir", type=click.Path(file_okay=False), default="out",
              help="Output directory for results.")
def main(gw_file, ra, dec, z, ttime, outdir):
    """Run a minimal association analysis and print results."""
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    assoc = Association(gw_file, {"ra": ra, "dec": dec, "z": z, "time": ttime})
    results = assoc.compute_odds()

    print("=== gw-assoc CLI results ===")
    for k, v in results.items():
        print(f"{k:>10s}: {v}")

    fig_path = out / "skymap.png"
    assoc.plot_skymap(str(fig_path))
    print(f"Saved {fig_path}")

    with open(out / "results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved {out/'results.json'}")

if __name__ == "__main__":
    main()
