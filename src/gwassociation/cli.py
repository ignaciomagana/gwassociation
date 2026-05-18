from pathlib import Path
import json
import click
import numpy as np

from gwassociation.analysis.spatial import skymap_overlap_integral
from gwassociation.io.skymap import load_gw_skymap


def _jsonify(obj):
    if isinstance(obj, dict):
        return {k: _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, tuple):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


@click.group(invoke_without_command=True)
@click.pass_context
def main(ctx):
    """GWAssociation CLI tools."""
    if ctx.invoked_subcommand is None:
        ctx.invoke(analyze)


@main.command(name="analyze")
@click.option("--gw-file", type=click.Path(exists=True, dir_okay=False), required=True,
              help="Path to primary GW skymap file (FITS).")
@click.option("--secondary-skymap", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Optional secondary skymap (e.g., EM localization) for coincidence analysis.")
@click.option("--ra", type=float, default=None, help="Transient RA [deg] (required without --secondary-skymap).")
@click.option("--dec", type=float, default=None, help="Transient Dec [deg] (required without --secondary-skymap).")
@click.option("--z", type=float, default=None, help="Transient redshift.")
@click.option("--z-err", type=float, default=None, help="Redshift uncertainty.")
@click.option("--time", "ttime", type=float, default=None,
              help="Transient time (GPS or MJD). Required without --secondary-skymap.")
@click.option("--secondary-time", type=float, default=None,
              help="Event time for the secondary skymap (GPS).")
@click.option("--gw-time", type=float, default=None, help="GW event time (GPS).")
@click.option("--model", type=click.Choice(['kilonova', 'grb', 'afterglow']),
              default='kilonova', help="EM counterpart model.")
@click.option("--out", "outdir", type=click.Path(file_okay=False), default="out",
              help="Output directory for results.")
@click.option("--verbose", is_flag=True, help="Verbose output.")
def analyze(gw_file, secondary_skymap, ra, dec, z, z_err, ttime,
            secondary_time, gw_time, model, outdir, verbose):
    from gwassociation import Association
    from gwassociation.plots import plot_association_summary

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    if secondary_skymap is None:
        missing = [name for name, value in (("RA", ra), ("Dec", dec), ("time", ttime)) if value is None]
        if missing:
            raise click.UsageError(f"{', '.join(missing)} required unless --secondary-skymap is provided.")

    if gw_time is None:
        gw_time = ttime - 86400 if ttime is not None else 0.0

    transient_payload = None
    if ra is not None and dec is not None and ttime is not None:
        transient_payload = {'ra': ra, 'dec': dec, 'z': z, 'z_err': z_err, 'time': ttime, 'gw_time': gw_time}

    assoc = Association(gw_file, transient_payload, secondary_skymap=secondary_skymap, secondary_event_time=secondary_time)
    results = assoc.compute_odds(em_model=model)

    if verbose:
        click.echo(f"P(Associated): {results['confidence']:.1%}")
    else:
        click.echo(f"P(Associated) = {results['confidence']:.1%}")
        click.echo(f"Decision: {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")

    try:
        if transient_payload:
            assoc.plot_skymap(str(out / "skymap.png"))
        plot_association_summary(results, str(out / "association_summary.png"))
    except Exception:
        pass

    (out / "results.json").write_text(json.dumps(_jsonify(results), indent=2, default=str))


@main.command(name="pairwise-lensing-associations")
@click.option("--skymap-dir", type=click.Path(exists=True, file_okay=False), required=True,
              help="Directory containing GW FITS skymaps.")
@click.option("--pattern", default="*.fits*", show_default=True, help="Glob pattern for skymap files.")
@click.option("--top-k", type=int, default=10, show_default=True, help="Number of strongest pairs to print.")
@click.option("--out", "output_json", type=click.Path(dir_okay=False), default="pairwise_associations.json",
              show_default=True, help="Output JSON file path.")
@click.option("--download-gracedb", is_flag=True,
              help="Attempt GraceDB download for missing superevent FITS before computing associations.")
@click.option("--gracedb-ids", default="", help="Comma-separated superevent IDs used with --download-gracedb.")
def pairwise_lensing_associations(skymap_dir, pattern, top_k, output_json, download_gracedb, gracedb_ids):
    skymap_path = Path(skymap_dir)
    skymap_path.mkdir(parents=True, exist_ok=True)

    if download_gracedb:
        try:
            from ligo.gracedb.rest import GraceDb
        except Exception as exc:
            raise click.UsageError(f"GraceDB download requested but ligo.gracedb is unavailable: {exc}")
        client = GraceDb()
        ids = [item.strip() for item in gracedb_ids.split(",") if item.strip()]
        if not ids:
            raise click.UsageError("--gracedb-ids is required when --download-gracedb is set.")
        for superevent_id in ids:
            target = skymap_path / f"{superevent_id}.fits"
            if target.exists():
                continue
            fetched = False
            for remote_name in (
                "bayestar.multiorder.fits",
                "Bilby.offline0.multiorder.fits",
                "Bilby.multiorder.fits",
                "Bilby.multiorder.fits,0",
            ):
                try:
                    data = client.files(superevent_id, remote_name).data
                    target.write_bytes(data)
                    fetched = True
                    break
                except Exception:
                    continue
            if not fetched:
                click.echo(f"warning: unable to fetch skymap for {superevent_id}")

    files = sorted([p for p in skymap_path.glob(pattern) if p.is_file()])
    if len(files) < 2:
        raise click.UsageError("Need at least two skymap files to compute pairwise associations.")

    maps = {p.stem.split(".")[0]: load_gw_skymap(str(p)) for p in files}
    names = sorted(maps)
    overlaps = {}
    for i, name_i in enumerate(names):
        for name_j in names[i + 1:]:
            overlap = skymap_overlap_integral(maps[name_i], ext_skymap=maps[name_j])
            overlaps[f"{name_i},{name_j}"] = float(overlap)

    Path(output_json).write_text(json.dumps(overlaps, indent=2))
    click.echo(f"Wrote {len(overlaps)} pairwise overlaps to {output_json}")

    ranked = sorted(overlaps.items(), key=lambda kv: kv[1], reverse=True)[:max(top_k, 0)]
    for pair, value in ranked:
        click.echo(f"{pair}: {value:.6e}")


if __name__ == "__main__":
    main()
