import os
import json
import pathlib
import click
import numpy as np

@click.command()
@click.option("--gw-file", type=click.Path(exists=True, dir_okay=False), required=True,
              help="Path to GW skymap file (FITS).")
@click.option("--ra", type=float, required=True, help="Transient RA [deg].")
@click.option("--dec", type=float, required=True, help="Transient Dec [deg].")
@click.option("--z", type=float, default=None, help="Transient redshift.")
@click.option("--z-err", type=float, default=None, help="Redshift uncertainty.")
@click.option("--time", "ttime", type=float, required=True, help="Transient time (GPS or MJD).")
@click.option("--gw-time", type=float, default=None, help="GW event time (GPS).")
@click.option("--model", type=click.Choice(['kilonova', 'grb', 'afterglow']), 
              default='kilonova', help="EM counterpart model.")
@click.option("--out", "outdir", type=click.Path(file_okay=False), default="out",
              help="Output directory for results.")
@click.option("--verbose", is_flag=True, help="Verbose output.")
def main(gw_file, ra, dec, z, z_err, ttime, gw_time, model, outdir, verbose):
    """Run GW-EM association analysis"""
    
    # Only import here to avoid issues if package not fully installed
    from gw_assoc import Association
    from gw_assoc.plots import plot_association_summary
    
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    
    # Set GW time to slightly before transient time if not provided
    if gw_time is None:
        gw_time = ttime - 86400  # Default to 1 day before transient
    
    # Create association object
    # Note: gw_time is passed separately and handled by Association class
    assoc = Association(gw_file, {
        'ra': ra, 
        'dec': dec, 
        'z': z, 
        'z_err': z_err,
        'time': ttime,
        'gw_time': gw_time  # Association class will handle this
    })
    
    # Compute odds with proper model
    results = assoc.compute_odds(em_model=model)
    
    if verbose:
        print("\n=== GW-EM Association Analysis Results ===")
        print(f"GW File: {gw_file}")
        print(f"Transient: RA={ra:.3f}°, Dec={dec:.3f}°")
        if z is not None:
            print(f"Redshift: z={z:.4f} ± {z_err:.4f}" if z_err else f"Redshift: z={z:.4f}")
        print(f"Time: {ttime:.2f} (transient), {gw_time:.2f} (GW)")
        print(f"EM Model: {model}")
        print(f"\nOverlap Integrals:")
        print(f"  Spatial (I_Ω):  {results['I_omega']:.3e}")
        print(f"  Distance (I_DL): {results['I_dl']:.3e}")  
        print(f"  Temporal (I_t):  {results['I_t']:.3e}")
        print(f"\nStatistics:")
        print(f"  Bayes Factor:    {results['bayes_factor']:.3e}")
        print(f"  Posterior Odds:  {results['posterior_odds']:.3e}")
        print(f"  Log₁₀ Odds:      {results['log_posterior_odds']:.2f}")
        print(f"  P(Associated):   {results['confidence']:.1%}")
        print(f"\nDecision: {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")
    else:
        print(f"P(Associated) = {results['confidence']:.1%}")
        print(f"Decision: {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")
    
    # Generate plots
    try:
        fig_path = out / "skymap.png"
        assoc.plot_skymap(str(fig_path))
        if verbose:
            print(f"\nSaved skymap: {fig_path}")
        
        # Additional plots if we have the plotting module
        try:
            summary_path = out / "association_summary.png"
            plot_association_summary(results, str(summary_path))
            if verbose:
                print(f"Saved summary: {summary_path}")
        except Exception as e:
            if verbose:
                print(f"Warning: Could not generate summary plot: {e}")
                
    except Exception as e:
        if verbose:
            print(f"Warning: Could not generate plots: {e}")
    
    # Save results to JSON
    with open(out / "results.json", "w") as f:
        # Convert numpy types to native Python types for JSON serialization
        json_results = {}
        for k, v in results.items():
            if isinstance(v, (np.float32, np.float64)):
                json_results[k] = float(v)
            elif isinstance(v, (np.int32, np.int64)):
                json_results[k] = int(v)
            elif isinstance(v, np.ndarray):
                json_results[k] = v.tolist()
            elif v == np.inf:
                json_results[k] = "inf"
            elif v == -np.inf:
                json_results[k] = "-inf"
            elif isinstance(v, dict):
                # Handle nested dictionaries
                json_results[k] = {
                    kk: float(vv) if isinstance(vv, (np.float32, np.float64)) else vv
                    for kk, vv in v.items()
                }
            else:
                json_results[k] = v
        
        json.dump(json_results, f, indent=2, default=str)
    
    if verbose:
        print(f"Saved results: {out/'results.json'}")
        print("\n=== Analysis Complete ===")

if __name__ == "__main__":
    main()