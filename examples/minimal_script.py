from gw_assoc import Association

GW_FILE = "Bilby.offline1.multiorder.fits.gz"  # adjust path if needed
TRANSIENT = {"ra": 255.72, "dec": 21.52, "z": 0.03136, "time": 1234567890}

def main():
    assoc = Association(GW_FILE, TRANSIENT)
    results = assoc.compute_odds()   # potato numbers
    print("=== gw_assoc MVP results ===")
    for k, v in results.items():
        print(f"{k:>10s}: {v}")
    assoc.plot_skymap("skymap_demo.png")
    print("Saved skymap_demo.png")

if __name__ == "__main__":
    main()
