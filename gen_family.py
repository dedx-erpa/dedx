"""Generate the Stage-2 cold gap-table family (te=0.025) over a wi grid.
Each table is generated in a fresh subprocess (clean isolation).
Usage: python gen_family.py [np]
"""
import subprocess, sys, os, time
GRID = [("0.0","tct_g0000"),("0.1","tct_g0100"),("0.2","tct_g0200"),
        ("0.4","tct_g0400"),("0.8","tct_g0800")]
TE = "0.025"
NP = sys.argv[1] if len(sys.argv) > 1 else "27"
if __name__ == '__main__':
    for wi, fn in GRID:
        t0 = time.time()
        print("=== %s.dat  wi=%s Ha  (%.1f eV) ==="%(fn, wi, float(wi)*27.2114), flush=True)
        r = subprocess.run([sys.executable, "gen_table.py", fn+".dat", TE, wi, NP],
                           cwd=os.path.dirname(os.path.abspath(__file__)))
        print("   %s rc=%d in %.0fs"%(fn, r.returncode, time.time()-t0), flush=True)
    print("FAMILY_DONE", flush=True)
