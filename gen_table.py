"""Generate one eRPA loss table (parametrized np for multi-core boxes).
Usage: python gen_table.py <outfile> <te_eV> <wi_Ha> [np]
Mirrors dief.__main__ grids (30 dens 1e19-1e29, 50 E 1e-3..1e2 MeV).
"""
import sys, dief
if __name__ == '__main__':
    fn = sys.argv[1]
    te = float(sys.argv[2])
    wi = float(sys.argv[3])
    npp = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    dief.tabdedx(fn, ts=[te], dr=[1e19, 1e29, 30], er=[1e-3, 1e2, 50],
                 wi=wi, np=npp)
