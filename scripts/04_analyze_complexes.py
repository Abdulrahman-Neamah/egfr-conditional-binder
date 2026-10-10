"""
Structural checks on predicted EGFR-binder complexes (Boltz-2 or BoltzGen output PDB files).

Usage:
    python 04_analyze_complexes.py "results/**/*.pdb" --out complex_summary.csv

For every 2-chain PDB (EGFR domain III + one binder) it reports:
  - binder chain, number of EGFR residues within 4.5 A of the binder
  - which hotspots (default Q384, H409, K465) are contacted
  - contacts with residues that differ between human and mouse EGFR (and the risky ones,
    BLOSUM62 <= 0; in practice only S418G mattered)
  - nearest binder Asp/Glu to the ring nitrogens of His409 (candidate pH-sensitive charge pair)
  - smallest binder-target heavy-atom distance (values below ~2.2 A hint at overlaps)

Target chain = the chain with --target-len residues (191 for 5SX4 domain III). Target residues
are renumbered to PDB numbering (first residue = 312), so Boltz-2 (1..191) and BoltzGen outputs
can be compared directly. Only numpy and pandas are needed. Predicted structures only.
"""
import argparse
import collections
import glob

import numpy as np
import pandas as pd

FIRST_RESIDUE = 312
HOTSPOTS = {384: "Q384", 409: "H409", 465: "K465"}
HIS = 409
# human -> mouse differences in the domain III stretch: PDB number -> (human, mouse, BLOSUM62)
DIFFS = {324: ("S", "T", 1), 328: ("D", "N", 1), 337: ("N", "Y", -2), 340: ("S", "A", 1),
         353: ("R", "K", 2), 359: ("H", "R", 0), 366: ("Q", "R", 1), 369: ("D", "E", 2),
         388: ("E", "D", 2), 390: ("R", "W", -3), 418: ("S", "G", 0), 420: ("D", "N", 1),
         443: ("K", "R", 2), 460: ("S", "P", -1), 461: ("G", "N", 0), 467: ("I", "M", 1),
         468: ("S", "N", 1), 471: ("G", "A", 0), 473: ("N", "K", 0), 474: ("S", "D", 0),
         478: ("T", "V", 0), 479: ("G", "N", 0), 480: ("Q", "H", 0), 483: ("H", "N", 1),
         484: ("A", "P", -1), 488: ("P", "S", -1), 502: ("H", "C", -3)}
RISKY = {n for n, (_, _, s) in DIFFS.items() if s <= 0}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def load(path):
    res = collections.OrderedDict()
    for line in open(path):
        if line.startswith("ATOM"):
            key = (line[21], int(line[22:26]), line[26])
            r = res.setdefault(key, {"name": line[17:20].strip(), "atoms": {}})
            r["atoms"][line[12:16].strip()] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return res


def dmat(a, b):
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))


def coords(res):
    return np.array(list(res["atoms"].values()))


def analyse(path, target_len, cutoff):
    st = load(path)
    counts = collections.Counter(k[0] for k in st)
    tch = [c for c, n in counts.items() if n == target_len]
    if len(counts) != 2 or len(tch) != 1:
        return None
    tch = tch[0]
    bch = [c for c in counts if c != tch][0]
    trows = [(k, v) for k, v in st.items() if k[0] == tch]
    target = {FIRST_RESIDUE + i: v for i, (_, v) in enumerate(trows)}
    binder = [(k[1], v) for k, v in st.items() if k[0] == bch]
    bxyz = np.vstack([coords(v) for _, v in binder])

    pairs = {n: int((dmat(coords(v), bxyz) <= cutoff).sum()) for n, v in target.items()}
    contact = {n for n, c in pairs.items() if c}
    mind = min(dmat(coords(v), bxyz).min() for v in target.values())
    # nearest binder Asp/Glu to His ring nitrogens
    his_n = np.array([target[HIS]["atoms"][a] for a in ("ND1", "NE2")])
    best = (None, None)
    for num, v in binder:
        if v["name"] in ("ASP", "GLU"):
            names = ("OD1", "OD2") if v["name"] == "ASP" else ("OE1", "OE2")
            d = dmat(np.array([v["atoms"][a] for a in names if a in v["atoms"]]), his_n).min()
            if best[0] is None or d < best[0]:
                best = (d, f"{ONE[v['name']]}{num}")
    return dict(
        file=path, target_chain=tch, binder_chain=bch, binder_len=len(binder),
        contact_residues=len(contact),
        hotspots_touched=",".join(lab for n, lab in HOTSPOTS.items() if n in contact) or "none",
        mouse_diff_contacts=",".join(f"{DIFFS[n][0]}{n}{DIFFS[n][1]}" for n in sorted(contact & set(DIFFS))) or "-",
        risky_mouse_diff_contacts=",".join(f"{DIFFS[n][0]}{n}{DIFFS[n][1]}" for n in sorted(contact & RISKY)) or "-",
        nearest_acid_to_his409=best[1] or "-",
        acid_to_his409_A=round(best[0], 1) if best[0] is not None else None,
        min_heavy_atom_distance_A=round(float(mind), 1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pattern", help='glob for PDB files, e.g. "results/**/*.pdb"')
    ap.add_argument("--target-len", type=int, default=191)
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--out", default="complex_summary.csv")
    a = ap.parse_args()
    files = sorted(glob.glob(a.pattern, recursive=True))
    rows = [r for r in (analyse(f, a.target_len, a.cutoff) for f in files) if r]
    if not rows:
        raise SystemExit(f"No usable 2-chain PDB files matched {a.pattern!r}")
    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    pd.set_option("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 60)
    print(df.drop(columns=["target_chain"]).to_string(index=False))
    print(f"\nWrote {a.out} ({len(df)} complexes)")
