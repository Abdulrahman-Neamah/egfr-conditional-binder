"""
Check whether predicted binders would clash with the REST of the EGFR extracellular region.

Designs were made against domain III only (residues 312-502). The lab measures binding to the full
extracellular region, so a binder that sits where domains I, II or IV are would not work there.

Usage:
    python check_full_receptor_clashes.py --receptor 1YY9_EGFR_only.pdb "complexes/**/*.pdb" --out clashes.csv

--receptor : full EGFR extracellular chain (e.g. chain A of PDB 1YY9, PDB numbering ~1-614)
complexes  : predicted 2-chain complexes (EGFR domain III target + binder), Boltz-2 or BoltzGen output

Method: superpose the complex's domain III onto the same residues of the full receptor (CA atoms matched by
PDB number, Kabsch), move the binder into the receptor frame, then measure its distance to receptor atoms
OUTSIDE residues 312-502 (those were visible during design). Heavy-atom distance < --clash A counts as a clash.
Limits: one receptor conformation only (the extended, antibody-bound one in 1YY9); predicted structures.
"""
import argparse
import collections
import glob

import numpy as np
import pandas as pd

FIRST = 312            # PDB number of the first residue of the 191-residue design target
LAST = 502


def load_pdb(path, chain=None):
    res = collections.OrderedDict()
    for l in open(path):
        if l.startswith("ATOM") and (chain is None or l[21] == chain):
            key = (l[21], int(l[22:26]), l[26])
            r = res.setdefault(key, {"atoms": {}})
            r["atoms"][l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return res


def kabsch(P, Q):
    """Rotation R and translation t that best map P onto Q (least squares)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def domain(n):
    return "I" if n <= 165 else "II" if n <= 310 else "IV" if n > LAST else "III/design target"


def analyse(path, receptor, rec_ca, target_len, clash):
    st = load_pdb(path)
    counts = collections.Counter(k[0] for k in st)
    tch = [c for c, n in counts.items() if n == target_len]
    if len(counts) != 2 or len(tch) != 1:
        return None
    tch = tch[0]
    bch = [c for c in counts if c != tch][0]
    trows = [v for k, v in st.items() if k[0] == tch]
    tnum = {FIRST + i: v["atoms"]["CA"] for i, v in enumerate(trows)}   # renumber to PDB numbering
    common = [n for n in tnum if n in rec_ca]
    if len(common) < 100:
        return dict(file=path, error=f"only {len(common)} matching CA atoms")
    R, t = kabsch(np.array([tnum[n] for n in common]), np.array([rec_ca[n] for n in common]))
    rms = float(np.sqrt(((np.array([tnum[n] for n in common]) @ R.T + t - np.array([rec_ca[n] for n in common])) ** 2).sum(1).mean()))
    binder = np.vstack([np.array(list(v["atoms"].values())) for k, v in st.items() if k[0] == bch]) @ R.T + t
    out_atoms, out_labels = [], []
    for (ch, num, ic), v in receptor.items():
        if num < FIRST or num > LAST:
            for name, xyz in v["atoms"].items():
                out_atoms.append(xyz)
                out_labels.append(num)
    out_atoms, out_labels = np.array(out_atoms), np.array(out_labels)
    D = np.sqrt(((binder[:, None, :] - out_atoms[None, :, :]) ** 2).sum(-1))
    nearest = int(out_labels[D.min(0).argmin()])
    return dict(file=path, align_CA=len(common), align_rmsd_A=round(rms, 2),
                min_dist_to_other_domains_A=round(float(D.min()), 1),
                clashing_binder_atoms=int((D.min(1) < clash).sum()),
                binder_atoms_within_4A=int((D.min(1) < 4.0).sum()),
                nearest_receptor_residue=nearest, nearest_domain=domain(nearest), error="")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pattern")
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--target-len", type=int, default=191)
    ap.add_argument("--clash", type=float, default=2.5)
    ap.add_argument("--out", default="clashes.csv")
    a = ap.parse_args()
    rec = load_pdb(a.receptor, chain="A")
    rec_ca = {k[1]: v["atoms"]["CA"] for k, v in rec.items() if "CA" in v["atoms"]}
    print(f"Receptor: {len(rec)} residues, numbers {min(rec_ca)}-{max(rec_ca)}")
    rows = [r for r in (analyse(f, rec, rec_ca, a.target_len, a.clash) for f in sorted(glob.glob(a.pattern, recursive=True))) if r]
    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    pd.set_option("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 70)
    print(df.to_string(index=False))
