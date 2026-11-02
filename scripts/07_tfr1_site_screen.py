"""
Screen every surface patch of human TfR1 for usability as a design epitope.

A patch is usable only if ALL of these hold:
  1. enough exposed apolar area to grip                  (from the patch map)
  2. conserved in mouse, for the cross-reactivity goal   (from the patch map)
  3. NOT buried by the TfR1 homodimer partner            (11QC, author-determined dimer)
  4. NOT inside the transferrin footprint                (3S9L, biological assembly 1)

Criteria 3 and 4 are the ones that are easy to get wrong, because both structures must be
expanded to their biological assemblies first - see 08_biological_assembly_check.py.
In 3S9L the deposited chains B and D are NOT part of assembly 1, so contacts with them are
crystal packing; the real partner of chain A is a symmetry copy built from REMARK 350.

Input : results/tfr1_patch_map.csv (centre, domain, apolar_area, apolar_frac, conserved_frac)
Output: results/tfr1_site_screen.csv

Usage : python 07_tfr1_site_screen.py
Needs : biopython, numpy, pandas, requests (network access to files.rcsb.org).
"""
import os
import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

DIMER_ENTRY, TF_ENTRY = "11QC", "3S9L"
TF_SYM = (np.diag([-1.0, 1.0, -1.0]), np.array([231.676, 0.0, 0.0]))  # 3S9L REMARK 350 BIOMT op 2
CONTACT, NEAR, MIN_CONS = 4.5, 8.0, 0.75
APICAL, HELICAL = range(189, 384), range(607, 761)


def fetch(entry):
    path = f"{entry}.pdb"
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{entry}.pdb", timeout=60)
        r.raise_for_status()
        open(path, "w").write(r.text)
    return path


def chain_coords(model, cid, op=None):
    c = np.array([a.coord for r in model[cid] if is_aa(r, standard=True) for a in r])
    return c if op is None else c @ op[0].T + op[1]


def main():
    dim = PDBParser(QUIET=True).get_structure("d", fetch(DIMER_ENTRY))[0]
    tfm = PDBParser(QUIET=True).get_structure("t", fetch(TF_ENTRY))[0]

    partner = chain_coords(dim, "B")                     # real dimer partner (11QC is DIMERIC as deposited)
    transferrin = np.vstack([chain_coords(tfm, "C"),     # chain A's own transferrin
                             chain_coords(tfm, "C", TF_SYM)])  # and the symmetry copy
    patches = pd.read_csv("results/tfr1_patch_map.csv").set_index("center")

    def mind(res, pts):
        rc = np.array([a.coord for a in res])
        return float(np.min(np.linalg.norm(pts[None, :, :] - rc[:, None, :], axis=2)))

    def domain(n):
        return "apical" if n in APICAL else ("helical" if n in HELICAL else "protease-like")

    rows = []
    for label in patches.index:
        num = int(label[1:])
        try:
            rd, rt = dim["A"][num], tfm["A"][num]
        except KeyError:
            continue
        aa = seq1(rd.get_resname())
        if aa != label[0]:
            raise SystemExit(f"numbering mismatch at {label}: structure has {aa}{num}")
        p = patches.loc[label]
        d_dim, d_tf = mind(rd, partner), mind(rt, transferrin)
        rows.append(dict(
            residue=label, domain=domain(num),
            apolar_area=p.apolar_area, apolar_frac=p.apolar_frac,
            conserved_frac_mouse=p.conserved_frac,
            dist_to_dimer_partner_A=round(d_dim, 2),
            dist_to_transferrin_A=round(d_tf, 2),
            dimer_interface=d_dim < CONTACT,
            transferrin_contact=d_tf < CONTACT,
            usable=bool(d_dim >= NEAR and d_tf >= NEAR and p.conserved_frac >= MIN_CONS)))

    df = pd.DataFrame(rows).sort_values(["usable", "apolar_area"], ascending=[False, False])
    os.makedirs("results", exist_ok=True)
    df.to_csv("results/tfr1_site_screen.csv", index=False)
    print(f"{len(df)} patches | {int(df.dimer_interface.sum())} at the dimer interface | "
          f"{int(df.transferrin_contact.sum())} in the transferrin footprint | "
          f"{int(df.usable.sum())} usable")
    print("\nTop usable patches by exposed apolar area:")
    cols = ["residue", "domain", "apolar_area", "apolar_frac", "conserved_frac_mouse",
            "dist_to_dimer_partner_A", "dist_to_transferrin_A"]
    print(df[df.usable].head(12)[cols].to_string(index=False))
    print("\nHighest-scoring patches that had to be REJECTED:")
    print(df[~df.usable].head(8)[cols + ["dimer_interface", "transferrin_contact"]].to_string(index=False))


if __name__ == "__main__":
    main()
