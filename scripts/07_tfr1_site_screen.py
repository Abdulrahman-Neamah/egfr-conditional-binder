"""
Screen every exposed surface residue of human TfR1 for usability as a design hotspot.

A hotspot is usable only if ALL of these hold:
  1. the hotspot residue is itself solvent-exposed      (RSA >= 0.25)
  2. the hotspot residue is itself conserved in mouse   (strict identity)
  3. no residue of its surface patch is buried by the TfR1 homodimer partner (>= 8 A)
  4. no residue of its surface patch lies in the transferrin footprint       (>= 8 A)

Criterion 1 is the one this project learned the hard way. An earlier version of this
screen ranked patches by the apolar area of the whole neighbourhood and never checked the
exposure of the residue being named as a hotspot. Two of the five residues handed to the
generator for the "Site P" runs, K189 (RSA 0.14) and V192 (RSA 0.15), are mostly buried in
their own protomer, so no binder could engage them - a plausible part of why confident
designs kept landing somewhere other than the requested site. See docs/lessons_challenge4.md
items 33 and 34.

Criteria 3 and 4 aggregate over the whole patch, not the centre residue, because a binder
covers the neighbourhood rather than one side chain. Both require the BIOLOGICAL ASSEMBLY,
not the deposited chains - see 08_biological_assembly_check.py. Patch mean conservation is
reported but NOT gated: requiring it punishes a good conserved hotspot for having variable
fringes, and (worse) can pass a hotspot that differs in mouse itself because its neighbours
do not.

Output: results/tfr1_site_screen.csv
Usage : python 07_tfr1_site_screen.py
Needs : biopython, numpy, pandas, scipy, requests (network access to RCSB and UniProt).
"""
import os
import numpy as np
import pandas as pd
import requests
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

TF_ENTRY = "3S9L"                                                    # Tf-TfR1 complex
TF_SYM = (np.diag([-1.0, 1.0, -1.0]), np.array([231.676, 0.0, 0.0]))  # REMARK 350 BIOMT op 2
HUMAN, MOUSE = "P02786", "Q62351"
RSA_MIN, CLEAR, CONTACT, PATCH_R = 0.25, 8.0, 4.5, 10.0

# Tien et al. 2013, theoretical maximum accessible surface area
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197,
              L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)


def fetch_pdb(entry):
    p = f"{entry}.pdb"
    if not os.path.exists(p):
        r = requests.get(f"https://files.rcsb.org/download/{entry}.pdb", timeout=60)
        r.raise_for_status(); open(p, "w").write(r.text)
    return p


def uniprot(acc):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30)
    r.raise_for_status(); return "".join(r.text.splitlines()[1:])


def conservation():
    """human position (1-based) -> True if identical in mouse."""
    a = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1,
                              substitution_matrix=substitution_matrices.load("BLOSUM62"))
    h, m = uniprot(HUMAN), uniprot(MOUSE)
    aln = a.align(h, m)[0]
    cons = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1):
            cons[s1 + k + 1] = h[s1 + k] == m[s2 + k]
    ident = sum(cons.values())
    print(f"human {HUMAN} vs mouse {MOUSE}: {ident}/{len(h)} identical ({ident/len(h):.1%})")
    return cons


def main():
    model = PDBParser(QUIET=True).get_structure("t", fetch_pdb(TF_ENTRY))[0]
    cons = conservation()

    # --- exposure, computed on the ISOLATED protomer so interface residues still count ---
    import copy
    solo = copy.deepcopy(model["A"])
    for r in list(solo):
        if not is_aa(r, standard=True):
            solo.detach_child(r.id)
    ShrakeRupley(n_points=960).compute(solo, level="R")
    res = [r for r in solo if is_aa(r, standard=True)]
    rsa = {r.id[1]: r.sasa / MAXASA[seq1(r.get_resname())] for r in res}

    # --- biological assembly partners (REMARK 350 BIOMOLECULE 1: chains A + C, 2 ops) ---
    def coords(cid, op=None):
        c = np.array([a.coord for r in model[cid] if is_aa(r, standard=True) for a in r])
        return c if op is None else c @ op[0].T + op[1]
    partner = cKDTree(coords("A", TF_SYM))                       # the real dimer partner
    transferrin = cKDTree(np.vstack([coords("C"), coords("C", TF_SYM)]))

    def mind(r, tree):
        return float(min(tree.query(np.array([a.coord for a in r]))[0]))

    d_dim = {r.id[1]: mind(model["A"][r.id[1]], partner) for r in res}
    d_tf = {r.id[1]: mind(model["A"][r.id[1]], transferrin) for r in res}

    # --- patches: exposed residues whose CB (CA for Gly) is within PATCH_R of the centre ---
    centres = [r for r in res if rsa[r.id[1]] >= RSA_MIN]
    print(f"{len(res)} resolved residues | {len(centres)} exposed centres (RSA >= {RSA_MIN})")
    def anchor(r):
        return (r["CB"] if "CB" in r else r["CA"]).coord
    ids = [r.id[1] for r in centres]
    tree = cKDTree(np.array([anchor(r) for r in centres]))

    rows = []
    for r in centres:
        n = r.id[1]
        members = [ids[j] for j in tree.query_ball_point(anchor(r), PATCH_R)]
        rows.append(dict(
            residue=f"{seq1(r.get_resname())}{n}", resnum=n, rsa=round(rsa[n], 3),
            n_members=len(members),
            d_dimer_centre=round(d_dim[n], 2), d_tf_centre=round(d_tf[n], 2),
            d_dimer_patch_min=round(min(d_dim[m] for m in members), 2),
            d_tf_patch_min=round(min(d_tf[m] for m in members), 2),
            conserved_centre=bool(cons.get(n, False)),
            cons_patch_mean=round(np.mean([cons.get(m, False) for m in members]), 3),
            dimer_interface=min(d_dim[m] for m in members) < CONTACT,
            transferrin_contact=min(d_tf[m] for m in members) < CONTACT,
            usable=bool(rsa[n] >= RSA_MIN and cons.get(n, False)
                        and min(d_dim[m] for m in members) >= CLEAR
                        and min(d_tf[m] for m in members) >= CLEAR),
            members=" ".join(f"{seq1(model['A'][m].get_resname())}{m}" for m in sorted(members))))

    df = pd.DataFrame(rows).sort_values(["usable", "rsa"], ascending=[False, False])
    os.makedirs("results", exist_ok=True)
    df.to_csv("results/tfr1_site_screen.csv", index=False)
    print(f"\n{len(df)} patches | {int(df.dimer_interface.sum())} at the dimer interface | "
          f"{int(df.transferrin_contact.sum())} in the transferrin footprint | "
          f"{int(df.usable.sum())} usable")
    cols = ["residue", "rsa", "conserved_centre", "cons_patch_mean",
            "d_dimer_patch_min", "d_tf_patch_min"]
    print("\nUsable hotspots, most exposed first:")
    print(df[df.usable].head(15)[cols].to_string(index=False))

    print("\nThe Site P residues handed to the generator in the Challenge 4 runs:")
    for label in ["K189", "V192", "V380", "L381", "K382"]:
        n = int(label[1:])
        if n not in rsa:
            print(f"  {label}: not resolved"); continue
        hit = df[df.resnum == n]
        if hit.empty:
            print(f"  {label}: RSA {rsa[n]:.3f} -> REJECTED, below the {RSA_MIN} exposure cutoff; "
                  f"unreachable for a binder and should never have been a hotspot")
        else:
            x = hit.iloc[0]
            print(f"  {label}: RSA {x.rsa:.3f} | conserved {x.conserved_centre} | "
                  f"patch mean cons {x.cons_patch_mean} | usable {x.usable}")


if __name__ == "__main__":
    main()
