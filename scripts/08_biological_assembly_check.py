"""
Check whether residues of interest sit at a biological-assembly interface.

Why this exists
  Three times in this competition I designed against a surface that is buried in the
  real complex. A deposited PDB file is an asymmetric unit, which is NOT always the
  biological unit:
    - the deposited chains can be MORE than the biological unit (crystal neighbours), and
    - the true partner of a chain can be a SYMMETRY COPY that is not deposited at all.
  REMARK 350 states the author-determined biological unit and the transforms that build
  it. This script applies those transforms and reports, for each residue you name, the
  minimum distance to every other copy in the assembly. Anything under ~5 A is buried in
  the real complex and must not be used as a design hotspot.

Usage
  python 08_biological_assembly_check.py --entry 11QC --chain A --residues 189,192,380,381,382
  python 08_biological_assembly_check.py --entry 3S9L --chain A --residues 211,212,641,644,647

Needs: biopython, numpy, requests (and network access to files.rcsb.org).
"""
import argparse, os, re, sys
import numpy as np
import requests
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

CONTACT = 5.0          # angstrom; below this a residue is at an interface
NEAR = 8.0             # angstrom; below this it is on the rim


def fetch(entry, cache="."):
    path = os.path.join(cache, f"{entry}.pdb")
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{entry}.pdb", timeout=60)
        r.raise_for_status()
        open(path, "w").write(r.text)
    return path


def parse_remark350(path):
    """Return [{'chains': [...], 'transforms': [(3x3, 3)], 'note': str}] per biomolecule."""
    bios, cur = [], None
    for line in open(path):
        if not line.startswith("REMARK 350"):
            continue
        body = line[11:].strip()
        if body.startswith("BIOMOLECULE:"):
            if cur:
                bios.append(cur)
            cur = {"id": body.split(":")[1].strip(), "chains": [], "transforms": [], "note": ""}
        elif cur is None:
            continue
        elif "BIOLOGICAL UNIT:" in body:
            cur["note"] = body.split("BIOLOGICAL UNIT:")[1].strip()
        elif body.startswith("APPLY THE FOLLOWING TO CHAINS:"):
            cur["chains"] += [c.strip() for c in body.split(":")[1].split(",") if c.strip()]
        elif body.startswith("AND CHAINS:"):
            cur["chains"] += [c.strip() for c in body.split(":")[1].split(",") if c.strip()]
        elif re.match(r"BIOMT[123]\s", body):
            n = int(body[5])
            vals = [float(x) for x in body.split()[2:6]]
            if n == 1:
                cur["transforms"].append([vals])
            else:
                cur["transforms"][-1].append(vals)
    if cur:
        bios.append(cur)
    for b in bios:
        b["transforms"] = [(np.array([t[0][:3], t[1][:3], t[2][:3]]),
                            np.array([t[0][3], t[1][3], t[2][3]])) for t in b["transforms"]]
    return bios


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entry", required=True, help="PDB entry id, e.g. 11QC")
    ap.add_argument("--chain", required=True, help="author chain the residues belong to")
    ap.add_argument("--residues", required=True, help="comma-separated author residue numbers")
    ap.add_argument("--biomolecule", default="1", help="which BIOMOLECULE to build (default 1)")
    args = ap.parse_args()

    path = fetch(args.entry)
    model = PDBParser(QUIET=True).get_structure(args.entry, path)[0]
    deposited = sorted(ch.id for ch in model)
    bios = parse_remark350(path)
    bio = next((b for b in bios if b["id"] == args.biomolecule), None)
    if bio is None:
        sys.exit(f"No BIOMOLECULE {args.biomolecule} in REMARK 350 of {args.entry}")

    print(f"{args.entry}: deposited chains {deposited}")
    print(f"BIOMOLECULE {bio['id']} ({bio['note']}): chains {bio['chains']}, "
          f"{len(bio['transforms'])} transform(s)")
    extra = [c for c in deposited if c not in bio["chains"]]
    if extra:
        print(f"  NOTE: deposited chains {extra} are NOT part of this biological unit. "
              f"Contacts with them are crystal packing, not biology.")
    if len(bio["transforms"]) > 1:
        print(f"  NOTE: the assembly needs {len(bio['transforms'])} copies. A chain's true partner "
              f"may be a symmetry copy that is not in the deposited file.")

    # build every copy of the assembly except the query chain's identity copy
    others = {}
    for ci, (R, t) in enumerate(bio["transforms"]):
        for cid in bio["chains"]:
            if cid not in [c.id for c in model]:
                continue
            coords = np.array([a.coord for r in model[cid] if is_aa(r, standard=True) for a in r])
            if not len(coords):
                continue
            moved = coords @ R.T + t
            if ci == 0 and cid == args.chain:
                continue                       # the query chain itself
            others[f"{cid}_op{ci+1}"] = moved
    for cid in extra:                          # report crystal neighbours separately
        coords = np.array([a.coord for r in model[cid] if is_aa(r, standard=True) for a in r])
        if len(coords):
            others[f"{cid}_NOT_IN_ASSEMBLY"] = coords

    print(f"\n{'residue':10s}" + "".join(f"{k:>22s}" for k in others))
    verdicts = []
    for num in [int(x) for x in args.residues.split(",")]:
        try:
            res = model[args.chain][num]
        except KeyError:
            print(f"{num:<10d}  absent from chain {args.chain}")
            continue
        rc = np.array([a.coord for a in res])
        row, best_real = f"{seq1(res.get_resname())}{num:<9d}", 1e9
        for k, pts in others.items():
            d = float(np.min(np.linalg.norm(pts[None, :, :] - rc[:, None, :], axis=2)))
            row += f"{d:22.1f}"
            if "NOT_IN_ASSEMBLY" not in k:
                best_real = min(best_real, d)
        flag = ("BURIED IN ASSEMBLY" if best_real < CONTACT else
                "rim" if best_real < NEAR else "free surface")
        verdicts.append((num, seq1(res.get_resname()), best_real, flag))
        print(row + f"   -> {flag}")

    print("\nVerdict (use only 'free surface' residues as design hotspots):")
    for num, aa, d, flag in verdicts:
        print(f"  {aa}{num}: nearest assembly partner {d:.1f} A -> {flag}")


if __name__ == "__main__":
    main()
