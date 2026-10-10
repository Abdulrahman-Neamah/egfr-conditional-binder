"""
Build Boltz-2 input files to test binders against MOUSE and HUMAN EGFR domain III.

Usage (free Colab CPU runtime is enough; needs internet for UniProt):
    python 03_make_boltz2_inputs.py designs.csv

designs.csv must have columns: name,sequence   (keep it local; it is git-ignored)

For every design it writes two small FASTA files (chain A = EGFR domain III, chain B = binder):
    boltz2_inputs/<name>_mouse.fasta and boltz2_inputs/<name>_human.fasta
Submit ONE complex per Boltz-2 job. Putting several complexes in one job lets binders hop
between EGFR copies and gives meaningless scores (see docs/lessons.md).

Human target = the 191-residue EGFR domain III chain of PDB 5SX4 (PDB numbering 312-502).
Mouse target = the matching stretch of mouse EGFR (UniProt Q01279), found by sliding the human
chain along it without gaps (requires >= 85% identity, otherwise the script stops).
"""
import csv
import os
import sys

import requests

HUMAN = ("VCNGIGIGEFKDSLSIDATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRT"
         "DLHAFENLEIIRGRTKQHGQFSLAVVSLDITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRG"
         "ENSCKATGQVCHALCSPEGCWGPEPRDCVSH")
FIRST_RESIDUE = 312          # PDB (author) number of the first residue in HUMAN
HOTSPOTS = (384, 409, 465)   # residues the designs were built around

if len(sys.argv) != 2:
    raise SystemExit(__doc__)
with open(sys.argv[1], newline="") as fh:
    designs = {r["name"]: r["sequence"].strip().upper() for r in csv.DictReader(fh)}
if not designs:
    raise SystemExit("No designs found in the CSV (need columns: name,sequence).")

# ---- fetch mouse EGFR and find the matching stretch ----
r = requests.get("https://rest.uniprot.org/uniprotkb/Q01279.fasta", timeout=30)
r.raise_for_status()
mouse_full = "".join(r.text.splitlines()[1:])
n = len(HUMAN)
best_ident, best_off = max(
    (sum(a == b for a, b in zip(HUMAN, mouse_full[off:off + n])), off)
    for off in range(len(mouse_full) - n + 1))
MOUSE = mouse_full[best_off:best_off + n]
print(f"Mouse EGFR {len(mouse_full)} aa; best ungapped match {best_ident}/{n} ({best_ident / n:.0%})")
if best_ident / n < 0.85:
    raise SystemExit("Match too weak; there may be an insertion/deletion. Stop and check.")

diffs = [(i + FIRST_RESIDUE, a, b) for i, (a, b) in enumerate(zip(HUMAN, MOUSE)) if a != b]
print(f"{len(diffs)} human->mouse differences: " + ", ".join(f"{a}{num}{b}" for num, a, b in diffs))
for h in HOTSPOTS:
    i = h - FIRST_RESIDUE
    print(f"  hotspot {HUMAN[i]}{h}: mouse {MOUSE[i]} ({'conserved' if HUMAN[i] == MOUSE[i] else 'DIFFERENT'})")

os.makedirs("boltz2_inputs", exist_ok=True)
for name, seq in designs.items():
    for species, target in (("mouse", MOUSE), ("human", HUMAN)):
        with open(f"boltz2_inputs/{name}_{species}.fasta", "w") as fh:
            fh.write(f">A\n{target}\n>B\n{seq}\n")
print(f"Wrote {2 * len(designs)} FASTA files in boltz2_inputs/")
