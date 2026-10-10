"""
Write a BoltzGen design specification for the EGFR Challenge 1 target (Colab-ready).

Why this script exists
  BoltzGen's YAML uses mmCIF *label* numbering (starting at 1) and *label* chain IDs,
  NOT the author (PDB) numbering used by the epitope script. This script converts
  the hotspots, checks that each converted residue is the amino acid you expect,
  and writes egfr_spec.yaml.

Run in Colab:
    !pip install -q biopython requests
    then paste this into a cell and run it.
It needs internet access (downloads 5SX4.cif if it is not already in the folder).
"""
import os
import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# ---------------- settings ----------------
PDB_ID = "5SX4"
AUTH_CHAIN = "M"                        # EGFR chain found by the epitope script
HOTSPOTS_AUTH = {384: "GLN", 409: "HIS", 417: "VAL", 465: "LYS", 469: "ASN"}
BINDER_LENGTH = "60..100"               # random length range sampled per design
BINDER_ID = "Z"                         # any chain ID not used by the target
OUT_YAML = "egfr_spec.yaml"
# ------------------------------------------

cif_path = f"{PDB_ID}.cif"
if not os.path.exists(cif_path):
    r = requests.get(f"https://files.rcsb.org/download/{PDB_ID}.cif", timeout=60)
    r.raise_for_status()
    open(cif_path, "w").write(r.text)

d = MMCIF2Dict(cif_path)
cols = ["_atom_site.group_PDB", "_atom_site.auth_asym_id", "_atom_site.label_asym_id",
        "_atom_site.auth_seq_id", "_atom_site.label_seq_id", "_atom_site.label_comp_id",
        "_atom_site.pdbx_PDB_model_num"]
first_model = d["_atom_site.pdbx_PDB_model_num"][0]

label_chains, auth_to_label = set(), {}
for grp, ach, lch, aseq, lseq, comp, model in zip(*(d[c] for c in cols)):
    if grp == "ATOM" and ach == AUTH_CHAIN and model == first_model and lseq != ".":
        label_chains.add(lch)
        auth_to_label[int(aseq)] = (int(lseq), comp)

if len(label_chains) != 1:
    raise SystemExit(f"Expected one label chain for auth chain {AUTH_CHAIN}, found {label_chains}")
label_chain = label_chains.pop()
if label_chain == BINDER_ID:
    BINDER_ID = "Y"
print(f"Auth chain {AUTH_CHAIN} = label chain {label_chain}; "
      f"{len(auth_to_label)} residues modelled")

label_hotspots = []
for auth_num, expected in HOTSPOTS_AUTH.items():
    if auth_num not in auth_to_label:
        raise SystemExit(f"Residue {auth_num} not found in chain {AUTH_CHAIN}")
    lab, comp = auth_to_label[auth_num]
    ok = "OK" if comp == expected else f"MISMATCH (expected {expected})"
    print(f"  author {auth_num} {comp} -> label_seq_id {lab}   {ok}")
    if comp != expected:
        raise SystemExit("Residue type mismatch: stop and check the numbering before designing.")
    label_hotspots.append(lab)

hot_str = ",".join(str(x) for x in sorted(label_hotspots))
yaml_text = f"""entities:
  # Designed binder; length is sampled between the two numbers
  - protein:
      id: {BINDER_ID}
      sequence: {BINDER_LENGTH}

  # Target: EGFR domain III from {PDB_ID}, chain {AUTH_CHAIN} (label chain {label_chain})
  - file:
      path: {cif_path}
      include:
        - chain:
            id: {label_chain}
      binding_types:
        - chain:
            id: {label_chain}
            binding: {hot_str}
"""
open(OUT_YAML, "w").write(yaml_text)
print(f"\nWrote {OUT_YAML}:\n")
print(yaml_text)
print("Next: !boltzgen check " + OUT_YAML + "  and view the output in a structure viewer.")
