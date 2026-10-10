"""
EGFR hotspot picker for the Anthropic x Adaptyv Challenge 1 (Colab-ready).

What it does
  1. Downloads human (P00533) and mouse (Q01279) EGFR from UniProt and aligns them.
  2. Downloads antibody-EGFR complex structures from the RCSB PDB.
  3. Finds EGFR residues that contact the antibody (the epitope), measures how
     exposed they are on free EGFR, and checks whether each one is conserved in mouse.
  4. Ranks conserved, exposed epitope residues as hotspot candidates and lists
     conserved acidic residues near the epitope (places where a histidine on the
     binder could form a pH-sensitive contact).

Run in Colab:
    !pip install -q biopython pandas requests
    then paste this file into a cell (or upload and `%run egfr_hotspots.py`).

Before you rely on it:
  - Check on rcsb.org that each ID in PDB_IDS really is EGFR + the antibody you want.
    Add a panitumumab-EGFR complex ID after verifying it (contacts shared by both
    antibodies are the strongest consensus epitope).
  - UniProt accession IDs are checked by the printed sequence lengths (~1210 aa).
"""
import requests
import pandas as pd
from Bio import Align
from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, Select, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

try:
    from Bio.Data.PDBData import residue_max_acc
    MAX_ACC = residue_max_acc["Sander"]
except Exception:  # older Biopython
    from Bio.PDB.DSSP import residue_max_acc
    MAX_ACC = residue_max_acc["Sander"]

# ---------------- settings ----------------
PDB_IDS = ["1YY9", "5SX4"]  # cetuximab and panitumumab complexes (verify IDs on rcsb.org)
HUMAN_ACC, MOUSE_ACC = "P00533", "Q01279"
CONTACT_CUTOFF = 4.5      # angstrom, antibody atom to EGFR atom
MIN_RSASA = 0.20          # relative exposure needed on free EGFR
ACID_RADIUS = 8.0         # angstrom, how near the epitope an Asp/Glu must be
TOP_N = 8                 # hotspots to report per structure
# ------------------------------------------


def uniprot_seq(acc):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30)
    r.raise_for_status()
    return "".join(r.text.splitlines()[1:])


def make_aligner(mode):
    a = Align.PairwiseAligner()
    a.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score, a.mode = -10, -0.5, mode
    return a


GLOBAL, LOCAL = make_aligner("global"), make_aligner("local")


def index_map(a_seq, b_seq, aligner):
    """0-based index in a_seq -> 0-based index in b_seq, for aligned positions."""
    aln = aligner.align(a_seq, b_seq)[0]
    m = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1):
            m[s1 + k] = s2 + k
    return m


def chain_info(chain):
    res = [r for r in chain if r.id[0] == " " and is_aa(r, standard=True)]
    return res, "".join(seq1(r.get_resname()) for r in res)


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


human, mouse = uniprot_seq(HUMAN_ACC), uniprot_seq(MOUSE_ACC)
print(f"Human EGFR: {len(human)} aa | Mouse EGFR: {len(mouse)} aa")
print(f"Global identity human vs mouse: "
      f"{sum(human[i] == mouse[j] for i, j in index_map(human, mouse, GLOBAL).items()) / len(human):.1%} of human residues")
h2m = index_map(human, mouse, GLOBAL)
blosum = Align.substitution_matrices.load("BLOSUM62")

all_rows = []
for pid in PDB_IDS:
    cif = requests.get(f"https://files.rcsb.org/download/{pid}.cif", timeout=60)
    cif.raise_for_status()
    open(f"{pid}.cif", "w").write(cif.text)
    model = MMCIFParser(QUIET=True).get_structure(pid, f"{pid}.cif")[0]

    # pick the EGFR chain by alignment score to human EGFR
    scores = {}
    for ch in model:
        res, s = chain_info(ch)
        if len(s) >= 50:
            scores[ch.id] = LOCAL.score(human, s)
    egfr_id = max(scores, key=scores.get)
    print(f"\n{pid}: chain alignment scores vs human EGFR: {scores} -> EGFR chain {egfr_id}")
    if scores[egfr_id] < 500:
        print("  WARNING: no chain looks like EGFR; check this PDB ID.")
        continue

    egfr = model[egfr_id]
    egfr_res, egfr_seq = chain_info(egfr)
    ab_atoms = [a for ch in model if ch.id != egfr_id and ch.id in scores and scores[ch.id] < 500
                for r in chain_info(ch)[0] for a in r]
    ns_ab = NeighborSearch(ab_atoms)

    # exposure on free EGFR (antibody excluded)
    ShrakeRupley().compute(egfr, level="R")
    p2h = index_map(egfr_seq, human, LOCAL)   # PDB chain index -> human UniProt index

    # epitope acidic-residue search
    acidic_near = set()
    for i, r in enumerate(egfr_res):
        if seq1(r.get_resname()) in "DE":
            if any(ns_ab.search(a.coord, ACID_RADIUS) for a in r):
                acidic_near.add(i)

    for i, r in enumerate(egfr_res):
        contacts = sum(len(ns_ab.search(a.coord, CONTACT_CUTOFF)) for a in r)
        if contacts == 0 and i not in acidic_near:
            continue
        aa = seq1(r.get_resname())
        h_idx = p2h.get(i)
        m_idx = h2m.get(h_idx) if h_idx is not None else None
        m_aa = mouse[m_idx] if m_idx is not None else "-"
        all_rows.append(dict(
            pdb=pid, chain=egfr_id, resnum=r.id[1], aa=aa,
            uniprot_pos=None if h_idx is None else h_idx + 1,
            mouse_aa=m_aa,
            conserved=(m_aa == aa),
            similar=(m_aa != "-" and blosum[aa][m_aa] > 0),
            rsasa=round(r.sasa / MAX_ACC[r.get_resname()], 2),
            contacts=contacts,
            near_acidic=(i in acidic_near),
        ))

    io = PDBIO()
    io.set_structure(model)
    io.save(f"{pid}_EGFR_only.pdb", ChainSelect(egfr_id))
    print(f"  saved {pid}_EGFR_only.pdb (chain {egfr_id}, {len(egfr_res)} residues)")

df = pd.DataFrame(all_rows)
df.to_csv("egfr_epitope_table.csv", index=False)

for pid, g in df.groupby("pdb"):
    epi = g[g.contacts > 0]
    good = epi[epi.conserved & (epi.rsasa >= MIN_RSASA)].sort_values("contacts", ascending=False)
    print(f"\n=== {pid}: {len(epi)} epitope residues, "
          f"{int(epi.conserved.sum())} conserved in mouse, {len(good)} conserved + exposed ===")
    print(good.head(TOP_N)[["resnum", "aa", "uniprot_pos", "mouse_aa", "rsasa", "contacts"]].to_string(index=False))
    print("Hotspot residue numbers (PDB numbering, for a hotspot setting):",
          ",".join(str(x) for x in sorted(good.head(TOP_N).resnum)))

    print("\nEpitope residues that DIFFER in mouse (avoid making these key contacts):")
    diff = epi[~epi.conserved]
    print(diff[["resnum", "aa", "mouse_aa", "similar"]].to_string(index=False) if len(diff) else "  none")

    acid = g[g.near_acidic & g.conserved & (g.rsasa >= MIN_RSASA) & g.aa.isin(list("DE"))]
    print(f"\nConserved, exposed Asp/Glu within {ACID_RADIUS} A of the antibody footprint "
          f"(candidate partners for a binder histidine):")
    print(acid[["resnum", "aa", "uniprot_pos", "rsasa"]].to_string(index=False) if len(acid) else "  none")

# ---- consensus across antibody complexes (matched by UniProt position) ----
epi_all = df[(df.contacts > 0) & df.uniprot_pos.notna()].copy()
n_struct = epi_all.pdb.nunique()
if n_struct > 1:
    epi_all["label"] = epi_all.pdb + ":" + epi_all.resnum.astype(str)
    cons = epi_all.groupby("uniprot_pos").agg(
        structures=("pdb", "nunique"), aa=("aa", "first"),
        mouse_aa=("mouse_aa", "first"), conserved=("conserved", "first"),
        mean_rsasa=("rsasa", "mean"), total_contacts=("contacts", "sum"),
        pdb_numbering=("label", ", ".join)).reset_index()
    shared = cons[cons.structures == n_struct]
    best = shared[shared.conserved & (shared.mean_rsasa >= MIN_RSASA)] \
        .sort_values("total_contacts", ascending=False)
    cons.to_csv("egfr_consensus_epitope.csv", index=False)
    print(f"\n=== Consensus: {len(shared)} residues contacted by all {n_struct} antibodies; "
          f"{len(best)} of them conserved in mouse and exposed ===")
    print(best[["uniprot_pos", "aa", "mouse_aa", "mean_rsasa", "total_contacts", "pdb_numbering"]]
          .to_string(index=False))
    print("Prefer these as hotspots; residues seen by only one antibody are in egfr_consensus_epitope.csv.")

print("\nSaved: egfr_epitope_table.csv and <PDB>_EGFR_only.pdb. "
      "Use the EGFR-only PDB as the design target and keep PDB numbering for hotspots.")
