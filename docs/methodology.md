> Design sequences are intentionally not included in this repository (they will be published after the competition results). Metrics below are from predicted structures only; nothing has been tested in the lab.

# EGFR Challenge 1: methodology details

## Target and hotspots
- Target: EGFR domain III, from the panitumumab complex (PDB 5SX4, chain M); the cetuximab complex (1YY9) was used to define the shared epitope.
- Epitope residues (4.5 A from either antibody) were ranked by contact count, exposure on free EGFR, and conservation between human (P00533) and mouse (Q01279) EGFR.
- Hotspots used for design: Q384, H409, K465 (PDB numbering). All three are conserved in mouse. The mouse-different residue S418 (Gly in mouse) was avoided as a key contact.

## Design
- BoltzGen (run on Tamarind Bio), Protein task, binder length 70-100, 100 designs requested (96 returned), budget 20, binding-site filter on. 34 of 96 passed BoltzGen filters.
- Finalists with atomic overlaps (<2.0 A) or antibody-variable-domain-like sequence motifs were removed, then the automatic novelty check (3/4 or higher) was applied after upload.

## Cross-species and control checks (Boltz-2)
- Shortlisted designs were refolded with Boltz-2 against the matching mouse EGFR domain III stretch (MSA for the target, none for the binder), 5 models each.
- Negative control: a shuffled copy of d32 (same composition, scrambled order) scored ipTM 0.32-0.73, ipSAE 0.01-0.36, pDockQ2 0.01-0.07 against mouse EGFR, versus ipTM about 0.95 and ipSAE about 0.92 for the real designs.
- Note: an attempt to fold three designs in one six-chain job was discarded because the binders hopped between EGFR copies.

## pH-dependence hypothesis
- PROPKA predicts H409 has a pKa of about 6.3 (roughly 38% protonated at pH 6.5 versus about 7% at pH 7.4, from Henderson-Hasselbalch). H409 is conserved in mouse and is the most-contacted hotspot.
- In d32 (E69) and d31 (D83), an acidic binder residue sits about 3.3-3.9 A from H409 in the predicted complexes. The idea is a charge pair that strengthens when H409 is protonated. d32_A62E adds a second Glu nearby.
- Limits: structure predictors do not model protonation, so pH selectivity is untested. Most of the other designs have no specific pH mechanism; they are included for human/mouse binding. Any selectivity is expected to be modest.

## Designs submitted
| Design | Length | BoltzGen ipTM / ipSAE / PAE | Additional evidence |
|---|---|---|---|
| egfr_d32 | 91 | 0.62 / 0.40 / 3.7 | Boltz-2 on human AND mouse EGFR: ipTM 0.95-0.96, ipSAE 0.93-0.94 (5/5 models, both species); E69 3.3-3.9 A from H409 |
| egfr_d27 | 90 | 0.55 / 0.23 / 5.0 | Boltz-2 on mouse: ipTM 0.88-0.91, ipSAE 0.76-0.82 (5/5 models); no contact with S418G; no acid near H409 |
| egfr_d31 | 91 | 0.58 / 0.34 / 4.1 | Boltz-2 on mouse: ipTM 0.66-0.89, ipSAE 0.33-0.84 (3/5 models strong, 2/5 weaker); D83 3.3-3.5 A from H409 in the strong models |
| egfr_d32_A62E | 91 | (parent d32: 0.62 / 0.40 / 3.7) | Boltz-2 on mouse: ipTM 0.95-0.96, ipSAE 0.92-0.93 (5/5 models); d32 plus a second Glu near H409 (within 4.5 A of H409 in 2/5 models) |
| egfr_d33 | 91 | 0.54 / 0.25 / 5.1 | BoltzGen filters and novelty check only (no Boltz-2 check) |
| egfr_d38 | 91 | 0.52 / 0.23 / 5.0 | BoltzGen filters and novelty check only (no Boltz-2 check) |
| egfr_d21 | 90 | 0.50 / 0.18 / 6.0 | BoltzGen filters and novelty check only (no Boltz-2 check) |
| egfr_d96 | 98 | 0.50 / 0.18 / 5.8 | BoltzGen filters and novelty check only (no Boltz-2 check) |
| egfr_d20 | 90 | 0.47 / 0.16 / 6.3 | BoltzGen filters and novelty check only (no Boltz-2 check) |
| egfr_d73 | 98 | 0.47 / 0.15 / 5.8 | BoltzGen filters and novelty check only (no Boltz-2 check) |
| egfr_d39 | 91 | 0.46 / 0.17 / 6.6 | BoltzGen filters and novelty check only (no Boltz-2 check) |
