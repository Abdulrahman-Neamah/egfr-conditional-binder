# Conditional EGFR binder design: Anthropic x Adaptyv Challenge 1

Computational workflow for the first challenge of the
[Anthropic x Adaptyv protein design competition](https://proteinbase.com/competitions/anthropic-adaptyv-2026):
design a de novo protein that binds **human and mouse EGFR** at the cetuximab/panitumumab site
and binds at **pH 6.5 but not pH 7.4**.

**Status:** designs submitted; lab results pending. **Design sequences are not included here** and will be
added after the results are published. Everything below is based on *predicted* structures and scores.
Nothing has been tested in the lab, and predictors do not model pH.

## Workflow

1. **Epitope and hotspots** (`scripts/01\_epitope\_hotspots.py`): find EGFR residues within 4.5 A of cetuximab
(PDB 1YY9) and panitumumab (PDB 5SX4), measure exposure on free EGFR, and check conservation between
human (UniProt P00533) and mouse (Q01279). Chosen hotspots: **Q384, H409, K465** (PDB numbering), all
conserved in mouse. The mouse-different residue S418 (Gly in mouse) was avoided as a key contact.
2. **Design** (`scripts/02\_boltzgen\_spec.py`, then the Tamarind Bio web interface): BoltzGen, protein-binder
task, 70-100 residues, binding-site filter on. 96 designs returned, 34 passed BoltzGen's filters.
3. **Cross-species and control checks** (`scripts/03\_make\_boltz2\_inputs.py`, Boltz-2 on Tamarind): one
complex per job against the matching mouse (and for the lead design, human) EGFR domain III,
five models each, plus a shuffled-sequence negative control.
4. **Structural analysis** (`scripts/04\_analyze\_complexes.py`): hotspot contacts, contacts with
human/mouse-different residues, overlaps, and the nearest binder Asp/Glu to His409.
5. **Full-receptor clash check** (`scripts/05\_check\_full\_receptor\_clashes.py`): designs were made against domain III only, but the lab
measures the whole extracellular region. Each predicted binder was moved onto domain III of full-length sEGFR (PDB 1YY9, which is in
the tethered conformation) and checked against domains I, II and IV.
6. **Novelty check** (competition platform, automatic): designs needed a score of 3/4 or higher.

## Key results (predicted)

|Check|Result|
|-|-|
|Lead design, Boltz-2 on human and on mouse EGFR|ipTM 0.95-0.96, ipSAE 0.93-0.94 in all 10 models; binding mode overlaps 88% between species|
|Shuffled-sequence control (mouse EGFR)|ipTM 0.32-0.73, ipSAE 0.01-0.36, pDockQ2 0.01-0.07|
|Three other designs on mouse EGFR|one consistent (ipTM 0.88-0.90), one strong in 3 of 5 models (0.66-0.89), one strong (0.94-0.95) but it failed the novelty check|
|H409 pKa (PROPKA)|about 6.3, so about 38% protonated at pH 6.5 versus about 7% at pH 7.4|
|Full-receptor clash check (1YY9, tethered)|0 clashing atoms in 20 BoltzGen finalists and 30 Boltz-2 models; closest approach to domains I, II, IV 19.5-26.1 A|
|Novelty check|11 of 14 uploaded designs scored 3/4 and were submitted; 3 scored 2/4|

Full tables: [`results/`](results/). Methods and per-design evidence: [`docs/methodology.md`](docs/methodology.md).

### pH hypothesis, and its limits

His409 is the most-contacted hotspot and has a predicted pKa near 6.3. In two designs an acidic binder residue
sits about 3.3-3.9 A from His409 in the predicted complexes, a possible charge pair that strengthens when
His409 is protonated. Structure predictors do not model protonation, so this is untested, and one
histidine can give at most a modest pH effect. Most submitted designs have no specific pH mechanism.

## Challenge 2 (TNF-alpha): negative result

I also tried the second challenge (a conditional binder of the TNF-alpha trimer: binds at pH 7.4, releases at pH 6.0, human and mouse).
Five BoltzGen runs on Tamarind (trimer and dimer targets, wide and compact hotspot sets, and an unrestricted diagnostic) produced no
usable design (best ipTM 0.46, ipSAE 0.11 with a valid target), so nothing was submitted. One run that looked good was invalid:
the binders had filled the position of the omitted third protomer. See `results/tnf\_boltzgen\_runs.csv`, `docs/lessons.md` (items 10-14)
and `scripts/06\_tnf\_epitope\_analysis.py` (epitope mapping and surface-patch ranking). I include it because the failure modes are instructive.

## Reproducing

```bash
pip install -r requirements.txt
python scripts/01\_epitope\_hotspots.py                  # downloads 1YY9, 5SX4, UniProt sequences
python scripts/02\_boltzgen\_spec.py                     # converts hotspots to BoltzGen numbering
cp designs\_template.csv designs.csv                    # add your own sequences (git-ignored)
python scripts/03\_make\_boltz2\_inputs.py designs.csv    # writes one FASTA per complex
python scripts/04\_analyze\_complexes.py "your\_results/\*\*/\*.pdb"
```

Scripts need internet for RCSB/UniProt. The design and Boltz-2 steps were run through
[Tamarind Bio](https://www.tamarind.bio) (free plan: 10 jobs/month). They have **not** been tested end to end
by anyone other than the author; expect small fixes (see `docs/lessons.md` for problems I hit).

## What did not work

See [`docs/lessons.md`](docs/lessons.md): BindCraft was too slow on a free Colab T4, JAX compatibility errors,
a six-chain Boltz-2 job that gave meaningless scores, three designs failing the novelty check, and more.

## Notes

* Tools: [BoltzGen](https://github.com/HannesStark/boltzgen), [Boltz-2](https://github.com/jwohlwend/boltz),
[PROPKA](https://github.com/jensengroup/propka), Biopython. Each has its own license; please check them
before reuse.
* The workflow and scripts were developed with help from Claude (Anthropic).
* Competition terms (publication license, effect on patentability) apply to the designs, not to this code.

## License

Code: MIT (see `LICENSE`).

