# What did not work (and what I changed)

A record of dead ends, so others can skip them.

1. **BindCraft on free Colab was too slow.** With a free T4, one trajectory took about 48 minutes, and
   2 trajectories finished in roughly 3.5 hours with 0 accepted designs. The authors suggest expecting
   hundreds of trajectories, so this was not feasible inside a one-week window. I switched to BoltzGen.
2. **Colab's newer JAX broke the old BindCraft/ColabDesign notebook.** `jax.lib.xla_bridge` no longer
   exists and `jnp.clip` no longer accepts `a_min`/`a_max`. Small monkey-patches got past both errors, but
   another compatibility error may follow, and they did not fix the speed problem.
3. **Residue numbering.** BoltzGen YAML specs use mmCIF *label* numbering (from 1) and label chain IDs, not
   the author (PDB) numbering used in papers. `scripts/02_boltzgen_spec.py` converts and checks each residue type.
   (For 5SX4 chain M, label numbering is author numbering minus 306; Tamarind's web form renumbered the chain 1-191 itself.)
4. **Running several complexes in one Boltz-2 job is invalid.** Three binders and three EGFR copies in one
   six-chain job let binders hop between EGFR copies: overall ipTM was 0.33-0.38 and meaningless per design.
   Run one complex (two chains, Copies = 1) per job.
5. **Boltz-2 scores need a control.** A shuffled-sequence copy of the best design still placed itself on the
   hotspots in most models (contact patterns alone prove nothing) but scored far lower (ipTM up to 0.73,
   ipSAE up to 0.36, pDockQ2 at most 0.07) than real designs (ipTM about 0.95). The control is an easy
   negative (an unfolded sequence), so it only shows that scores are not trivially high.
6. **Novelty filter.** Three of the top designs scored 2/4 on the competition's novelty check (below the
   required 3/4), including one of the strongest by Boltz-2. Three conservative surface-mutation variants
   of one of them (8, 14 and 20 substitutions at positions more than 6 A from the target) failed too, so
   surface changes do not seem to move the score. I do not know what the score measures.
7. **A second glutamate near His409 did not reliably reach it.** The A62E variant of the best design kept its
   predicted binding mode and scores, but the new Glu was within 4.5 A of His409 in only 2 of 5 models.
8. **Free-tier limits.** Tamarind's free plan allows 10 jobs per month, so jobs have to be budgeted
   (design batches, human/mouse checks, controls).
9. **Check the rest of the receptor.** Designs were generated against domain III alone, but binding is measured on the full extracellular region.
   `scripts/05_check_full_receptor_clashes.py` found no clashes with domains I, II or IV in the tethered 1YY9 conformation
   (closest approach 19.5 A or more). I had first described 1YY9 as the extended conformation; its own coordinates show the
   domain II tether loop touching domain IV (2.6 A), i.e. tethered. The extended (ligand-bound) conformation was not tested.

---

# Challenge 2 (TNF-alpha): a negative result

Goal: a de novo binder of the soluble human TNF-alpha trimer that also binds mouse TNF-alpha at pH 7.4 and releases at pH 6.0.
No designs were submitted. Run summaries are in `results/tnf_boltzgen_runs.csv`.

10. **Pick hotspots that fit the binder.** The first run asked for 8 hotspots on two chains up to 36 A apart; an 80-120 residue
    binder is only about 30-35 A across. The EGFR hotspots that worked were all within about 15 A.
11. **A partial target creates a loophole.** Designing against only chains A and C of the trimer gave good-looking scores
    (23 of 60 designs with ipTM 0.5 or more) but in 19 of 20 finalists 94-100% of the binder atoms lay within 4 A of the
    omitted chain B: the binders had occupied the third protomer's place. Always check clashes with every part of the real target.
12. **Site polarity matters.** Mapping the exposed surface showed the antibody-binding groove core is polar and charged
    (the E23 patch is 3% apolar, E110 10%, Q67 15%), while the best patches on the protomer have 250-310 A^2 of apolar area.
    Targeting the greasy halves of the groove instead (chain A F144/A145/S147 or chain C S86/Y87/Q88/V91) still produced nothing.
13. **A diagnostic run answered the question.** With no hotspots and the binding-site filter off, 80 designs against the trimer
    still gave a best ipTM of 0.45 and ipSAE of 0.04. So the target, not the site choice, was the limit for BoltzGen in this free-tier setup
    (about 50-110 minutes per run). On EGFR the same tool gave 34 passing designs out of 96.
14. **Set a stop rule before spending runs.** I set one (no usable designs after the unrestricted run means stop) and followed it.
    Scripts used: `scripts/06_tnf_epitope_analysis.py`.
