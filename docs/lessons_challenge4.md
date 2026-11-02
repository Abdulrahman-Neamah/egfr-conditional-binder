
---

# Challenge 4 (TfR1, pH-sensitive BBB-shuttle binder): a negative result

**Goal.** A de novo protein that binds human TfR1 at pH 7.4, releases at pH 6.0, and also
binds mouse TfR1. Goals ranked hardest first by the organizers: (a) pH switch,
(b) mouse cross-reactivity, (c) strong human binding at a functional epitope.
Epitope hint, not a rule: the apical domain, around Tyr211.

**Outcome.** Six runs, two generators, two fold classes, five sites. No design engaged the
requested epitope with credible confidence. **Nothing was submitted.** The reasoning for
not submitting is in item 32.

Run table: `results/tfr1_runs.csv`. Surface map: `results/tfr1_patch_map.csv`.
Full site screen: `results/tfr1_site_screen.csv`. Notebooks (target and spec embedded, no
designs): `notebooks/`. Reusable scripts: `scripts/07`, `scripts/08`, `scripts/09`.

---

20. **Read the organizers' selection prompt early.** It was published with the challenge.
    It ranks the goals hardest first, says a weak but clearly pH-sensitive binder beats a
    strong insensitive one, and notes that one titratable site changes affinity at most
    tenfold per pH unit - so a credible switch needs several titratable sites and a stated
    mechanism. It also says screens are read *within a class*, and that single-sequence
    predictors score antibodies and nanobodies lower. That last clause matters: a gate
    calibrated on miniproteins is the wrong gate for an Ig fold (item 29).

21. **A deposited PDB file is not necessarily the biological unit. Check REMARK 350 before
    choosing any hotspot.** This cost me three separate mistakes, and it is the single most
    useful thing in this folder. `scripts/08_biological_assembly_check.py` automates it.
    Two distinct traps:
    - *Deposited chains that are not in the assembly.* In 3S9L (transferrin-TfR1), chains B
      and D sit in the asymmetric unit but are **not** part of biological assembly 1
      (chains A + C plus a symmetry copy). Measuring against chain D made Tyr211 look like a
      transferrin contact at 3.9 A. Against its real partners it is 41.6 A away - a free
      surface. I reported the wrong conclusion to myself for about ten minutes before the
      assembly expansion corrected it.
    - *The real partner may not be in the file at all.* In 3S9L the true dimer partner of
      chain A is a symmetry copy built from the REMARK 350 transform. Naively using chain B
      put N317 and H318 at 45-48 A when they are actually **3.3-3.7 A** from the real partner.

22. **I nearly designed into the dimer interface three times.** The two highest-scoring
    patches on the whole receptor - S644 (477 A^2 apolar, 100% mouse-conserved) and W641
    (366 A^2, the greasiest at 0.56 apolar fraction) - are 4.9 A and 3.5 A from the other
    subunit, and S644 is also 3.1 A from transferrin. The patch map had been computed on a
    monomer, so buried area was scored as exposed. A conserved histidine I was briefly
    excited about, H318, turned out to be dimer-buried at 3.5 A in **both** structures.
    *Any surface property computed on a monomer is wrong for an oligomer.*

23. **The transferrin footprint is purely helical and protease-like - the apical domain is
    free.** Against the correct assembly, TfR1 chain A contacts its own transferrin through
    L619, V622, R623, N626, R629, Q640, Y643, S644, R646, G647, D648, F650, R651, T658,
    G661-D667, V670, N758, plus R121/Y123/D125. **No apical residues.** So the hinted
    epitope is non-competitive with transferrin after all, and my earlier note that I
    avoided the helical domain "from memory" was the right call reached the wrong way -
    guessing happened to agree with the data, which is not a method.

24. **The real obstacle was a squeeze between two goals, visible in the data.** On the
    apical domain the greasy patch is not conserved and the conserved patches are not greasy:

    | patch | apolar area | apolar fraction | mouse-conserved |
    |---|---|---|---|
    | Y211 (the hint) | 324 A^2 | 0.56 | **0.29** |
    | K189 | 239 A^2 | 0.42 | 0.75 |
    | F321 | 211 A^2 | 0.47 | 0.86 |
    | K382 | 178 A^2 | 0.33 | **1.00** |

    For comparison, the GLP-1R groove that worked in Challenge 3 was 360-440 A^2 apolar, and
    the TNF-alpha sites that failed were 24-133 A^2. The usable TfR1 patches sit in the
    middle band, and nothing in the apical domain clears the bar that has predicted success
    so far. `results/tfr1_site_screen.csv` now scores every exposed residue on all four
    criteria (see items 33 and 34 for the corrected method): of 230 exposed centres, 74 have
    a patch member at the dimer interface, 46 have one in the transferrin footprint, and 53
    survive.

25. **`boltzgen check` verifies cropped-target indices**, by writing `<spec>.cif` with the
    binding-site residues at B-factor 80. Reading that file confirmed indices count from the
    first residue of each chain. Note that **Tamarind's web form uses author numbering
    instead** - it reports "Valid range: 122-760" and accepts 189,192,380-382 directly.
    Two conventions for the same job; do not carry one into the other.

26. **Crop honestly, and prove it.** A crop that removes a segment next to the epitope leaves
    a hole the generator will fill, producing a binder that cannot exist on the intact
    protein. `scripts/09_build_target_crop.py` keeps the core plus everything within a
    radius of the hotspots, writes each contiguous piece as its own chain, and **reports
    anything near the site it could not keep**. My first crop silently omitted 19 residues
    within 15 A of a candidate site, one of them at 8.5 A.

27. **Verify hotspot identity against the page, not just your own parse.** Tamarind's
    sequence panel labels each 20-residue block with the number of its **last** residue, not
    its first. Reading it the other way made K189/V192 look like K188/K192 and briefly
    convinced me the whole selection was off by one. Confirming letters from two independent
    sources (my parse of the file, and the live form's own display) resolved it in a minute.

28. **"Bound with high confidence" and "bound where I asked" are different questions.**
    BoltzGen's `bindsite_under_NRMSD` columns give the fraction of the *requested* residues a
    design actually touches. In run 4 the top five designs scored up to iiPTM 0.511 while
    engaging **zero** of the five requested residues - confident binders at the wrong
    epitope, which is also why 0 of 50 passed the binding-site filter. A submission built on
    those would have failed the organizers' own co-fold check of the stated epitope.

29. **Fold class is a real lever, and our own output hinted at it.** In miniprotein mode the
    top designs came back with VHH framework 1 (`EVQLVESGGGLVQ...GGSLR`) - the model reached
    for an Ig fold unprompted. Miniproteins press a rigid helical face into a groove and need
    greasy concave sites; antibodies and nanobodies reach flat polar surfaces with a
    protruding CDR3. Every real TfR1 brain shuttle (8D3, OX26, Denali's ETV) is
    antibody-based, not a miniprotein. Switching to nanobody mode moved site engagement from
    0 to 4 of 5 residues - the first time anything touched the epitope. Nanobody mode also
    samples 7 distinct VHH frameworks, which fixed a separate problem: the miniprotein run
    had collapsed onto a single scaffold family, so its 50 designs were barely diverse.

30. **A high score with a degenerate sequence is an artifact. Test it with a control.** The
    nanobody run's best design scored iiPTM 0.862 with min PAE 2.26 - by far our best number
    on this target - with CDR1 `GGGGGGYGMA` and CDR3 `AQGVSAGGGGGGG`. Poly-glycine has no
    side chains, so no hydrophobic burial, no side-chain hydrogen bonds, no salt bridges; it
    is also maximally flexible, so it pays an entropic cost on binding that the predictor
    does not charge it for. Real nanobody CDR3s are rich in aromatics and charged residues.
    **The clean test:** rerun identically with glycine banned (`Omit Amino Acids = CG`). The
    CDRs came back genuinely chemical (`ATTSSTSEEALKDVSSWDY`, with W, Y, D, E, K) and the
    score collapsed to iiPTM 0.337, min PAE 16.19, site engagement back to zero. One
    variable, one control: the 0.862 was bought entirely with the artifact. Watch for the
    substituted version too - a few designs came back poly-alanine instead.

31. **Set a stop rule and keep a budget.** The gate was >=3 designs at ipTM >= 0.5; it was
    never met. Free hosted compute resets monthly and is finite - 50 designs cost 2 of 10
    monthly jobs here, so "more designs" is a budget decision, not a free one. New tools can
    also simply be broken: RFD4-Proteina looked the best-matched generator available
    (its "Structures to Also Bind" field expresses cross-species binding as a *constraint*,
    which is goal (b) directly), but its hotspot selector never enabled on any target, so it
    could not be used. Reported upstream.

32. **Why no submission.** The honest summary is that I could not find a site satisfying all
    three goals with the tools available. Submitting 20 designs whose declared epitope they
    demonstrably do not touch would have failed the organizers' own epitope check, and the
    selection prompt states that unverifiable claims weigh little. A no-submission costs
    nothing; a submission that contradicts its own method statement costs credibility. The
    negative result, with the screen that produced it, is the more useful artifact.

33. **Check that a hotspot is solvent-exposed before naming it. I did not, and it probably
    cost the campaign.** The five residues given to the generator for Site P were chosen on
    the *apolar area of the surrounding patch* - the neighbourhood around K189 totals
    239 A^2 - without ever asking how exposed the named residues themselves were. Computing
    relative solvent accessibility on the isolated protomer (Shrake-Rupley, Tien et al. 2013
    maximum ASA) gives:

    | residue | RSA |
    |---|---|
    | **K189** | **0.143** |
    | **V192** | **0.154** |
    | V380 | 0.371 |
    | L381 | 0.389 |
    | K382 | 0.452 |

    K189 and V192 are ~14-15% exposed: mostly buried in their own protomer, and unreachable
    for a binder. Two of five hotspots were therefore impossible to satisfy. That is a
    plausible part of the failure mode we kept seeing - high-confidence designs that engaged
    none of the requested residues (item 28) - because the generator was asked to contact
    things it could not reach and did the next best thing. `07_tfr1_site_screen.py` now gates
    on RSA >= 0.25, and `09_build_target_crop.py` prints each hotspot's RSA and refuses to
    pass over a buried one silently.

34. **Aggregate distances over the patch, conservation over the residue.** Two screens of
    the same surface disagreed on 57 of 230 patches, purely over how to aggregate: judging
    the centre residue alone gave 73 usable patches, disqualifying a patch if ANY member
    failed gave 20. Neither is right.
    - For **distance**, aggregate over the whole patch. A binder covers the neighbourhood,
      so a member buried at the dimer interface really does block it. 20 of the 57 flips were
      of this kind, and the centre-only screen missed every one.
    - For **conservation**, judge the centre. Requiring the patch *mean* punishes a perfectly
      good conserved hotspot for having variable fringes (35 of the 57 flips), and - worse -
      it can pass a hotspot that differs in mouse itself because its neighbours are conserved.
      E272 and T536 both passed that way. Report the patch mean; do not gate on it.

    The screen in this folder uses: RSA >= 0.25 and identity at the centre, minimum distance
    over all patch members, and patch mean conservation reported only.

35. **An independent reimplementation is worth more than a re-read.** The corrected screen
    was written twice, from the structures up, by two different agents working from the same
    spec but no shared code. They agree on every RSA to three decimals, on all 230 exposed
    centres, and on the 74/46 interface counts. That agreement is why the numbers above are
    trustworthy; re-reading my own script would have reproduced my own blind spot, as it did
    for three runs. The second implementation also produced evidence I had not: deposited
    chains B and D of 3S9L make 183 and 157 heavy-atom contacts with chain A, against 1422
    for the true symmetry partner - a quantitative reason to distrust the asymmetric unit,
    where I had only the REMARK 350 annotation.
