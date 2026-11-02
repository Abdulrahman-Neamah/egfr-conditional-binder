"""
Crop a large receptor down to a design target WITHOUT creating empty space next to the site.

Why this exists
  Cropping speeds up design, but if you cut away a segment that sits next to your epitope,
  the generator happily fills the hole it left and you get a binder that cannot exist on
  the intact protein. This script keeps a named core region plus EVERY other residue of the
  chain within --radius of the chosen hotspots, splits what it keeps into contiguous
  segments, and writes each segment as its own chain so no backbone is built across a gap.

  It also reports anything within --check of the hotspots that it could not keep, so you
  can see whether the crop is honest before spending GPU time.

Residue numbering
  Output chains are renumbered from 1. BoltzGen YAML counts residues from the first residue
  of each chain, so the printed map gives the index to put in `binding:`. Verify with
  `boltzgen check <spec>.yaml`, which writes <spec>.cif with binding-site residues at
  B-factor 80. (Tamarind's web form takes author numbering instead - do not mix the two.)

Usage
  python 09_build_target_crop.py --pdb 11QC.pdb --chain A --core 189-383 \
      --hotspots 189,192,380,381,382 --radius 20 --out TfR1_crop.pdb
"""
import argparse, collections, csv, math

AA3 = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
       'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
ONE = dict(zip(sorted(AA3), "ARNDCQEGHILKMFPSTWYV"))
ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--core", required=True, help="first-last author numbering, e.g. 189-383")
    ap.add_argument("--hotspots", required=True, help="comma-separated author residue numbers")
    ap.add_argument("--radius", type=float, default=20.0, help="context radius, angstrom")
    ap.add_argument("--check", type=float, default=15.0, help="warn about omissions within this")
    ap.add_argument("--min-segment", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    res = collections.OrderedDict()
    for L in open(a.pdb):
        if L.startswith("ATOM") and L[21] == a.chain and L[17:20].strip() in AA3:
            res.setdefault(int(L[22:26]), {"name": L[17:20].strip(), "lines": []})["lines"].append(L)

    def atoms(n):
        return [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in res[n]["lines"]]

    hot = [int(x) for x in a.hotspots.split(",")]
    missing_hot = [h for h in hot if h not in res]
    if missing_hot:
        raise SystemExit(f"hotspots absent from chain {a.chain}: {missing_hot}")
    hot_atoms = [p for h in hot for p in atoms(h)]

    lo, hi = (int(x) for x in a.core.split("-"))
    keep = {n for n in res if lo <= n <= hi}
    for n in res:
        if n not in keep and any(math.dist(p, q) <= a.radius for p in atoms(n) for q in hot_atoms):
            keep.add(n)

    ks = sorted(keep)
    segs, cur = [], [ks[0]]
    for x, y in zip(ks, ks[1:]):
        (cur.append(y) if y == x + 1 else (segs.append(cur), cur := [y]))
    segs.append(cur)
    segs = [s for s in segs if len(s) >= a.min_segment]
    segs.sort(key=lambda s: (0 if s[0] <= lo <= s[-1] else 1, s[0]))

    out, index, rows = [], {}, []
    for si, seg in enumerate(segs):
        cid = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[si]
        for i, n in enumerate(seg, 1):
            index[n] = (cid, i)
            out += [L[:21] + cid + f"{i:4d}" + L[26:] for L in res[n]["lines"]]
        out.append("TER\n")
        rows.append(dict(new_chain=cid, source_chain=a.chain, first=seg[0], last=seg[-1],
                         residues=len(seg),
                         note="core region (hotspots here)" if seg[0] <= lo <= seg[-1]
                              else "neighbouring segment kept for context"))
    open(a.out, "w").writelines(out + ["END\n"])
    with open(a.out.replace(".pdb", "_chain_map.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    omitted = []
    for n in res:
        if n in keep:
            continue
        d = min(math.dist(p, q) for p in atoms(n) for q in hot_atoms)
        if d < a.check:
            omitted.append((n, ONE[res[n]["name"]], round(d, 1)))

    print(f"kept {sum(len(s) for s in segs)} residues in {len(segs)} segments -> {a.out}")
    for r in rows:
        print("   ", r)
    print("\nhotspot author -> crop (chain, index, expected letter):")
    for n in hot:
        cid, i = index[n]
        print(f"   {ONE[res[n]['name']]}{n} -> chain {cid} index {i}")
    print(f"   BoltzGen `binding:` for chain {index[hot[0]][0]} = "
          f"{','.join(str(index[n][1]) for n in sorted(hot))}")
    print(f"\nresidues within {a.check} A of the hotspots that were NOT kept: "
          f"{omitted if omitted else 'none - crop is clean'}")
    if omitted:
        print("   WARNING: increase --radius; a binder may fill the space these would occupy.")


if __name__ == "__main__":
    main()
