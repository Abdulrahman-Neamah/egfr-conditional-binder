"""Inter-chain backbone clash check for a predicted complex.

Two polypeptide backbones cannot occupy the same space. Any predicted complex whose
chains share backbone atoms is not a complex, however confident the docking score.
ipTM-family scores measure how consistently two chains are placed relative to each
other; none of them contains a steric term, so an interpenetrating pose can score well.

Usage:  python 10_interchain_backbone_clash.py <structure.pdb|.cif> <chainA> <chainB>

Written after the rank-1 Challenge 4 design (iiPTM 0.862, min PAE 2.26 A) turned out to
have its two poly-glycine CDR loops threaded through the target: 248 backbone-atom pairs
under 2.5 A, closest 0.41 A. Atoms left at the origin are treated as unplaced and skipped.
"""
import sys, warnings; warnings.filterwarnings('ignore')
import numpy as np
from scipy.spatial import distance
from Bio.PDB import MMCIFParser, PDBParser

path = sys.argv[1]; c1, c2 = sys.argv[2], sys.argv[3]
P = MMCIFParser if path.endswith(('.cif', '.mmcif')) else PDBParser
m = P(QUIET=True).get_structure('x', path)[0]
BB = ('N', 'CA', 'C', 'O')
g = lambda c: np.array([a.coord for r in m[c] for a in r
                        if a.get_id() in BB and not np.allclose(a.coord, 0, atol=1e-3)])
D = distance.cdist(g(c1), g(c2))
print(f'{path}  {c1} vs {c2}')
print(f'  closest backbone atoms : {D.min():.2f} A')
for cut in (2.0, 2.5, 3.0):
    print(f'  backbone pairs < {cut} A : {(D < cut).sum()}')
print('  VERDICT:', 'INTERPENETRATING - not a physical complex'
      if D.min() < 2.2 else 'no backbone interpenetration')
