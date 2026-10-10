"""
TNF-alpha epitope mapping for the Anthropic x Adaptyv Challenge 2 (negative-result record; see docs/lessons.md).

Usage:   python 06_tnf_epitope_analysis.py DATA_DIR
DATA_DIR must contain 1TNF.pdb (human TNF trimer), 3WD5.pdb1 (adalimumab complex) and 4G3Y.pdb1
(infliximab complex), all downloaded from files.rcsb.org/download/.

What it does (numpy/pandas only):
  1. aligns human (soluble domain, UniProt P01375 77-233) and mouse TNF-alpha;
  2. fits the TNF protomer of each Fab complex onto every protomer of the 1TNF trimer and moves the Fab
     with it, so the epitope can be read on BOTH protomers it spans (the complexes contain one protomer);
  3. reports contacted residues (4.5 A), mouse conservation and exposure on the free trimer;
  4. ranks surface patches (10 A radius) by exposed apolar area x mouse conservation.
Writes tnf_epitope_table.csv and tnf_patch_map.csv in the working directory. Predicted/structural analysis only.
"""
import sys, pickle
import numpy as np, pandas as pd, collections, itertools
three={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
HUMAN=("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQV"
       "LFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVF"
       "QLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
MOUSE=("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTHTVSRFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")
def load(path,model=None):
    """returns OrderedDict[(chain,resnum)] -> dict(name, atoms{name:xyz}, elem{name:el})"""
    res=collections.OrderedDict(); m=0; use=True
    for l in open(path):
        if l.startswith("MODEL"): m=int(l.split()[1]); use=(model is None or m==model)
        if l.startswith(("ATOM","HETATM")) and use and l.startswith("ATOM"):
            k=(l[21],int(l[22:26]))
            r=res.setdefault(k,{"name":l[17:20].strip(),"atoms":{},"el":{}})
            n=l[12:16].strip(); r["atoms"][n]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            r["el"][n]=(l[76:78].strip() or n[0])
    return res
def chain(res,ch): return collections.OrderedDict((k[1],v) for k,v in res.items() if k[0]==ch)
def kabsch(P,Q):
    pc,qc=P.mean(0),Q.mean(0); H=(P-pc).T@(Q-qc); U,_,Vt=np.linalg.svd(H); d=np.sign(np.linalg.det(Vt.T@U.T))
    R=Vt.T@np.diag([1,1,d])@U.T; return R,qc-R@pc
def xyz(res): return np.array([a for v in res.values() for a in v["atoms"].values()])
def dmat(A,B): return np.sqrt(((A[:,None,:]-B[None,:,:])**2).sum(-1))
MAXASA={'A':106,'R':248,'N':157,'D':163,'C':135,'Q':198,'E':194,'G':84,'H':184,'I':169,'L':164,'K':205,'M':188,'F':197,'P':136,'S':130,'T':142,'W':227,'Y':222,'V':142}
RAD={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
def sasa_per_residue(res,probe=1.4,npts=96):
    keys=list(res.keys()); atoms=[];owner=[];rad=[]
    for i,k in enumerate(keys):
        for n,x in res[k]["atoms"].items(): atoms.append(x);owner.append(i);rad.append(RAD.get(res[k]["el"][n][:1],1.7))
    atoms=np.array(atoms);rad=np.array(rad)+probe;owner=np.array(owner)
    idx=np.arange(npts)+0.5; phi=np.arccos(1-2*idx/npts); th=np.pi*(1+5**0.5)*idx
    sph=np.c_[np.cos(th)*np.sin(phi),np.sin(th)*np.sin(phi),np.cos(phi)]
    D=dmat(atoms,atoms); out=np.zeros(len(keys))
    for i in range(len(atoms)):
        nb=np.where((D[i]<rad[i]+rad)&(np.arange(len(atoms))!=i))[0]
        pts=atoms[i]+rad[i]*sph
        free=np.ones(npts,bool)
        for j in nb: free&=(((pts-atoms[j])**2).sum(1)>=rad[j]**2)
        out[owner[i]]+=4*np.pi*rad[i]**2*free.mean()
    return {keys[i]:out[i] for i in range(len(keys))}

DATA=sys.argv[1] if len(sys.argv)>1 else '.'
# ---------- human/mouse alignment (global, identity scoring with gap penalty) ----------
def nw(a,b,match=2,mis=-1,gap=-3):
    n,m=len(a),len(b); S=np.zeros((n+1,m+1)); T=np.zeros((n+1,m+1),int)
    S[:,0]=np.arange(n+1)*gap; S[0,:]=np.arange(m+1)*gap; T[1:,0]=1; T[0,1:]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            c=[S[i-1,j-1]+(match if a[i-1]==b[j-1] else mis),S[i-1,j]+gap,S[i,j-1]+gap]; k=int(np.argmax(c)); S[i,j]=c[k]; T[i,j]=k
    i,j=n,m; mp={}
    while i>0 or j>0:
        if i>0 and j>0 and T[i,j]==0: mp[i]=j; i-=1; j-=1
        elif i>0 and (j==0 or T[i,j]==1): i-=1
        else: j-=1
    return mp   # human index (1-based) -> mouse index (1-based)
h2m=nw(HUMAN,MOUSE)
ident=sum(HUMAN[i-1]==MOUSE[j-1] for i,j in h2m.items())
print(f"aligned {len(h2m)} residues, identical {ident} ({ident/len(h2m):.0%})  [page: 124 of 156]")
unal=[i for i in range(1,158) if i not in h2m]; print("human residues without a mouse partner:",unal,[HUMAN[i-1] for i in unal])
# ---------- target trimer ----------
tnf=load(DATA+"/1TNF.pdb"); chains={c:chain(tnf,c) for c in "ABC"}
sas=sasa_per_residue(tnf)
rsasa={k:sas[k]/MAXASA[three[tnf[k]["name"]]] for k in tnf}
# ---------- place each Fab onto each 1TNF protomer ----------
def fab_sets(path,model,tch):
    st=load(path,model); t=chain(st,tch)
    fab={ch:chain(st,ch) for ch in "HL"}
    return t,fab
sources={"adalimumab (3WD5)":fab_sets(DATA+"/3WD5.pdb1",None,"A"),"infliximab (4G3Y)":fab_sets(DATA+"/4G3Y.pdb1",1,"C")}
rows=[]; summary={}
for name,(t,fab) in sources.items():
    fx=np.vstack([xyz(fab[c]) for c in "HL"])
    for X in "ABC":
        common=[n for n in t if n in chains[X] and "CA" in t[n]["atoms"] and "CA" in chains[X][n]["atoms"]]
        P=np.array([t[n]["atoms"]["CA"] for n in common]); Q=np.array([chains[X][n]["atoms"]["CA"] for n in common])
        R,tr=kabsch(P,Q); rms=np.sqrt(((P@R.T+tr-Q)**2).sum(1).mean())
        moved=fx@R.T+tr
        per={}
        for Y in "ABC":
            for n,v in chains[Y].items():
                c=int((dmat(np.array(list(v["atoms"].values())),moved)<=4.5).sum())
                if c: per[(Y,n)]=c
        by=collections.Counter(); 
        for (Y,n),c in per.items(): by[Y]+=1
        summary[(name,X)]=(len(common),rms,dict(by),per)
        print(f"{name}: Fab placed on protomer {X}: CA fit {len(common)} residues, rmsd {rms:.2f} A; epitope residues per chain {dict(by)}")

prev={"A":"C","B":"A","C":"B"}      # neighbour protomer contacted by a Fab placed on X (from the placements above)
def site(name):
    prim=collections.defaultdict(list); neigh=collections.defaultdict(list)
    for X in "ABC":
        per=summary[(name,X)][3]
        for (Y,n),c in per.items():
            if Y==X: prim[n].append(c)
            elif Y==prev[X]: neigh[n].append(c)
    return prim,neigh
rows=[]
for name in ("adalimumab (3WD5)","infliximab (4G3Y)"):
    prim,neigh=site(name)
    for role,d in (("primary protomer",prim),("neighbour protomer",neigh)):
        for n,cs in sorted(d.items()):
            aa=HUMAN[n-1]; m=MOUSE[h2m[n]-1] if n in h2m else "-"
            r=np.mean([rsasa[(Y,n)] for Y in "ABC" if (Y,n) in rsasa])
            rows.append(dict(antibody=name.split()[0],role=role,residue=f"{aa}{n}",placements=len(cs),mean_contacts=round(np.mean(cs),1),mouse=m,conserved=(aa==m),rsasa_free_trimer=round(float(r),2)))
df=pd.DataFrame(rows); df.to_csv("tnf_epitope_table.csv",index=False)
pd.set_option("display.width",200)
for ab in ("adalimumab","infliximab"):
    for role in ("primary protomer","neighbour protomer"):
        g=df[(df.antibody==ab)&(df.role==role)]
        print(f"\n=== {ab}, {role}: {len(g)} residues contacted in >=1 placement; seen in all 3 placements: {int((g.placements==3).sum())}; conserved in mouse: {int(g.conserved.sum())}/{len(g)}")
        print(g[g.placements>=2].sort_values("mean_contacts",ascending=False)[["residue","placements","mean_contacts","mouse","conserved","rsasa_free_trimer"]].to_string(index=False))


# ---------- surface patch map on protomer A ----------
A=chain(tnf,"A"); cons={n:(n in h2m and HUMAN[n-1]==MOUSE[h2m[n]-1]) for n in range(1,158)}
apolar=set("AVLIMFWYC"); ca={n:A[n]["atoms"]["CA"] for n in A}; ns=sorted(A)
ex={n:rsasa[("A",n)] for n in ns}; asa={n:ex[n]*MAXASA[HUMAN[n-1]] for n in ns}; prow=[]
for ctr in ns:
    if ex[ctr]<0.2: continue
    pat=[n for n in ns if np.linalg.norm(ca[n]-ca[ctr])<=10 and ex[n]>=0.2]
    ap=sum(asa[n] for n in pat if HUMAN[n-1] in apolar); tot=sum(asa[n] for n in pat)
    prow.append(dict(center=f"{HUMAN[ctr-1]}{ctr}",exposed_res=len(pat),apolar_area=int(ap),apolar_frac=round(ap/tot,2),conserved_frac=round(float(np.mean([cons[n] for n in pat])),2),score=round(ap*float(np.mean([cons[n] for n in pat])))))
pd.DataFrame(prow).sort_values("score",ascending=False).to_csv("tnf_patch_map.csv",index=False)
print("wrote tnf_epitope_table.csv and tnf_patch_map.csv")
