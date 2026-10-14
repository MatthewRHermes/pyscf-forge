import sys, ctypes, itertools, collections
import numpy as np
from pyscf.fci import cistring
from pyscf.csf_fci.csfstring import CSFTransformer
lib = ctypes.CDLL (sys.argv[1])
class Str3 (ctypes.Structure):
    _fields_ = [('dconf', ctypes.c_uint64), ('sconf', ctypes.c_uint64), ('spin', ctypes.c_uint64)]
lib.csf_EaiEbj.restype = ctypes.c_double
lib.csf_EaiEbj.argtypes = [ctypes.POINTER (Str3), ctypes.POINTER (Str3)] + [ctypes.c_uint]*4
norb = int (sys.argv[2]) if len (sys.argv) > 2 else 6

def Emat_string (norb, ne, p, q):
    strs = cistring.make_strings (range (norb), ne)
    n = len (strs); M = np.zeros ((n, n))
    for k, s in enumerate (strs):
        s = int (s)
        if not (s >> q) & 1: continue
        s1 = s ^ (1 << q); sg = (-1) ** bin (s1 & ((1 << q) - 1)).count ('1')
        if (s1 >> p) & 1: continue
        sg *= (-1) ** bin (s1 & ((1 << p) - 1)).count ('1')
        s2 = s1 | (1 << p)
        M[cistring.str2addr (norb, ne, s2), k] += sg
    return M

def run (na, nb, smult, tally, examples):
    t = CSFTransformer (norb, na, nb, smult)
    ncsf = t.ncsf
    if ncsf < 2: return
    U = np.stack ([t.vec_csf2det (np.eye (ncsf)[k]).ravel () for k in range (ncsf)], axis=1)
    Ia, Ib = np.eye (t.ndeta), np.eye (t.ndetb)
    Ecache = {}
    def E (p, q):
        if (p,q) not in Ecache:
            Ecache[(p,q)] = np.kron (Emat_string (norb, na, p, q), Ib) + np.kron (Ia, Emat_string (norb, nb, p, q))
        return Ecache[(p,q)]
    EEcache = {}
    def EE (a, i, b, j):
        k = (a,i,b,j)
        if k not in EEcache: EEcache[k] = U.T @ E (a,i) @ E (b,j) @ U
        return EEcache[k]
    d, s, sp = t.csfaddrs2str (np.arange (ncsf))
    d = np.maximum (d, 0)
    strs = []
    for k in range (ncsf):
        x = Str3 (int (d[k]), int (s[k]), int (sp[k]) if s[k] else 0)
        strs.append (x)
    occ = np.zeros ((ncsf, norb), dtype=int)
    for k in range (ncsf):
        dd = int (d[k]); ss = int (s[k]); free = [o for o in range (norb) if not (dd >> o) & 1]
        for o in range (norb):
            if (dd >> o) & 1: occ[k,o] = 2
        for ii, o in enumerate (free):
            if (ss >> ii) & 1: occ[k,o] = 1
    for kb, kk in itertools.permutations (range (ncsf), 2):
        diff = occ[kb] - occ[kk]
        if np.abs (diff).sum () != 4: continue
        As = [o for o in range (norb) for _ in range (max (diff[o], 0))]
        Is = [o for o in range (norb) for _ in range (max (-diff[o], 0))]
        a, b = As; i, j = Is
        if a == b or i == j: continue  # pair hops go to csf_EaiEaj; not classified here
        for (a_, i_, b_, j_) in ((a,i,b,j), (a,j,b,i)):
            # canonicalize as csf_EaiEbj does
            B, K, A_, I_, B_, J_ = kb, kk, a_, i_, b_, j_
            if min (A_, B_) > min (I_, J_): B, K, A_, I_, B_, J_ = K, B, I_, A_, J_, B_
            if A_ > B_: A_, I_, B_, J_ = B_, J_, A_, I_
            T = lib.csf_EaiEbj (ctypes.byref (strs[B]), ctypes.byref (strs[K]), A_, I_, B_, J_)
            c1 = EE (A_, I_, B_, J_)[B, K]; c2 = EE (A_, J_, B_, I_)[B, K]
            r = min (I_, J_, B_)
            if r == I_: kind, Ue, Le = 'nolink', c1, 0.0
            elif r == J_: kind, Ue, Le = 'phph', -0.5*c2, c1 + 0.5*c2
            else: kind, Ue, Le = 'pphh', 0.5*(c1+c2), 0.5*(c1-c2)
            cands = [('ok', Ue+Le), ('L_sign', Ue-Le), ('U_sign', -Ue+Le), ('both_sign', -Ue-Le),
                     ('L_missing', Ue), ('U_missing', Le)]
            verdict = 'other'
            for name, v in cands:
                if abs (T - v) < 1e-8: verdict = name; break
            nA = -occ[B][A_]; nB_ = -occ[B][B_]; nI = occ[K][I_]; nJ = occ[K][J_]
            order = ''.join (sorted ([(A_,'a'),(I_,'i'),(B_,'b'),(J_,'j')]) [k][1] for k in range (4))
            gap = int (sum (occ[K][o]==1 and occ[B][o]==1 for o in range (min (A_,I_,B_,J_)+1, max (A_,I_,B_,J_))))
            key = (kind, order, (nA, nI, nB_, nJ), 'i>j' if I_ > J_ else 'i<j')
            tally[key][verdict] += 1
            o4 = sorted ([A_,I_,B_,J_])
            gaps = tuple (int (sum (occ[K][o]==1 and occ[B][o]==1 for o in range (lo+1, hi))) for lo, hi in zip (o4[:-1], o4[1:]))
            below = int (sum (occ[K][o]==1 for o in range (0, o4[0])))
            above = int (sum (occ[K][o]==1 for o in range (o4[-1]+1, norb)))
            records.append (dict (kind=kind, order=order, occ=(int(nA),int(nI),int(nB_),int(nJ)), ij=key[3],
                gaps=gaps, below=below, above=above, smult=smult, verdict=verdict, T=T, U=Ue, L=Le,
                bra=t.printable_csfstring (B), ket=t.printable_csfstring (K), op=(A_,I_,B_,J_)))
            if verdict != 'ok' and len (examples[key]) < 2:
                examples[key].append ((smult, (na,nb), t.printable_csfstring (B), t.printable_csfstring (K),
                                       (A_,I_,B_,J_), round (T,6), round (Ue,6), round (Le,6)))
records = []
tally = collections.defaultdict (collections.Counter); examples = collections.defaultdict (list)
for n in range (2, 2*norb-1):
    for smult in range (1 + n%2, min (n, 2*norb-n) + 2, 2):
        s2 = smult - 1
        na, nb = (n+s2)//2, (n-s2)//2
        run (na, nb, smult, tally, examples)
tot = collections.Counter ()
for key in sorted (tally):
    c = tally[key]; tot.update (c)
    if set (c) != {'ok'}:
        print (key, dict (c))
        for e in examples[key]: print ('      e.g. smult,ne,bra,ket,(a,i,b,j),C,Uexact,Lexact =', e)
print ('ALL-OK classes:', sum (1 for k in tally if set (tally[k]) == {'ok'}), 'of', len (tally))
print ('totals:', dict (tot))

import pickle
pickle.dump (records, open (sys.argv[3], 'wb'))
