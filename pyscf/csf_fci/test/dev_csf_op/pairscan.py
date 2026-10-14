'''pairscan.py lib.so [norb] : graphical csf_EaiEaj (via csf_EaiEbj) vs csf_EaiEbj_ri for every pair hop.
Prints classes (ni, nj, gap parity, ...) and ratio statistics.'''
import sys, ctypes, itertools, collections, pickle
import numpy as np
from pyscf.csf_fci.csfstring import CSFTransformer
lib = ctypes.CDLL (sys.argv[1]); norb = int (sys.argv[2]) if len (sys.argv) > 2 else 6
class Str3 (ctypes.Structure):
    _fields_ = [('dconf', ctypes.c_uint64), ('sconf', ctypes.c_uint64), ('spin', ctypes.c_uint64)]
for fn in ('csf_EaiEbj', 'csf_EaiEbj_ri'):
    getattr (lib, fn).restype = ctypes.c_double
    getattr (lib, fn).argtypes = [ctypes.POINTER (Str3), ctypes.POINTER (Str3)] + [ctypes.c_uint]*4
def run2S (lbl):
    out = []; x = 0
    for ch in lbl[::-1]:
        x += {'u':1, 'd':-1}.get (ch, 0); out.append (x)
    return out
recs = []
for n in range (2, 2*norb-1):
    for smult in range (1 + n%2, min (n, 2*norb-n) + 2, 2):
        s2 = smult - 1; na, nb = (n+s2)//2, (n-s2)//2
        t = CSFTransformer (norb, na, nb, smult)
        if t.ncsf < 2: continue
        d, s, sp = t.csfaddrs2str (np.arange (t.ncsf)); d = np.maximum (d, 0)
        strs = [Str3 (int (d[k]), int (s[k]), int (sp[k]) if s[k] else 0) for k in range (t.ncsf)]
        lbl = [t.printable_csfstring (k) for k in range (t.ncsf)]
        occ = np.array ([[{'0':0,'2':2}.get (c, 1) for c in L[::-1]] for L in lbl])
        for kb, kk in itertools.permutations (range (t.ncsf), 2):
            diff = occ[kb] - occ[kk]
            if np.abs (diff).sum () != 4 or not (np.abs (diff) == 2).any (): continue
            As = [o for o in range (norb) for _ in range (max (diff[o], 0))]
            Is = [o for o in range (norb) for _ in range (max (-diff[o], 0))]
            a, b = As; i, j = Is
            if a == b and i == j: continue
            if a != b: continue   # mirror (i==j) cases are the same function with bra/ket swapped
            T = lib.csf_EaiEbj (ctypes.byref (strs[kb]), ctypes.byref (strs[kk]), a, i, b, j)
            R = lib.csf_EaiEbj_ri (ctypes.byref (strs[kb]), ctypes.byref (strs[kk]), a, i, b, j)
            lo, hi = min (i, j), max (i, j)
            sk, sb = run2S (lbl[kk]), run2S (lbl[kb])
            gap = sum (1 for o in range (lo+1, hi) if occ[kk][o] == 1)
            recs.append (dict (T=T, R=R, ni=int (occ[kk][lo]), nj=int (occ[kk][hi]), gap=gap,
                a_pos='below' if a < lo else ('between' if a < hi else 'above'),
                Sk=(sk[lo-1] if lo else 0, sk[lo], sk[hi-1], sk[hi]), Sb=(sb[lo-1] if lo else 0, sb[lo], sb[hi-1], sb[hi]),
                bra=lbl[kb], ket=lbl[kk], op=(a, i, b, j), smult=smult))
def verdict (r):
    T, R = r['T'], r['R']
    if abs (T - R) < 1e-8: return 'ok'
    if abs (R) < 1e-8: return 'should_be_0'
    return 'x%.4f' % (T / R)
for r in recs: r['v'] = verdict (r)
pickle.dump (recs, open ('/tmp/claude-1000/csfwork/pairrecs.pkl', 'wb'))
c = collections.defaultdict (collections.Counter)
for r in recs: c[(r['ni'], r['nj'])][r['v']] += 1
for k, v in sorted (c.items ()): print (k, dict (v.most_common (8)))
