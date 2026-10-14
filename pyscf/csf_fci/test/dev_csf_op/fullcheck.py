'''fullcheck.py [lib.so ...] : count wrong pspace elements vs exact, by excitation type.'''
import sys, os, ctypes, collections, pickle, hashlib
import numpy as np
from pyscf import gto, ao2mo
import pyscf.csf_fci.csf as csfmod
from pyscf.csf_fci import csf_solver
from pyscf.csf_fci.csfstring import CSFTransformer
NORB = int (os.environ.get ('FULLCHECK_NORB', 6))
rng = np.random.default_rng (11)
h = rng.standard_normal ((2, NORB, NORB)); h1e = h + h.transpose (0,2,1)
npair = NORB*(NORB+1)//2
g2e = ao2mo.restore (1, rng.standard_normal (npair*(npair+1)//2), NORB)
if os.environ.get ('FULLCHECK_G2E') == '0': g2e[:] = 0
if os.environ.get ('FULLCHECK_H1S') == '0': h1e[1] = h1e[0]
TAG = 'g%s_s%s' % (os.environ.get ('FULLCHECK_G2E', '1'), os.environ.get ('FULLCHECK_H1S', '1'))
mol = gto.M (atom='H 0 0 0; H 0 0 1', verbose=0)
sol = csf_solver (mol, smult=1)
def systems ():
    out = []
    for n in range (1, 2*NORB):
        for smult in range (1 + n%2, min (n, 2*NORB-n) + 2, 2):
            s2 = smult - 1
            for m in sorted ({s2, s2-2 if s2 >= 2 else s2}):
                out.append ((n, smult, ((n+m)//2, (n-m)//2)))
    return out
def occ_of (t):
    d, s, sp = t.csfaddrs2str (np.arange (t.ncsf)); d = np.maximum (d, 0)
    occ = np.zeros ((t.ncsf, NORB), dtype=int)
    for k in range (t.ncsf):
        dd = int (d[k]); ss = int (s[k]); free = [o for o in range (NORB) if not (dd >> o) & 1]
        for o in range (NORB):
            if (dd >> o) & 1: occ[k,o] = 2
        for ii, o in enumerate (free):
            if (ss >> ii) & 1: occ[k,o] = 1
    return occ
def types_of (occ):
    diff = occ[:,None,:] - occ[None,:,:]
    nmove = np.abs (diff).sum (-1) // 2
    pair = (np.abs (diff) == 2).any (-1)
    typ = np.full (nmove.shape, '>2', dtype=object)
    typ[nmove == 0] = 'spin'; typ[nmove == 1] = 'single'
    typ[(nmove == 2) & pair] = 'pairhop'; typ[(nmove == 2) & ~pair] = 'double'
    np.fill_diagonal (typ, 'hdiag')
    return typ
def reference ():
    CACHE = os.environ.get ('FULLCHECK_CACHE', '/tmp/claude-1000/fullcheck_cache_norb%d_%s.pkl' % (NORB, TAG))
    if os.path.exists (CACHE): return pickle.load (open (CACHE, 'rb'))
    ref = {}
    for key in systems ():
        n, smult, ne = key
        t = CSFTransformer (NORB, ne[0], ne[1], smult)
        if t.ncsf == 0: continue
        h2eff = sol.absorb_h1e (h1e, g2e, NORB, ne, .5)
        ndet = t.ndeta*t.ndetb
        hdet = np.zeros ((ndet, ndet))
        for i in range (ndet):
            c = np.zeros (ndet); c[i] = 1
            hdet[i] = sol.contract_2e (h2eff, c.reshape (t.ndeta, t.ndetb), NORB, ne).ravel ()
        U = np.stack ([t.vec_csf2det (np.eye (t.ncsf)[k]).ravel () for k in range (t.ncsf)], axis=1)
        ref[key] = (U.T @ hdet @ U, types_of (occ_of (t)), [t.printable_csfstring (k) for k in range (t.ncsf)])
    pickle.dump (ref, open (CACHE, 'wb'))
    return ref
def run (libpath):
    if not (getattr (csfmod.libcsf, '_name', None) == libpath): csfmod.libcsf = ctypes.CDLL (libpath)
    out = {}
    for key in systems ():
        n, smult, ne = key
        t = CSFTransformer (NORB, ne[0], ne[1], smult)
        if t.ncsf == 0: continue
        _, h0 = sol.pspace (h1e, g2e, NORB, ne, smult=smult, csf_addr=np.arange (t.ncsf, dtype=np.int32), npsp=t.ncsf)
        out[key] = h0
    return out
def score (res, ref, verbose=False):
    cnt = collections.Counter (); tot = collections.Counter (); ex = {}
    for key, (R, typ, lbl) in ref.items ():
        A = res[key]
        bad = np.abs (np.tril (A - R)) > 1e-8
        low = np.tril (np.ones_like (bad))
        for tp in ('hdiag', 'spin', 'single', 'pairhop', 'double'):
            m = (typ == tp) & (low > 0)
            tot[tp] += int (m.sum ()); cnt[tp] += int ((bad & m).sum ())
            if verbose and (bad & m).any () and tp not in ex:
                i, j = np.argwhere (bad & m)[0]
                ex[tp] = (key, lbl[i], lbl[j], float (A[i,j]), float (R[i,j]))
    return cnt, tot, ex
if __name__ == '__main__':
    ref = reference ()
    libs = sys.argv[1:] or [os.path.join (os.path.dirname (os.path.abspath (__file__)), '../../../lib/libcsf.so')]
    for lp in libs:
        cnt, tot, ex = score (run (lp), ref, verbose=True)
        print ('%s  wrong/total: ' % os.path.basename (os.path.dirname (lp) or lp) +
               '  '.join ('%s %d/%d' % (tp, cnt[tp], tot[tp]) for tp in ('hdiag','spin','single','pairhop','double')))
        for tp, e in ex.items (): print ('    e.g. %-7s %s' % (tp, e))
