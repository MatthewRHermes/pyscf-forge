'''timing.py lib1.so lib2.so : pspace wall time for two builds on random CI problems (norb 10-14).'''
import sys, time, ctypes, numpy as np
from pyscf import gto, ao2mo
import pyscf.csf_fci.csf as csfmod
from pyscf.csf_fci import csf_solver
from pyscf.csf_fci.csfstring import CSFTransformer
mol = gto.M (atom='H 0 0 0; H 0 0 1', verbose=0); sol = csf_solver (mol, smult=1)
rng = np.random.default_rng (3)
for norb, ne, smult, npsp in ((10, (5,5), 1, 400), (12, (6,6), 1, 400), (12, (7,5), 3, 400), (14, (7,7), 1, 400)):
    h = rng.standard_normal ((norb, norb)); h1e = h + h.T
    npair = norb*(norb+1)//2; g2e = ao2mo.restore (1, rng.standard_normal (npair*(npair+1)//2), norb)
    t = CSFTransformer (norb, ne[0], ne[1], smult)
    addr = np.sort (rng.choice (t.ncsf, size=min (npsp, t.ncsf), replace=False)).astype (np.int32)
    res = []
    for lp in sys.argv[1:]:
        csfmod.libcsf = ctypes.CDLL (lp)
        t0 = time.perf_counter (); _, h0 = sol.pspace (h1e, g2e, norb, ne, smult=smult, csf_addr=addr, npsp=len (addr)); dt = time.perf_counter () - t0
        res.append ((dt, h0))
    print ('norb %2d ne %s smult %d ncsf %7d npsp %d: ' % (norb, ne, smult, t.ncsf, len (addr)) +
           '  '.join ('%s %.3fs' % (lp.split ('/')[-1], r[0]) for lp, r in zip (sys.argv[1:], res)) +
           '  max|diff| %.2e' % np.abs (res[0][1] - res[1][1]).max (), flush=True)
