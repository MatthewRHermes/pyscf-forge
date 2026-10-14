'''count4.py [lib.so ...] : wrong off-diagonal pspace elements in test_csf_4elec, per (smult, m).'''
import sys, os, ctypes, numpy as np
T = os.path.join (os.path.dirname (os.path.abspath (__file__)), '../test_csf_4elec.py')
ns = {'__name__': 'x'}
exec (compile (open (T).read (), T, 'exec'), ns)
ns['setUpModule'] ()
import pyscf.csf_fci.csf as csf
sol = ns['sol']; nel = (ns['neleci'], ns['nelec'])
libs = sys.argv[1:] or [os.path.join (os.path.dirname (os.path.abspath (__file__)), '../../../lib/libcsf.so')]
for so in libs:
    csf.libcsf = ctypes.CDLL (so)
    out = []
    for smult in range (1, ns['smult_lim']):
        n = sum (nel[smult % 2]); s2 = smult - 1
        for m in range (-s2, s2+1, 2):
            ne = [(n+m)//2, (n-m)//2]
            addr, h0 = sol.pspace (ns['h1e'], ns['g2e'], ns['norb'], ne, smult=smult)
            ref = ns['get_h2mat_ref'] (ne, smult)[addr,:][:,addr]
            d = np.abs (np.tril (h0 - ref, -1))
            out.append ('s%d m%+d: %4d/%4d' % (smult, m, (d > 1e-8).sum (), len (addr)*(len (addr)-1)//2))
    print (os.path.basename (os.path.dirname (so)) + '/' + os.path.basename (so), ' | '.join (out))
