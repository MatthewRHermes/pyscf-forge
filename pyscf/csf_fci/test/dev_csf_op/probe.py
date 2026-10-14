'''probe.py [srcdir] : detect code that reads the formally-undefined pair intermediate spin.
See README.md.  Run under env.sh.'''
import sys, os, re, shutil, subprocess, tempfile, collections, ctypes
import numpy as np
HERE = os.path.dirname (os.path.abspath (__file__))
sys.path.insert (0, HERE)
import fullcheck as fc
import pyscf.csf_fci.csf as csfmod
_ARGS = [x for x in sys.argv[1:] if not x.startswith ('--')]
SRC = _ARGS[0] if _ARGS else os.path.join (HERE, '../../../lib/csf')
TOP = r'''
#include <stdio.h>
static long softhits[8000];
static char ctx_call[256] = "";
static char ctx_cgc[512] = "";
static int ctx_nex = 3;
static void _soft_hit (int line){
    softhits[line]++;
    if (softhits[line] <= ctx_nex){ fprintf (stderr, "EX line %d | %s | %s\n", line, ctx_call, ctx_cgc); }
}
#define SOFTASSERT(c) do { if (!(c)) { _soft_hit (__LINE__); } } while (0)
static void _ctx_arrays (const char * fn, unsigned int * k, unsigned int * b, unsigned int n,
                         int i1, int i2, int i3, int i4, int n1, int n2, int n3, int n4){
    int o = snprintf (ctx_cgc, 512, "%s idx=(%d,%d,%d,%d) n=(%d,%d,%d,%d) Sk=", fn, i1, i2, i3, i4, n1, n2, n3, n4);
    for (unsigned int x=0; x<=n && o<480; x++) o += snprintf (ctx_cgc+o, 512-o, "%u.", k[x]);
    o += snprintf (ctx_cgc+o, 512-o, " Sb=");
    for (unsigned int x=0; x<=n && o<500; x++) o += snprintf (ctx_cgc+o, 512-o, "%u.", b[x]);
}
static long w6jbad[8000];
double _w6j_impl (unsigned int j1, unsigned int j2, unsigned int j3, unsigned int j5, unsigned int j6);
static double w6j_logged (int line, unsigned int j1, unsigned int j2, unsigned int j3, unsigned int j5, unsigned int j6){
    if (((j1+j2+j3)%2) || ((j1+j5+j6)%2) || ((j2+j6+1)%2) || ((j5+j3+1)%2)){ w6jbad[line]++; if (w6jbad[line] <= ctx_nex) fprintf (stderr, "EX6J line %d | %s | %s | 6j(%u,%u,%u,1,%u,%u)\n", line, ctx_call, ctx_cgc, j1, j2, j3, j5, j6); }
    return _w6j_impl (j1, j2, j3, j5, j6);
}
void soft_get (long * a, long * b){ for (int i=0;i<8000;i++){ a[i] = softhits[i]; b[i] = w6jbad[i]; } }
void soft_set_nex (int n){ ctx_nex = n; }
void soft_reset (void){ for (int i=0;i<8000;i++){ softhits[i] = 0; w6jbad[i] = 0; } }
#define _get_wigner_6j_j41h(a,b,c,d,e) w6j_logged (__LINE__, a, b, c, d, e)
'''
PAT = re.compile (r'^(\s*)parity \+= (.*);(.*)$')
HOOKS = [
 ("    double xdiag = 1.0;\n    unsigned int twoS0b, twoS0k;\n    int offk = 0;\n    int offb = 0;\n    int parity = 0;\n    unsigned int i;\n",
  "    double xdiag = 1.0;\n    unsigned int twoS0b, twoS0k;\n    int offk = 0;\n    int offb = 0;\n    int parity = 0;\n    unsigned int i;\n    _ctx_arrays (\"X\", twoSk, twoSb, nspin, t, q, r, p, nt, nq, nr, np);\n"),
 ("    int parity = 0;\n    double xdiag = 1.0;\n    unsigned int twoS0k, twoS0b;\n",
  "    int parity = 0;\n    double xdiag = 1.0;\n    unsigned int twoS0k, twoS0b;\n    _ctx_arrays (\"1e\", twoSk, twoSb, nspin, q, p, 0, 0, nq, np, 0, 0);\n"),
 ("    Str3 brap, ketp;\n    _pad_Str3 (bra, &brap);\n    _pad_Str3 (ket, &ketp);\n    // CSF orthogonality\n    unsigned int p = a;\n",
  "    snprintf (ctx_call, 256, \"EaiEbj <%lu,%lu,%lu|%u'%u %u'%u|%lu,%lu,%lu>\", bra->dconf, bra->sconf, bra->spin, a, i, b, j, ket->dconf, ket->sconf, ket->spin);\n    Str3 brap, ketp;\n    _pad_Str3 (bra, &brap);\n    _pad_Str3 (ket, &ketp);\n    // CSF orthogonality\n    unsigned int p = a;\n"),
]
def softify (s):
    for o, n in HOOKS:
        if s.count (o) != 1: print ('WARNING hook not applied:', o[:60].replace ('\n', '|'), file=sys.stderr)
        s = s.replace (o, n)
    s = s.replace ('double _get_wigner_6j_j41h (', 'double _w6j_impl (')
    s = re.sub (r'\bassert ?\(', 'SOFTASSERT (', s)
    out = []
    for l in s.split ('\n'):
        m = PAT.match (l)
        if m: l = '%s{ int _pt = (int) (%s); SOFTASSERT ((_pt %% 2) == 0); parity += _pt; }%s' % m.groups ()
        out.append (l)
    lines = out; last = max (k for k, l in enumerate (lines) if l.startswith ('#include'))
    return '\n'.join (lines[:last+1]) + TOP + '#line %d "csf_op.c"\n' % (last+2) + '\n'.join (lines[last+1:])
def build (tmp, name, defs, soft):
    d = os.path.join (tmp, name); os.makedirs (d)
    for f in os.listdir (SRC):
        if f.endswith (('.c', '.h')): shutil.copy (os.path.join (SRC, f), d)
    if soft:
        fn = os.path.join (d, 'csf_op.c'); txt = softify (open (fn).read ()); open (fn, 'w').write (txt)
    out = os.path.join (d, 'libcsf.so')
    flags = ['-DNDEBUG'] if soft else []
    r = subprocess.run ([os.path.join (HERE, 'build.sh'), d, out] + flags + defs, capture_output=True, text=True)
    if r.returncode: raise RuntimeError (r.stdout + r.stderr)
    return out
def softrun (lp):
    lib = ctypes.CDLL (lp); csfmod.libcsf = lib
    lib.soft_reset ()
    fc.run (lp.replace ('libcsf.so', 'libcsf.so'))  # run () reloads; keep handle consistent
    a = np.zeros (8000, dtype=np.int64); b = np.zeros (8000, dtype=np.int64)
    lib2 = csfmod.libcsf
    lib2.soft_get (a.ctypes.data_as (ctypes.c_void_p), b.ctypes.data_as (ctypes.c_void_p))
    return a, b
MID = {'bad': '(x)', 'alt': '((x)>0?(x)-1:(x)+1)'}
if __name__ == '__main__':
    if '--examples' in sys.argv:
        sys.argv.remove ('--examples')
        SRC = sys.argv[1] if len (sys.argv) > 1 else SRC
        with tempfile.TemporaryDirectory (dir='/tmp/claude-1000') as tmp:
            for k in (1, 2):
                print ('######## S%d' % k, flush=True)
                lp = build (tmp, 'S%dbadex' % k, ['-DPAIR_MID%d(x)=%s' % (k, MID['bad'])], True)
                sys.stderr.flush (); softrun (lp); sys.stderr.flush ()
        sys.exit (0)
    ref = fc.reference ()
    with tempfile.TemporaryDirectory (dir='/tmp/claude-1000') as tmp:
        ctrl = softrun (build (tmp, 'ctrl', [], True))
        src_lines = open (os.path.join (SRC, 'csf_op.c')).read ().split ('\n')
        print ('== illegal reads (per source line of %s; control-subtracted)' % os.path.join (SRC, 'csf_op.c'))
        total = 0
        for k in (1, 2):
            a, b = softrun (build (tmp, 'S%dbad' % k, ['-DPAIR_MID%d(x)=%s' % (k, MID['bad'])], True))
            for kind, x, x0 in (('parity/assert', a, ctrl[0]), ('6j', b, ctrl[1])):
                for line in np.where (x != x0)[0]:
                    total += 1
                    print ('  S%d %-13s line %4d (%6d hits): %s' % (k, kind, line, x[line]-x0[line], src_lines[line-1].strip ()[:90]))
        print ('  TOTAL flagged (slot, line) pairs:', total)
        base = fc.run (build (tmp, 'normal', [], False))
        print ('== path dependence (elements changed when the intermediate takes its other legal value)')
        for k in (1, 2):
            alt = fc.run (build (tmp, 'S%dalt' % k, ['-DPAIR_MID%d(x)=%s' % (k, MID['alt'])], False))
            cnt = collections.Counter ()
            for key, (R, typ, lbl) in ref.items ():
                ch = np.abs (np.tril (alt[key] - base[key])) > 1e-10
                for tp in np.unique (typ[ch]): cnt[tp] += int ((ch & (typ == tp)).sum ())
            print ('  S%d: %s' % (k, dict (cnt) if cnt else 'none'))
        cnt, tot, ex = fc.score (base, ref)
        print ('== wrong vs exact: ' + '  '.join ('%s %d/%d' % (tp, cnt[tp], tot[tp]) for tp in ('hdiag','spin','single','pairhop','double')))
