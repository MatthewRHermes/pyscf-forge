# dev_csf_op: debugging harness for `pyscf/lib/csf/csf_op.c`

Scratch tooling from the csf_mem_reads debugging session (Claude-assisted). Not meant for merging.

* `build.sh [srcdir] [out.so] [extra gcc flags]` -- gcc build of libcsf (default: this worktree's
  `pyscf/lib/csf` -> `pyscf/lib/libcsf.so`). Avoids cmake configure (which clobbers core pyscf config.h).
* `env.sh cmd...` -- run with this worktree first in `PYSCF_EXT_PATH` (OMP_NUM_THREADS=1 unless set).
* `probe.py` -- "illegal read" probe.  Builds variants of csf_op.c in which the fabricated intermediate
  spin of a doubly-occupied orbital (inserted by `_get_spinindices1` = S1 / `_get_spinindices2` = S2) is
  (a) given the wrong parity, with every `parity += x` term and every 6j call checked for parity
  (asserts made non-fatal and counted per source line), and (b) given the other legal value
  (path-dependence test: a correct implementation can't depend on it).  Prints per-line hit counts and the
  number of pspace elements that change.
* `fullcheck.py [lib.so]` -- full pspace (all CSFs, hdiag on the diagonal) vs exact reference for every
  electron count/spin at norb=6 with random spin-dependent integrals; counts wrong elements by
  excitation type.  References are cached in `/tmp/claude-...` or $FULLCHECK_CACHE.
* `runtests.sh [lib.so]` -- the tests that must keep passing: test_csf_2elec, test_csf_spinless,
  test_csf_1eop, and test_csf_4elec with skips removed (only smult=4 is expected to pass).
