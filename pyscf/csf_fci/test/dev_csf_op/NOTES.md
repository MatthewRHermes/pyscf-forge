# csf_mem_reads: what was done and why (2026-10-14, Claude-assisted)

Goal (from MRH): remove every read of the formally undefined running spin "between" the two
electrons of a doubly-occupied orbital, without breaking test_csf_2elec, test_csf_spinless,
test_csf_1eop, or test_csf_4elec smult=4; then debug test_csf_4elec.

## Tools (this directory)

* `probe.py` -- the undefined intermediates are stored by `_get_spinindices1` (S1) and
  `_get_spinindices2` (S2) through the `PAIR_MID1/2` hooks. The probe rebuilds csf_op.c with them set
  to a wrong-parity value, turns every assert and every `parity += x` term into a per-line counter
  of parity violations, and logs 6j calls with impossible arguments. It also sets them to the other
  legal value and counts pspace elements that change (path dependence). Positive control:
  `probe.py <copy with CSF_PPHH_GRAPHICAL=1>` still finds the 14 reads.
* `fullcheck.py` -- all-CSF pspace (with hdiag) vs exact `U^T H_det U` for every electron count and
  spin, random spin-dependent integrals, norb=6 (FULLCHECK_NORB=7 also works; use
  FULLCHECK_STREAM=1 for norb=8).
* `scan2e.py`/`scan2e_rec.py`/`scan_pair.py`/`pairscan.py` -- per-call comparisons of `csf_EaiEbj`.
  The Fierz split used: r=i (no linked graph) U=c1; phph (r=j) U=-c2/2, L=c1+c2/2;
  pphh (r=b) U=(c1+c2)/2 (singlet channel), L=(c1-c2)/2 (triplet channel); c1=<E_ai E_bj>, c2=<E_aj E_bi>.
* `-DCSF_EAIEBJ_PART=1|2` returns only the linked|unlinked term of `csf_EaiEbj`.
* `-DCSF_2E_RI=1` -- exact oracle: every Str3_link==2 element by resolution of the identity over
  the intermediate occupation's CSFs, using only `csf_Eai`.

## Baseline (ce1d5c8, norb=6 fullcheck, wrong/total)

hdiag 0/2429, spin-only 946/946, single 0/29910, pair hop 17166/27120, double 11829/41670.
probe: 14 (slot, line) illegal reads, all in csf_EaiEbj doubles and csf_EaiEaj pair hops.
Singles are 100% correct, which is the lever used below.

## Commits

1. `e0f9129` tooling + PAIR_MID hooks (no behavior change).
2. `0c666bd` record intermediate slots; DUMMY_RULE experiment. No convention for the intermediate
   (up, down, (anti)parallel to the facing chain, copy) fixes anything: the code that reads it is
   wrong regardless of the value. (The probe was accidentally blind between this commit and a27e629.)
3. `6c4efac` spin-only elements: the exchange loop started at q=p and added (pp|pp)*1 to every
   off-diagonal spin-only element (csf_EaiEaj(i==j) returns 1). 946 -> 0.
4. `8bf5cbb` `csf_EaiEbj_ri`: <bra|E_ai E_bj|ket> = sum_m <bra|E_ai|m><m|E_bj|ket> over the CSFs of
   the intermediate occupation (E is spin-free, so this is exact). With -DCSF_2E_RI=1 every test
   passes, including all of test_csf_4elec and test_csf.py / test_csf_symm.py with skips removed.
5. `a27e629` pair hops: on (i,j), E_ai E_aj and E_ji act identically (sqrt2 x singlet part, kill the
   triplet), and a closed pair can be moved between orbitals for free, so a pair hop is one
   `csf_Eai` between CSFs with one pair relocated (signs +1, +1, -1 for nj=1, ni=1, ni=nj=2,
   verified for every pair hop vs the oracle). The whole-pair move now checks that the couplings
   agree. No undefined reads remain in pair hops. 17166 -> 0.
6. `ac611f0` pphh: `csf_EaiEbj_ri_local` (the same identity, enumerating only the open shells inside
   both operator ranges; equal to the full enumeration on 202920 calls). Used for canonical pphh
   calls; `-DCSF_PPHH_GRAPHICAL=1` restores the graphical path. **Design decision for you to review.**
   probe 14 -> 0 illegal reads; all 4elec tests pass.
7. `01493f6` CGC_2e_X normalization: when two doubly-occupied operator orbitals are adjacent in
   the same chain, their "no CG coefficient" claims overlap on one index and one (2S+1) factor
   survives (the linked term came out sqrt(2S+1) too large; never the unlinked term). Counting
   claims per index and dividing the leftover fixes all remaining phph errors.
   fullcheck norb=6 and norb=7: 0 wrong in every category.
8. `da01fe1` fullcheck low-memory streaming mode. `FULLCHECK_STREAM=1 FULLCHECK_NORB=8`: all 71
   systems, 0 wrong out of ~8.38M lower-triangle elements (hdiag 36684, spin 37029, single 1444464,
   pair hop 1810116, double 5053860). This reaches spectator patterns norb=6 can't (up to 4 open
   spectators around a double). ~15 min total, largest system 36 s, memory well under 1 GB.

## State of the graphical pphh path (CSF_PPHH_GRAPHICAL=1), for future derivation

At norb=6 only 7820 / 40320 pphh calls are fully right; linked wrong in 28532, unlinked in 26648.
The linked term is exactly 0 in 12488 calls where it should not be; other ratios are scattered
(1/2, 1/12, 1/4, ~0.204, ~0.408, 1/3, 2, ...). The unlinked term is wrong by 2, -1, -1/2, -2, 2/3, ...
This is the code path (CGC_2e_X with qt_pp/offb=-2, CGC_1e pphh branch with twoSp1, the fabricated
`twoSk[nspin-1] = 1` pad) that read all 14 undefined intermediates. It looks like it needs
re-deriving rather than patching.

Cost of the local RI for pphh vs the (wrong) graphical pphh: +10-25% pspace wall time (timing.py).

## Broader regression check (567d5a9)

`pytest pyscf/csf_fci/test` (as committed) plus every other forge test that uses csf_solver
(siso/test_siso_modelspace, prop/dip_moment/test_{lpdft,mcpdft,cmspdft,sapdft}_pdm,
prop/trans_dip_moment/test_cmspdft_tdm, mcdcft/test_energy_h2), OMP=6:

* csf_mem (ce1d5c8 + its build): 30 failed, 49 passed -- including 6 subtests of
  test_lpdft_pdm.py (h2o L-PDFT dipoles, ftLDA and ftPBE), i.e. the 2e pspace bugs reach real
  MC-PDFT property calculations.
* csf_mem_reads: 49 passed, 6 skipped, 935 subtests, 0 failed.
