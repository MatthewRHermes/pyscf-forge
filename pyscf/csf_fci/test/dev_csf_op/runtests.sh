#!/bin/bash
# runtests.sh : run the must-stay-passing csf tests against this worktree's libcsf.so.
# test_csf_4elec is run (separately) with every @unittest.skip removed; only its smult<=3 subtests may fail.
H=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); T=$H/..
TMP=$(mktemp -d /tmp/claude-1000/runtests.XXXX)
sed "s/^\(\s*\)@unittest.skip('debug')/\1#SKIPREMOVED/" $T/test_csf_4elec.py > $TMP/test_csf_4elec_all.py
PT="python -m pytest -p no:warnings -p no:cacheprovider -q -rA"
DEV_OMP=${DEV_OMP:-6} conda run -n pyscf $H/env.sh $PT $T/test_csf_2elec.py $T/test_csf_spinless.py $T/test_csf_1eop.py > $TMP/a.txt 2>&1
(cd $TMP && DEV_OMP=${DEV_OMP:-6} conda run -n pyscf $H/env.sh $PT test_csf_4elec_all.py > $TMP/b.txt 2>&1)
fa=$(grep -E '^(SUBFAILED|FAILED|ERROR)' $TMP/a.txt | grep -v '^ERROR conda' | sed 's/ - .*//')
fb=$(grep -E '^(SUBFAILED|FAILED|ERROR)' $TMP/b.txt | grep -v '^ERROR conda' | sed 's/ - .*//' | grep -vE '^FAILED test_csf_4elec_all.py::KnownValues::test_(kernel|pspace)$' | grep -vE '^SUBFAILED\(smult=[123](, m=-?[0-9]+)?\) test_csf_4elec_all.py::KnownValues::test_(kernel|pspace)$')
echo "2elec/spinless/1eop: $(grep -E '[0-9]+ (passed|failed)' $TMP/a.txt | tail -1)"
echo "4elec (all):         $(grep -E '[0-9]+ (passed|failed)' $TMP/b.txt | tail -1)"
echo "4elec smult<=3 subfailures (allowed): $(grep -cE '^SUBFAILED\(smult=[123]' $TMP/b.txt)"
if [ -n "$(echo "$fa$fb" | tr -d "[:space:]")" ]; then echo "UNEXPECTED FAILURES:"; echo "$fa"; echo "$fb"; else echo "OK: no unexpected failures"; fi
cat $TMP/a.txt $TMP/b.txt | grep -E 'Assertion `|Segmentation|Fatal Python' | head -3
[ -n "$KEEP" ] && echo "kept $TMP" || rm -rf $TMP
