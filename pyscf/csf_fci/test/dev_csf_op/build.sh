#!/bin/bash
# build.sh [srcdir] [out.so] [extra gcc flags...]
N=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)
SRC=${1:-$N/pyscf/lib/csf}; OUT=${2:-$N/pyscf/lib/libcsf.so}; shift 2 2>/dev/null
PL=/home/herme068/anaconda3/envs/pyscf/lib/python3.12/site-packages/pyscf/lib
gcc -Wall -msse3 -O1 -g -fPIC -fopenmp -shared -I$SRC -I$PL -I$PL/deps/include "$@" $SRC/*.c -o $OUT -lm 2>&1 | grep -E 'error|warning: unused variable' | grep error
exit ${PIPESTATUS[0]}
