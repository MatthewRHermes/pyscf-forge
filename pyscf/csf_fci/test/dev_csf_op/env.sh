#!/bin/bash
N=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)
export PYSCF_EXT_PATH=$N:$PYSCF_EXT_PATH
export OMP_NUM_THREADS=${DEV_OMP:-1}
exec "$@"
