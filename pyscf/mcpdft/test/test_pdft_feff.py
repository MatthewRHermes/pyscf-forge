#!/usr/bin/env python
# Copyright 2014-2023 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Matthew Hennefarth <mhennefarth@uchicago.edu>

import numpy as np
from pyscf import gto, scf, fci, ao2mo, lib, mcscf
from pyscf.lib import temporary_env
from pyscf.mcscf import mc_ao2mo, newton_casscf
from pyscf.mcpdft import pdft_feff, _dms
from pyscf import mcpdft
import unittest


def setUpModule():
    global h2, lih
    h2 = scf.RHF(gto.M(atom='H 0 0 0; H 1.2 0 0', basis='6-31g',
                       output='/dev/null', verbose=0)).run()
    lih = scf.RHF(gto.M(atom='Li 0 0 0; H 1.2 0 0', basis='sto-3g',
                        output='/dev/null', verbose=0)).run()


def tearDownModule():
    global h2, lih
    h2.mol.stdout.close()
    lih.mol.stdout.close()
    del h2, lih


def get_feff_ref(mc, state=0, c_casdm1s=None, c_casdm2=None):
    nao, nmo = mc.mo_coeff.shape
    casdm1s = mc.make_one_casdm1s(mc.ci, state=state)
    casdm2 = mc.make_one_casdm2(mc.ci, state=state)
    dm1s = _dms.casdm1s_to_dm1s(mc, casdm1s)
    cascm2 = _dms.dm2_cumulant(casdm2, casdm1s)
    mo_cas = mc.mo_coeff[:, mc.ncore:][:, :mc.ncas]
    if c_casdm1s is None:
        c_dm1s = dm1s

    else:
        c_dm1s = _dms.casdm1s_to_dm1s(mc, c_casdm1s)

    if c_casdm2 is None:
        c_cascm2 = cascm2

    else:
        c_cascm2 = _dms.dm2_cumulant(c_casdm2, c_casdm1s)

    v1, v2_ao = pdft_feff.lazy_kernel(mc.otfnal, dm1s, cascm2, c_dm1s, c_cascm2, mo_cas)
    with temporary_env(mc._scf, _eri=ao2mo.restore(4, v2_ao, nao)):
        with temporary_env(mc.mol, incore_anyway=True):
            v2 = mc_ao2mo._ERIS(mc, mc.mo_coeff, method='incore')

    return v1, v2


def case(kv, mc):
    ncore, ncas, nelecas = mc.ncore, mc.ncas, mc.nelecas
    nao, nmo = mc.mo_coeff.shape
    nocc, nvir = ncore + ncas, nmo - ncore - ncas
    ngorb = ncore * ncas + nocc * nvir
    fcasscf = mcscf.CASSCF(mc._scf, ncas, nelecas)
    fcasscf.__dict__.update(mc.__dict__)

    feff1, feff2 = mc.get_pdft_veff(mc.mo_coeff, mc.ci, incl_coul=False, paaa_only=True)
    veff1, veff2 = mc.get_pdft_veff(mc.mo_coeff, mc.ci, incl_coul=False, paaa_only=True)

    with lib.temporary_env(fcasscf, get_hcore=lambda: feff1):
        g_all, _, _, hdiag_all = newton_casscf.gen_g_hop(fcasscf, mc.mo_coeff, mc.ci, feff2)

    g_numzero = np.abs(g_all) < 1e-8
    hdiag_all[g_numzero] = 1
    x0 = -g_all / hdiag_all
    xorb_norm = np.linalg.norm(x0[:ngorb])
    xci_norm = np.linalg.norm(x0[ngorb:])
    x0 = g_all * np.random.rand(*x0.shape) - 0.5
    x0[g_numzero] = 0
    x0[:ngorb] *= xorb_norm / np.linalg.norm(x0[:ngorb])
    x0[ngorb:] *= xci_norm / (np.linalg.norm(x0[ngorb:]) or 1)
    err_tab = np.zeros((0, 2))

    def seminum(x):
        uorb, ci1 = newton_casscf.extract_rotation(fcasscf, x, 1, mc.ci)
        mo1 = mc.rotate_mo(mc.mo_coeff, uorb)
        veff1_1, veff2_1 = mc.get_pdft_veff(mo=mo1, ci=ci1, incl_coul=False, paaa_only=True)
        return veff1 - veff1_1, veff2.papa - veff2_1.papa

    for ix, p in enumerate(range(20,25)):
        x1 = x0/(2**p)
        print(x1.shape)
        x1_norm = np.linalg.norm(x1)
        dg_test = np.dot(g_all, x1)
        dveff1, dveff2 = seminum(x1)
        dg_ref = np.concatenate(dveff1).sum() + np.concatenate(dveff2).sum()
        dg_err = abs((dg_test - dg_ref)/dg_ref)
        print("________")
        print(dg_test)
        print(dg_ref)
        print(dg_err)


class KnownValues(unittest.TestCase):

    def test_dvot(self):
        np.random.seed(1)
        for mol, mf in zip(("H2", "LiH"), (h2, lih)):
            for state, nel in zip(('Singlet', 'Triplet'), (2, (2, 0))):
                for fnal in ('tLDA,VWN3', 'ftLDA,VWN3', 'tPBE', 'ftPBE'):
                    mc = mcpdft.CASSCF(mf, fnal, 2, nel, grids_level=1).run()
                    with self.subTest(mol=mol, state=state, fnal=fnal):
                        case(self, mc)

    def test_feff_ao2mo(self):
        for mol, mf in zip(("H2", "LiH"), (h2, lih)):
            for state, nel in zip(('Singlet', 'Triplet'), (2, (2, 0))):
                for fnal in ('tLDA,VWN3', 'ftLDA,VWN3', 'tPBE', 'ftPBE'):
                    mc = mcpdft.CASSCF(mf, fnal, 2, nel, grids_level=1).run()
                    f1_test, f2_test = mc.get_pdft_feff(jk_pc=True)
                    f1_ref, f2_ref = get_feff_ref(mc)
                    f_test = [f1_test, f2_test.vhf_c, f2_test.papa,
                              f2_test.ppaa, f2_test.j_pc, f2_test.k_pc]
                    f_ref = [f1_ref, f2_ref.vhf_c, f2_ref.papa, f2_ref.ppaa,
                             f2_ref.j_pc, f2_ref.k_pc]
                    terms = ['f1', 'f2.vhf_c', 'f2.papa', 'f2.ppaa', 'f2.j_pc',
                             'f2.k_pc']
                    for test, ref, term in zip(f_test, f_ref, terms):
                        with self.subTest(mol=mol, state=state, fnal=fnal,
                                          term=term):
                            self.assertAlmostEqual(lib.fp(test),
                                                   lib.fp(ref), delta=1e-4)

    def test_sa_contract_feff_ao2mo(self):
        for mol, mf in zip(("H2", "LiH"), (h2, lih)):
            for state, nel in zip(['Singlet'], [2]):
                for fnal in ('tLDA,VWN3', 'ftLDA,VWN3', 'tPBE', 'ftPBE'):
                    mc = mcpdft.CASSCF(mf, fnal, 2, nel, grids_level=1).state_average_([0.5, 0.5]).run()

                    sa_casdm1s = _dms.make_weighted_casdm1s(mc)
                    sa_casdm2 = _dms.make_weighted_casdm2(mc)

                    f1_test, f2_test = mc.get_pdft_feff(jk_pc=True, c_casdm1s=sa_casdm1s, c_casdm2=sa_casdm2)
                    f1_ref, f2_ref = get_feff_ref(mc, c_casdm1s=sa_casdm1s, c_casdm2=sa_casdm2)
                    f_test = [f1_test, f2_test.vhf_c, f2_test.papa,
                              f2_test.ppaa, f2_test.j_pc, f2_test.k_pc]
                    f_ref = [f1_ref, f2_ref.vhf_c, f2_ref.papa, f2_ref.ppaa,
                             f2_ref.j_pc, f2_ref.k_pc]
                    terms = ['f1', 'f2.vhf_c', 'f2.papa', 'f2.ppaa', 'f2.j_pc',
                             'f2.k_pc']
                    for test, ref, term in zip(f_test, f_ref, terms):
                        with self.subTest(mol=mol, state=state, fnal=fnal,
                                          term=term):
                            self.assertAlmostEqual(lib.fp(test),
                                                   lib.fp(ref), delta=1e-4)


if __name__ == "__main__":
    print("Full Tests for pdft_feff")
    unittest.main()