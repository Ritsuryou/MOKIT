#!/usr/bin/env python

# This file is firstly created at 20240701.
# The reasons to create this file:
# 1) PySCF use tol=1e-4 by default in the Davidson step for wave function stability
#    check, which is insufficient for some difficult systems (1e-7 required). The
#    tol value cannot be modified in the master branch currently. So jxzou copies
#    and modifies some content from pyscf/scf/stability.py to this file.
# 2) The `stable=opt` functionality is moved from mokit/src/do_hf.f90 to this file.

# This file is for temporarily usage, which may be removed in ~2 years (when most
# people are used to new PySCF, say >= v2.7). Forcing all users to update to the
# newest PySCF is not realistic currently.

# Related pull request (https://github.com/pyscf/pyscf/pull/2299).
# If any PySCF developers think this action is inappropriate, please contact
# MOKIT developers to delete this file.

from pyscf import scf, lib
from pyscf.soscf import newton_ah
import numpy as np

MAX_CYCLE = 300

def _rotate_mo(mo_coeff, mo_occ, dx):
    dr = scf.hf.unpack_uniq_var(dx, mo_occ)
    u = newton_ah.expmat(dr)
    return np.dot(mo_coeff, u)


# This function can be used for either RHF or ROHF
def rhf_internal(mf, tol=1e-7, verbose=None):
    log = lib.logger.new_logger(mf, verbose)
    if mf.mol.spin == 0:
        g, hop, hdiag = newton_ah.gen_g_hop_rhf(mf, mf.mo_coeff, mf.mo_occ,
                                                with_symmetry=True)
    else:
        g, hop, hdiag = newton_ah.gen_g_hop_rohf(mf, mf.mo_coeff, mf.mo_occ,
                                                 with_symmetry=True)
    hdiag *= 2
    stable = True

    def precond(dx, e, x0):
        hdiagd = hdiag - e
        hdiagd[abs(hdiagd)<1e-8] = 1e-8
        return dx/hdiagd

    def hessian_x(x):
        return hop(x).real * 2

    x0 = np.zeros_like(g)
    x0[g!=0] = 1. / hdiag[g!=0]
    e, v = lib.davidson(hessian_x, x0, precond, tol=tol, max_cycle=MAX_CYCLE, verbose=log)
    if e < -1e-5:
        log.note(f'{mf.__class__} wavefunction has an internal instability.')
        mo = _rotate_mo(mf.mo_coeff, mf.mo_occ, v)
        stable = False
    else:
        log.note(f'{mf.__class__} wavefunction is stable in the internal '
                 'stability analysis')
        mo = mf.mo_coeff
    return mo, stable


def uhf_internal(mf, tol=1e-7, verbose=None):
    log = lib.logger.new_logger(mf, verbose)
    g, hop, hdiag = newton_ah.gen_g_hop_uhf(mf, mf.mo_coeff, mf.mo_occ,
                                            with_symmetry=True)
    hdiag *= 2
    stable = True

    def precond(dx, e, x0):
        hdiagd = hdiag - e
        hdiagd[abs(hdiagd)<1e-8] = 1e-8
        return dx/hdiagd

    def hessian_x(x):
        return hop(x).real * 2

    x0 = np.zeros_like(g)
    x0[g!=0] = 1. / hdiag[g!=0]
    e, v = lib.davidson(hessian_x, x0, precond, tol=tol, max_cycle=MAX_CYCLE, verbose=log)
    if e < -1e-5:
        log.note(f'{mf.__class__} wavefunction has an internal instability.')
        nocca = np.count_nonzero(mf.mo_occ[0]> 0)
        nvira = np.count_nonzero(mf.mo_occ[0]==0)
        mo = (_rotate_mo(mf.mo_coeff[0], mf.mo_occ[0], v[:nocca*nvira]),
              _rotate_mo(mf.mo_coeff[1], mf.mo_occ[1], v[nocca*nvira:]))
        stable = False
    else:
        log.note(f'{mf.__class__} wavefunction is stable in the internal '
                 'stability analysis')
        mo = mf.mo_coeff
    return mo, stable


def loop_soscf(mf):
    from mokit.lib.rwwfn import get_occ_from_na_nb, get_occ_from_na_nb2
    na, nb = mf.mol.nelec
    uhf = isinstance(mf, scf.uhf.UHF)
    if uhf:
        nif = mf.mo_coeff[0].shape[1]
    else:
        nif = mf.mo_coeff.shape[1]
    old_cyc = mf.max_cycle
    mf.max_cycle = 64

    def is_descending_np(arr):
        return np.all(arr[:-1] >= arr[1:])

    for i in range(10):
        run_kernel = not mf.converged
        if mf.converged and not is_descending_np(mf.mo_occ):
            if uhf:
                occ = get_occ_from_na_nb2(nif, na, nb)
            else:
                occ = get_occ_from_na_nb(nif, na, nb)
            mf.mo_occ = occ.copy()
            run_kernel = True
        if run_kernel:
            if not isinstance(mf, newton_ah._CIAH_SOSCF):
                mf = mf.newton()
            mf.kernel()
        else:
            break
    else:
        raise OSError('PySCF SOSCF failed after 10 attempts.')
    mf.max_cycle = old_cyc
    return mf


def hf_stable_opt_internal(mf):
    # in case that the input object `mf` is unconverged, let's check it first
    if not mf.converged:
        mf = loop_soscf(mf)

    rhf = isinstance(mf, scf.rhf.RHF) or isinstance(mf, scf.rohf.ROHF)
    uhf = isinstance(mf, scf.uhf.UHF)

    i = 0
    while (i < 10):
        i += 1
        if rhf:
            mo, stable = rhf_internal(mf, verbose=5)
        elif uhf:
            mo, stable = uhf_internal(mf, verbose=5)
        else:
            raise OSError('Only RHF/ROHF/UHF are supported currently.')

        if (stable):
            break
        else:
            mf.converged = False
            mf.mo_coeff = mo
            mf = loop_soscf(mf)
    if not stable:
        raise OSError('PySCF R(O)HF stable=opt failed after 10 attempts.')
    return mf

# alias
rhf_stable_opt_internal = hf_stable_opt_internal
rohf_stable_opt_internal = hf_stable_opt_internal
uhf_stable_opt_internal = hf_stable_opt_internal

