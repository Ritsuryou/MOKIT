
'''
automr utilities and APIs for automr workflow
'''

def find_antibonding_orb(mol, mo, i1=0, i2=0, i3=0, start_from_one=False,
                         ao_ovlp=None, popm='lowdin'):
    '''
    Construct antibonding orbitals for a set of bonding orbitals.
    mol: PySCF molecule object
    mo: molecular orbital coefficients
    k1, k2, k3: orbital indices
    start_from_one: whether the given orbital indices starts from 1 (Fortran
     convention). start_from_one=False means starting from 0 (Python convention).
    ao_ovlp: AO overlap integral matrix. If `None` is given, it will be calculated
     using mol.intor_symmetric('int1e_ovlp') below, otherwise it will be directly
     used.
    popm: population method, 'mulliken' or 'lowdin'
    Note:
    1) mo(:,i3:nif) will be updated, where mo(:,i3:i3+i2-i1) are generated anti-
     bonding orbitals, and mo(:,i3+i2-i1+1:nif) are remaining virtual orbitals.
    2) all MOs are still orthonormalized after calling this function.
    '''
    from mokit.lib.gaussian import get_bfirst_from_mol, get_ao_dip
    from mokit.lib.ortho import check_orthonormal
    from mokit.lib.wfn_analysis import find_antibonding_orbitals
    from mokit.lib.rwwfn import calc_diag_gross_pop, get_mo_center_from_pop

    if start_from_one is True:
        k1 = i1
        k2 = i2
        k3 = i3
    else:
        k1 = i1+1
        k2 = i2+1
        k3 = i3+1
    npair = k2 - k1 + 1
    nbf = mo.shape[0]
    nif = mo.shape[1]

    if k2<k1 or k3<=k2 or nif<k3+npair-1:
        print('i1= %d, i2= %d, i3= %d, nif= %d' %(i1, i2, i3, nif))
        raise ValueError('Wrong orbital indices.')

    natom = mol.natm
    bfirst = get_bfirst_from_mol(mol)

    if ao_ovlp is None:
        S = mol.intor_symmetric('int1e_ovlp')
    else:
        S = ao_ovlp
    pop = calc_diag_gross_pop(natom, nbf, npair, bfirst, S, mo[:,k1-1:k2], popm)
    mo_center = get_mo_center_from_pop(natom, npair, pop)
    center, ao_dip = get_ao_dip(mol, fix_center=True)
    new_mo = find_antibonding_orbitals(k1, k2, k3, natom, nbf, nif, bfirst,
                                       mo_center, S, ao_dip, mo)
    check_orthonormal(nbf, nif, new_mo, S)
    return new_mo


def find_antibonding_orb_in_fch(fchname, i1, i2, i3, start_from_one=False, popm='lowdin'):
    '''
    Construct antibonding orbitals for a set of bonding orbitals. The original MOs
    are stored in fchname.
    i1, i2, i3: orbital indices
    start_from_one: whether the given orbital indices starts from 1 (Fortran
     convention). start_from_one=False means starting from 0 (Python convention).
    popm: population method, 'mulliken' or 'lowdin'
    Note:
    1) mo(:,i3:nif) will be updated, where mo(:,i3:i3+i2-i1) are generated anti-
     bonding orbitals, and mo(:,i3+i2-i1+1:nif) are remaining virtual orbitals.
    2) all MOs are still orthonormalized after calling this function.
    '''
    from mokit.lib.fch2py import fch2py
    from mokit.lib.py2fch import py2fch
    from mokit.lib.rwwfn import load_mol_from_fch, read_nbf_and_nif_from_fch, read_eigenvalues_from_fch
    import shutil

    mol = load_mol_from_fch(fchname)
    nbf, nif = read_nbf_and_nif_from_fch(fchname)
    mo = fch2py(fchname, nbf, nif, 'a')
    new_mo = find_antibonding_orb(mol, mo, i1, i2, i3, start_from_one, None, popm)
    new_fch = fchname[0:fchname.rindex('.fch')]+'_a.fch'
    shutil.copyfile(fchname, new_fch)
    ev = read_eigenvalues_from_fch(fchname, nif, 'a')
    py2fch(new_fch, nbf, nif, new_mo, 'a', ev, False, False)


# update PySCF CASCI solver
def set_fcisolver(fcisolver, mol, mem, mult=None, csf=False, HardWFN=False,
                  CrazyWFN=False):
    if csf: # CSF-based CASCI solver
        from pyscf.csf_fci import csf_solver
        if mult is None:
            mult = mol.spin + 1
        fcisolver = csf_solver(mol, smult=mult)
        fcisolver.max_memory = mem # in MB
        fcisolver.conv_tol = 1e-8
        if CrazyWFN:
            fcisolver.pspace_size = 2000
            fcisolver.max_cycle = 800
        elif HardWFN:
            fcisolver.pspace_size = 1000
            fcisolver.max_cycle = 400
        else:
            fcisolver.pspace_size = 400
            fcisolver.max_cycle = 200
    else:   # determinant-based CASCI solver
        fcisolver.max_memory = mem # in MB
        if CrazyWFN:
            fcisolver.level_shift = 0.2
            fcisolver.pspace_size = 2800
            fcisolver.max_space = 100
            fcisolver.max_cycle = 1000
        elif HardWFN:
            fcisolver.pspace_size = 1400
            fcisolver.max_cycle = 500
        else:
            fcisolver.pspace_size = 500
            fcisolver.max_cycle = 250
    return fcisolver


# CASCI wrapper which is applicable to the ground state, or any one excited
# state, or multiple excited states
def casci_wrapper(mf, nacto, nacte, iroot=0, nstates=0, mult=None, natorb=True,
                  force_csf=False, force_fix_spin=False, HardWFN=False,
                  CrazyWFN=False, ci0=None):
    '''
    We wish to flexibly switch between determinant-based and CSF-based CASCI
    solver, so we design this wrapper for PySCF CASCI. Note: DF (density fitting)
    is considered in `mf` object, we do not need to check or specify DF here.
    iroot=0, nstates=0, the ground state
    iroot=1, nstates=0, the 1st excited state
    iroot=2, nstates=0, the 2nd excited state
    iroot=0, nstates=1, the lowest two states
    iroot=0, nstates=2, the lowest three states
    '''
    import math
    from pyscf import mcscf

    if nacto < 0:
        raise ValueError('nacto is supposed to be a non-negative integer.')
    if iroot < 0:
        raise ValueError('iroot is supposed to be a non-negative integer.')
    if mult is None:
        mult = mf.mol.spin + 1
    else:
        if mult < 1:
            raise ValueError('mult >= 1 is required.')
    if (iroot > 0) and (ci0 is not None):
        raise ValueError('`iroot > 0` and `ci0 != None` are incompatible!')
    if (iroot > 0) and (nstates > 0):
        raise ValueError('`iroot > 0` and `nstates > 0` are incompatible!')
    target_ss = 0.25*(mult*mult - 1)

    # mf.max_memory is in unit MB
    mem1 = math.floor(mf.max_memory*0.7 + 0.5)
    mem2 = math.floor(mf.max_memory*0.3 + 0.5)

    mc = mcscf.CASCI(mf, nacto, nacte)
    mc.max_memory = mem1
    mc.fcisolver = set_fcisolver(mc.fcisolver, mf.mol, mem2, mult, False, HardWFN,
                                 CrazyWFN)
    if iroot > 0:
        # modify mc.fcisolver.nroots here does not make any effect
        mc = mc.state_specific_(iroot)
    if nstates > 0:
        mc.fcisolver.nroots = nstates + 1
    mc.natorb = natorb
    mc.verbose = 5
    if ci0 is None:
        mc.kernel()
    else:
        mc.kernel(ci0=ci0)

    ssquare_diff = 1e-4
    try_csf = False
    if nstates == 0:
        ss = mc.fcisolver.spin_square(mc.ci, mc.ncas, mc.nelecas)
        if abs(ss[0] - target_ss) > ssquare_diff:
            try_csf = True
    else:
        for i in range(nstates):
            ss = mc.fcisolver.spin_square(mc.ci[i], mc.ncas, mc.nelecas)
            if abs(ss[0] - target_ss) > ssquare_diff:
                try_csf = True
                break

    if try_csf:
        print('Remark from MOKIT casci_wrapper: determinant-based CASCI solver does not lead')
        print('to the expected spin state.')
        use_csf = False
        use_fix_spin = False

        if force_fix_spin:
            if force_csf:
                raise ValueError('force_csf and force_fix_spin cannot both be True.')
            else:
                use_fix_spin = True
        else:
            try:
                from pyscf.csf_fci import csf_solver
                use_csf = True
                print('pyscf-forge detected.')
            except ImportError as err:
                if force_csf:
                    raise ImportError(
                        'pyscf.csf_fci is not available. Install the pyscf-forge package'
                        ' or set force_fix_spin=True'
                    ) from err
                else:
                    use_fix_spin = True
                    print('pyscf-forge not detected. Switch to fix_spin.')

        if use_csf:
            print('Using CSF-based CASCI solver pyscf.csf_fci')
            mc.fcisolver = set_fcisolver(mc.fcisolver, mf.mol, mem2, mult, True,
                                         HardWFN, CrazyWFN)
            # state_specific_ is set again since it is deactivated when coming here
            if iroot > 0:
                # modify mc.fcisolver.nroots here does not make any effect
                mc = mc.state_specific_(iroot)
            if nstates > 0:
                mc.fcisolver.nroots = nstates + 1
        if use_fix_spin:
            print('Using PySCF built-in mc.fix_spin_()')
            mc.fix_spin_(ss = target_ss)
        if natorb: # reset MO coefficients
            mc.mo_coeff = mf.mo_coeff.copy()
        mc.kernel()
    return mc


# state-specific CASSCF wrapper which is applicable to ground/excited state
# not applicable to multiple excited states
def casscf_wrapper(mf, nacto, nacte, iroot=0, mult=None, natorb=True, force_csf=False,
                   force_fix_spin=False, HardWFN=False, CrazyWFN=False):
    '''
    We wish to flexibly switch between determinant-based and CSF-based CASCI
    solver, so we design this wrapper for PySCF CASSCF.
    '''
    import math
    import numpy as np
    from pyscf import mcscf

    if mult is None:
        mult = mf.mol.spin + 1
    target_ss = 0.25*(mult*mult - 1)

    if CrazyWFN:
        mc = casci_wrapper(mf, nacto, nacte, 0, iroot+2, mult, False, force_csf,
                           force_fix_spin, False, True, None)
    elif HardWFN:
        mc = casci_wrapper(mf, nacto, nacte, 0, iroot+1, mult, False, force_csf,
                           force_fix_spin, True, False, None)
    else:
        mc = casci_wrapper(mf, nacto, nacte, 0, iroot, mult, False, force_csf,
                           force_fix_spin, False, False, None)

    if (HardWFN is False) and (CrazyWFN is False) and (iroot==0):
        ci0 = mc.ci
    else:
        ci0 = mc.ci[:iroot+1]

    # mf.max_memory is in unit MB
    mem1 = math.floor(mf.max_memory*0.7 + 0.5)
    mem2 = math.floor(mf.max_memory*0.3 + 0.5)

    mc = mcscf.CASSCF(mf, nacto, nacte)
    mc.max_memory = mem1
    mc.fcisolver = set_fcisolver(mc.fcisolver, mf.mol, mem2, mult, False, HardWFN,
                                 CrazyWFN)
    mc.max_cycle = 128
    if iroot > 0:
        mc = mc.state_specific_(iroot)
    mc.natorb = natorb
    mc.verbose = 5
    mc.kernel(ci0=ci0)
    if not mc.converged:
        raise OSError('CASSCF is not converged.')

    ssquare_diff = 1e-4
    ss = mc.fcisolver.spin_square(mc.ci, mc.ncas, mc.nelecas)
    if abs(ss[0] - target_ss) > ssquare_diff:
        print('Remark from MOKIT casscf_wrapper: determinant-based CASCI solver does not')
        print('lead to the expected spin state.')
        use_csf = False
        use_fix_spin = False

        if force_fix_spin:
            if force_csf:
                raise ValueError('force_csf and force_fix_spin cannot both be True.')
            else:
                use_fix_spin = True
        else:
            try:
                from pyscf.csf_fci import csf_solver
                use_csf = True
                print('pyscf-forge detected.')
            except ImportError as err:
                if force_csf:
                    raise ImportError(
                        'pyscf.csf_fci is not available. Install the pyscf-forge package'
                        ' or set force_fix_spin=True'
                    ) from err
                else:
                    use_fix_spin = True
                    print('pyscf-forge not detected. Switch to fix_spin.')

        if use_csf:
            print('Using CSF-based CASCI solver pyscf.csf_fci')
            mc.fcisolver = set_fcisolver(mc.fcisolver, mf.mol, mem2, mult, True,
                                         HardWFN, CrazyWFN)
            if iroot > 0:
                mc = mc.state_specific_(iroot)
        if use_fix_spin:
            print('Using PySCF built-in mc.fix_spin_()')
            mc.fix_spin_(ss = target_ss)
        # reset MO coefficients
        mc.mo_coeff = mf.mo_coeff.copy()
        mc.kernel()
    return mc

