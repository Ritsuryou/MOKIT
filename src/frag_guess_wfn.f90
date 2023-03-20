! written by jxzou at 20201121: perform fragment-guess wavefunction calculations
! updated by jxzou at 20210409: add support of MOROKUMA-EDA and GKS-EDA
! updated by jxzou at 20210414: automatically add $DFTTYP and $PCM into GAMESS .inp file
! updated by jxzou at 20210417: construct supermolecule MOs from direct sum of fragments MOs

! The fragment-guess wavefunction calculation in Gaussian has the following
! shortcomings:
! (1) cannot ensure the wavefunction stability of each fragment. This often
!     makes the calculations of transition-metal molecules take more SCF cycles,
!     since the initial guess is not as good as we expect.
! (2) do not support mixed methods, e.g. cannot use RHF/RDFT for some fragments
!     and UHF/UDFT for other fragments
! This utility aims to overcome these shortcomings.

module frag_info
 implicit none
 integer :: nfrag0 ! No. of fragments in .gjf file
 integer :: nfrag  ! No. of fragments (including total system) to be calculated
 ! For none-EDA and Morokuma-EDA calculations, nfrag = nfrag0
 ! For GKS-EDA calculations, nfrag = 2*nfrag0 + 1

 type :: frag
  integer :: charge = 0
  integer :: mult = 1
  integer :: natom = 0
  integer :: wfn_type = 0            ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
  integer, allocatable :: atm_map(:) ! map to parent system, i.e. supermolecule
  real(kind=8) :: e = 0d0            ! electronic energy
  real(kind=8) :: ssquare = 0d0      ! S(S+1)
  real(kind=8), allocatable :: coor(:,:)   ! size natom
  character(len=2), allocatable :: elem(:) ! size natom
  character(len=240) :: fname = ' '        ! gjf filename
  logical :: pos = .true.            ! .True./.False. for alpha/beta total spin
  logical :: noiter = .false.        ! .True./.False. for skipping SCF or not
  logical, allocatable :: ghost(:)   ! .True./.False. for ghost atoms or not, size natom
 end type frag

 type(frag), allocatable :: frags(:)
contains

! copy type frag
subroutine copy_frag(frag1, frag2)
 implicit none
 integer :: natom
 type(frag), intent(in) :: frag1
 type(frag), intent(out) :: frag2

 frag2%charge   = frag1%charge
 frag2%mult     = frag1%mult
 frag2%wfn_type = frag1%wfn_type
 frag2%e        = frag1%e
 frag2%ssquare  = frag1%ssquare
 frag2%fname    = frag1%fname
 frag2%pos      = frag1%pos
 frag2%noiter   = frag1%noiter
 natom = frag1%natom
 frag2%natom = natom
 if(natom > 0) then
  allocate(frag2%coor(3,natom), source=frag1%coor)
  allocate(frag2%elem(natom), source=frag1%elem)
  if(allocated(frag1%atm_map)) then
   allocate(frag2%atm_map(natom), source=frag1%atm_map)
  end if
  if(allocated(frag1%ghost)) then
   allocate(frag2%ghost(natom), source=frag1%ghost)
  end if
 end if
end subroutine copy_frag
end module frag_info

module theory_level
 implicit none
 integer :: mem, nproc ! mem is in MB
 integer :: wfn_type = 0  ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 integer :: eda_type = 0  ! 0/1/2/3/4 for none/Morokuma-EDA/LMO-EDA/GKS-EDA/SAPT
 integer :: disp_type = 0 ! 0/1/2 for no dispersion correction/GD3/GD3BJ
 character(len=9) :: method = ' '
 character(len=21) :: basis = ' '
 character(len=60) :: scrf  = ' '
 character(len=40) :: solvent = ' '
 character(len=40) :: solvent_gau = ' '
 logical :: sph = .true.  ! spherical harmonic/Cartesian functions
end module theory_level

program main
 implicit none
 integer :: i
 character(len=24) :: data_string
 character(len=240) :: gau_path, gjfname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine frag_guess_wfn: wrong command line&
                   & arguments!'
  write(6,'(/,A)') " Example 1 (in bash): frag_guess_wfn h2o_dimer.gjf >& h2o_dimer.out &"
  write(6,'(A,/)') " Example 2 (in dash): frag_guess_wfn h2o_dimer.gjf >h2o_dimer.out 2>&1 &"
  stop
 end if

 gjfname = ' '
 call getarg(1, gjfname)
 if(index(gjfname,'.gjf') == 0) then
  write(6,'(/,A)') "ERROR in subroutine frag_guess_wfn: '.gjf' suffix not found&
                  & in filename "//TRIM(gjfname)
  stop
 end if

 call require_file_exist(gjfname)
 call get_gau_path(gau_path)

 call fdate(data_string)
 write(6,'(A)') TRIM(data_string)

 call frag_guess_wfn(gau_path, gjfname)
 call gen_inp_of_frags()

 call fdate(data_string)
 write(6,'(/,A)') TRIM(data_string)
end program main 

! perform fragment-guess wavefunction calculations, one by one fragment
subroutine frag_guess_wfn(gau_path, gjfname)
 use util_wrapper, only: formchk, unfchk
 use frag_info, only: nfrag0, nfrag, frag, frags, copy_frag
 use theory_level
 use phys_cons, only: au2kcal
 implicit none
 integer :: i, j, k, m, fid, ifrag0, ifrag, iatom
 integer :: charge, mult, natom, maxfrag
 integer, allocatable :: cm(:), nuc(:)
 real(kind=8), allocatable :: coor(:,:), tmp_coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: buf, chkname, fchname, logname, basname
 character(len=240), intent(in) :: gau_path, gjfname
 character(len=1200) :: longbuf
 logical :: guess_read, stab_chk
 type(frag) :: tmp_frag1, tmp_frag2

 buf = ' '; longbuf = ' '

 call read_eda_type_from_gjf(gjfname, eda_type, stab_chk)
 if(.not. stab_chk) then
  write(6,'(/,A)') 'Warning from subroutine frag_guess_wfn: wave function stabili&
                   &ty check is turned'
  write(6,'(A,/)') 'off. Not recommended. You should know what you are doing.'
 end if

 call read_disp_ver_from_gjf(gjfname, disp_type)
 if(eda_type==4 .and. disp_type>0) then
  write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: dispersion correction&
                   & cannot be used in SAPT.'
  write(6,'(A)') 'It may be supported in DFT-SAPT, but this method is not&
                 & supported for frag_guess_wfn.'
  stop
 end if

 call check_sph_in_gjf(gjfname, sph)
 if(sph .and. eda_type==1) then
  write(6,'(A)') 'Warning in subroutine frag_guess_wfn: spherical harmonic&
                 & functions (5D 7F) cannot be'
  write(6,'(A)') 'used in Morokuma-EDA method in GAMESS. Automatically switch&
                 & to Cartesian functions (6D 10F).'
  sph = .false.
 end if

 call read_mem_and_nproc_from_gjf(gjfname, mem, nproc)
 write(6,'(A,I0,A)') '%mem=', mem, 'MB'
 write(6,'(A,I0)') '%nprocshared=', nproc
 call read_natom_from_gjf(gjfname, natom) ! read the total number of atoms
 write(6,'(A,I0)') 'natom=', natom

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: incomplete file '//TRIM(gjfname)
  write(6,'(A)') "Failed to locate the Route Section '#' in .gjf file."
  close(fid)
  stop
 end if

 longbuf = TRIM(buf)
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  longbuf = TRIM(longbuf)//' '//TRIM(buf)
 end do ! for while

 call lower(longbuf)
 call read_scrf_string_from_buf(longbuf, scrf)
 if(LEN_TRIM(scrf)>0 .and. eda_type==4) then
  write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: SCRF detected.'
  write(6,'(A)') 'Implicit solvent model is not supported in SAPT computations.'
  stop
 end if
 call determine_solvent_from_gau2gms(scrf, solvent)
 call read_nfrag_from_buf(longbuf, nfrag0)

 select case(eda_type)
 case(0) ! fragment guess wfn
  nfrag = nfrag0 + 1
  maxfrag = 998
 case(1) ! MOROKUMA-EDA
  nfrag = nfrag0 + 1
  maxfrag = 998
 case(2) ! LMO-EDA
  nfrag = 2*nfrag0 + 1
  maxfrag = 499
 case(3) ! GKS-EDA
  nfrag = 2*nfrag0 + 1
  maxfrag = 499
 case(4) ! SAPT
  if(nfrag0 /= 2) then
   write(6,'(A,I0)') 'ERROR in subroutine frag_guess_wfn: invalid nfrag0=',nfrag0
   write(6,'(A)') 'Only two fragments/monomers are allowed in SAPT.'
   stop
  end if
  nfrag = 3
  maxfrag = 3
 case default
  write(6,'(A,I0)') 'ERROR in subroutine frag_guess_wfn: invalid eda_type=',eda_type
  write(6,'(A,I0)') 'Only 0,1,2,3,4 are allowed.'
  close(fid)
  stop
 end select

 if(nfrag0 > maxfrag) then
  write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: nfrag0 > maxfrag, too large!'
  write(6,'(2(A,I0))') 'nfrag0=', nfrag0, ', maxfrag=', maxfrag
  close(fid)
  stop
 end if

 write(6,'(A,I0)') 'nfrag0=', nfrag0
 write(6,'(A,I0)') 'nfrag=', nfrag

 call read_method_and_basis_from_buf(longbuf, method, basis, wfn_type)
 write(6,'(A,I0)') 'method='//TRIM(method)//', basis='//TRIM(basis)//&
                  &', wfn_type=', wfn_type

! Strange thing: if the following 'if(' lines are uncommented, GKS-EDA will
!  signal unknown errors
! if(TRIM(method) == 'pbe0') then
!  write(6,'(/,A)') "ERROR in subroutine frag_guess_wfn: 'PBE0' string found in&
!                   & gjf file."
!  write(6,'(A)') "Be careful! The PBE0 functional in papers corresponds to the&
!                 & keyword 'PBE1PBE' in Gaussian."
!  write(6,'(A)') 'You should check what you want to calculate.'
!  stop
! end if

 if(index(basis,'gen') > 0) then
  close(fid)
  call record_gen_basis_in_gjf(gjfname, basname, .false.)
  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do ! for while
 else
  basname = ' '
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while
 allocate(cm(2+2*nfrag0),source=0)

 read(fid,*,iostat=i) cm
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: incomplete charges and&
                   & spin multiplicities'
  write(6,'(A)') 'detected in file '//TRIM(gjfname)
  close(fid)
  stop
 end if

 charge = cm(1); mult = cm(2)

 if(mult/=1 .and. eda_type==1) then
  write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: Morokuma-EDA can only&
                 & be applied to RHF.'
  write(6,'(A)') 'But the total spin is not singlet.'
  close(fid)
  stop
 end if

 ! compare the spin and method name, determine RHF/ROHF/UHF-type wfn
 if(mult == 1) then
  if(wfn_type == 0) then
   wfn_type = 1 ! RHF
  else if(wfn_type == 2) then
   write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: total spin is singlet.&
                  & But you specify ROHF/RODFT.'
   close(fid)
   stop
  end if
 else ! mult /= 1
  if(wfn_type == 0) then
   wfn_type = 3 ! UHF
  else
   if(wfn_type == 1) then
    write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: total spin is non-singlet.&
                   & But you specify RHF/RDFT.'
    stop
   end if
  end if
 end if

 if(eda_type==1 .and. wfn_type/=1) then
  write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: Morokuma-EDA can only&
                 & be applied to RHF. But'
  write(6,'(A,I0)') 'RHF-type wavefunction is not specified. wfn_type=',wfn_type
  close(fid)
  stop
 end if

 ! Initially, the restricted/unrestricted type of wfn of each fragment is
 !  inherited from the method of the total system.
 ! This type may be updated when reading Cartesian coordinates
 allocate(frags(nfrag))
 frags(:)%wfn_type = wfn_type

 forall(i = 1:nfrag0)
  frags(i)%charge = cm(2*i+1)
  frags(i)%mult = cm(2*i+2)
 end forall
 deallocate(cm)

 i = SUM(frags(1:nfrag0)%charge)
 if(i /= charge) then
  write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: sum of fragment&
                   & charges is not equal to'
  write(6,'(2(A,I0))') 'total charge. Total charge=', charge, &
                       ', sum(frag_charges)=', i
  write(6,'(A)') 'Wrong charges in file '//TRIM(gjfname)
  deallocate(frags)
  close(fid)
  stop
 end if

 if(ANY(frags(:)%mult==0) .or. ANY(frags(:)%mult==-1)) then
  write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: nonsense chare or spin&
                   & multiplicity'
  write(6,'(A)') 'Please check your file '//TRIM(gjfname)
  deallocate(frags)
  close(fid)
  stop
 end if

 if(mult==1 .and. ANY(frags(:)%mult/=1) .and. wfn_type==1) then
  write(6,'(/,A)') 'Remark: the total spin is singlet but some fragment is open&
                   &-shell.'
  write(6,'(A)') 'Automatically switching to wfn_type=3.'
  wfn_type = 3; frags(:)%wfn_type = 3
 end if

 j = 0 ! ne(alpha) - ne(beta)

 do i = 1, nfrag0, 1
  k = frags(i)%mult
  if(k > 0) then
   j = j + k - 1
  else
   j = j + k + 1
  end if
 end do ! for i

 if(mult-1 /= j) then
  write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: number of unpaired ele&
                   &ctrons calculated from'
  write(6,'(A)') 'fragments spin multiplicities is not equal to that calculate&
                 &d from the total spin'
  write(6,'(A)') 'multiplicity. Some spin multiplicity must be wrong.'
  write(6,'(A)') 'Check your spin in file '//TRIM(gjfname)
  deallocate(frags)
  close(fid)
  stop
 end if

 ! now make all spin of fragment positive, move negative(if any) into pos
 do i = 1, nfrag0, 1
  k = frags(i)%mult
  if(k < 0) then
   frags(i)%mult = -k
   frags(i)%pos = .false.
  end if
 end do ! for i

 allocate(coor(3,natom), source=0d0)
 allocate(elem(natom))
 elem = '  '; ifrag0 = 1

 do i = 1, natom, 1
  read(fid,'(A)',iostat=j) buf
  if(j/=0 .or. LEN_TRIM(buf)==0) exit

  j = index(buf,'(')
  k = index(buf,'=')
  m = index(buf,')')
  if(j*k*m == 0) then
   write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: wrong format in&
                    & file '//TRIM(gjfname)//'.'
   write(6,'(A)') 'Please read examples in Section 5.3.2 ~ 5.3.4 of MOKIT&
                  & manual and learn how to'
   write(6,'(A)') 'write a valid input file.'
   stop
  end if
  read(buf(1:j-1),*) elem(i)
  read(buf(k+1:m-1),*) ifrag

  if((i==1 .and. eda_type/=0 .and. ifrag/=1) .or. ifrag<ifrag0) then
   write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: error definition of fr&
                    &agment number.'
   write(6,'(A)') 'This task requires the definition of fragment number is mon&
                  &omer1, monomer2, monomer3,...'
   close(fid)
   stop
  end if

  if(ifrag > nfrag0) then
   write(6,'(A)') 'ERROR in subroutine frag_guess_wfn: the number of fragments&
                  & detected in Cartesian coordinates is'
   write(6,'(A)') 'larger than that found in charges and spin multiplicities.'
   write(6,'(2(A,I0))') 'ifrag = ', ifrag, ', nfrag0=', nfrag0
   close(fid)
   stop
  end if

  read(buf(m+1:),*) coor(1:3,i)
  frags(ifrag)%natom = frags(ifrag)%natom + 1
  iatom = frags(ifrag)%natom

  if(iatom == 1) then
   allocate(frags(ifrag)%coor(3,1))
   frags(ifrag)%coor(:,1) = coor(:,i)

   frags(ifrag)%elem = [elem(i)]
   frags(ifrag)%atm_map = [i]
  else
   allocate(tmp_coor(3,iatom-1))
   tmp_coor = frags(ifrag)%coor
   deallocate(frags(ifrag)%coor)
   allocate(frags(ifrag)%coor(3,iatom))
   frags(ifrag)%coor(:,1:iatom-1) = tmp_coor
   deallocate(tmp_coor)
   frags(ifrag)%coor(:,iatom) = coor(:,i)

   frags(ifrag)%elem = [frags(ifrag)%elem, elem(i)]
   frags(ifrag)%atm_map =[frags(ifrag)%atm_map, i]
  end if

  j = LEN_TRIM(buf)
  if(buf(j:j)=='r' .or. buf(j:j)=='R') then
   if(frags(ifrag)%mult == 1) then
    frags(ifrag)%wfn_type = 1 ! RHF
   else
    frags(ifrag)%wfn_type = 2 ! ROHF
   end if
  end if

  if(wfn_type==1 .and. frags(ifrag)%wfn_type/=1) then
   close(fid)
   write(6,'(/,A)') 'ERROR in subroutine frag_guess_wfn: conflicted type of wfn&
                    & between the total'
   write(6,'(A)') 'system and some fragment. Restricted closed shell wfn is spe&
                  &cified for the total'
   write(6,'(A)') 'system. But restricted open shell or unrestricted wfn is spe&
                  &cified for some fragment.'
   write(6,'(3(A,I0))') 'natom=', natom, ', i=', i, ', ifrag=', ifrag
   stop
  end if

  ifrag0 = ifrag ! update ifrag0
 end do ! for i

 close(fid)

 do i = 1, nfrag0, 1
  iatom = frags(i)%natom
  if(iatom == 0) then
   write(6,'(/,A,I0,A)') 'ERROR in subroutine frag_guess_wfn: the ',i,'-th&
                         & fragment has 0 atom.'
   write(6,'(A)') 'Please check your specification of fragments in file '//&
                   TRIM(gjfname)
   stop
  end if
  allocate(frags(i)%ghost(iatom))
  frags(i)%ghost = .false.
 end do ! for i

 ! this is the super molecule/total system
 frags(nfrag)%natom = natom
 allocate(frags(nfrag)%elem(natom), frags(nfrag)%coor(3,natom), nuc(natom), &
          frags(nfrag)%atm_map(natom))
 forall(i = 1:natom) frags(nfrag)%atm_map(i) = i
 call read_elem_and_coor_from_gjf(gjfname, natom, frags(nfrag)%elem, nuc, &
  frags(nfrag)%coor, frags(nfrag)%charge, frags(nfrag)%mult)
 deallocate(nuc)
 iatom = frags(nfrag)%natom
 allocate(frags(nfrag)%ghost(iatom))
 frags(nfrag)%ghost = .false.
 frags(nfrag)%noiter = .true.

 select case(eda_type)
 case(2,3) ! LMO-/GKS-EDA
  do i = 1, nfrag0, 1
   call gen_extend_bas_frag(frags(i), frags(nfrag), frags(i+nfrag0))
  end do ! for i
 case(4) ! SAPT
  call gen_extend_bas_frag(frags(1), frags(3), tmp_frag1)
  call gen_extend_bas_frag(frags(2), frags(3), tmp_frag2)
  call copy_frag(tmp_frag1, frags(1)) ! overwrite frags(1)
  call copy_frag(tmp_frag2, frags(2)) ! overwrite frags(2)
 end select

 ! generate these SCF .gjf files
 do i = 1, nfrag, 1
  write(frags(i)%fname,'(I3.3,A)') i, '-'//TRIM(gjfname)
  guess_read = .false.
!  if(eda_type==3 .and. i>nfrag0 .and. i<nfrag) guess_read = .true.
  call gen_gjf_from_type_frag(frags(i), guess_read, stab_chk, basname)
 end do ! for i

 j = index(frags(1)%fname, '.gjf', back=.true.)

 ! do SCF computations one by one
 do i = 1, nfrag, 1
  call do_scf_and_read_e(gau_path, gau_path, frags(i)%fname, frags(i)%noiter, &
                         frags(i)%e, frags(i)%ssquare)
  if(i < nfrag) then
   write(6,'(A,I3,A,F18.9,A,F6.2)') 'i=', i, ', frags(i)%e = ', frags(i)%e,&
                                    ', frags(i)%ssquare=', frags(i)%ssquare
  end if
  chkname = frags(i)%fname(1:j-1)//'.chk'
  fchname = frags(i)%fname(1:j-1)//'.fch'
  buf     = frags(i)%fname(1:j-1)//'.gjf'
  logname = frags(i)%fname(1:j-1)//'.log'
  call delete_file(logname)
  if(i < nfrag) call delete_file(buf)
  if(i==nfrag .and. eda_type==1) then ! MOROKUMA-EDA
   call formchk(chkname, fchname)
   call delete_files(2, [chkname, buf])
  end if
 end do ! for i

 select case(eda_type)
 case(2,3) ! For LMO/GKS-EDA, construct supermolecule MOs from direct sum of
           ! fragment MOs
  i = nfrag
  call direct_sum_frag_mo2super_mo(nfrag0, frags(1:nfrag0)%fname, &
        frags(1:nfrag0)%wfn_type, frags(1:nfrag0)%pos, frags(i)%fname, &
        frags(i)%wfn_type, stab_chk)
  ! or we can use the direct sum of MOs
  !call sum_frag_density_and_prt_into_fch(nfrag0, frags(nfrag0+1:2*nfrag0)%fname,&
  !                            frags(nfrag0+1:2*nfrag0)%pos, frags(nfrag)%fname)
  !k = 1
  !if(frags(nfrag)%wfn_type == 3) k = 0
  !call gen_no_using_density_in_fch(fchname, k)
  !call unfchk(fchname, chkname)
  !call modify_guess_only_in_gjf(frags(nfrag)%fname)
  frags(i)%noiter = .false.
  call do_scf_and_read_e(gau_path, gau_path, frags(i)%fname, frags(i)%noiter,&
                         frags(i)%e, frags(i)%ssquare)
  write(6,'(A,I3,A,F18.9,A,F6.2)') 'i=', i, ', frags(i)%e = ', frags(i)%e, &
                                   ', frags(i)%ssquare=', frags(i)%ssquare
  write(6,'(/,A)') 'If you are performing a GKS-EDA calculation using DFT or a &
                   &LMO-EDA calculation'
  write(6,'(A)') "using HF, then the following energy should be close to 'TOTAL&
                 & INTERACTION ENERGY'"
  write(6,'(A)') "of 'OWN BASIS SET' section in GAMESS output:"
  write(6,'(F10.2,A)') (frags(i)%e-SUM(frags(1:nfrag0)%e))*au2kcal, ' kcal/mol'

  write(6,'(/,A)') "and the following energy should be close to 'TOTAL INTERACT&
                   &ION ENERGY' of"
  write(6,'(A)') "'ALL BASIS SET' section in GAMESS output:"
  write(6,'(F10.2,A)') (frags(i)%e-SUM(frags(nfrag0+1:2*nfrag0)%e))*au2kcal,&
                       ' kcal/mol'

 case(4) ! For SAPT, add SCF density of two fragments to obtain approximate
         ! total SCF density
  call sum_frag_density_and_prt_into_fch(2, frags(1:2)%fname, frags(1:2)%pos, &
                                         frags(3)%fname)
  k = 1
  if(frags(3)%wfn_type == 3) k = 0
  call gen_no_using_density_in_fch(fchname, k)
  call unfchk(fchname, chkname)
  frags(3)%noiter = .false.
  call modify_guess_only_in_gjf(frags(3)%fname)
  call do_scf_and_read_e(gau_path, gau_path, frags(3)%fname, frags(3)%noiter, &
                         frags(3)%e, frags(3)%ssquare)
  write(6,'(A,F18.9,A,F6.2)') 'i=  3, frags(i)%e = ', frags(3)%e, &
                              ', frags(i)%ssquare=', frags(3)%ssquare
 end select
 call delete_file(frags(nfrag)%fname)
! call delete_files(2, [frags(nfrag)%fname, logname])

 if(eda_type == 4) then ! SAPT in PSI4
  write(6,'(A)') 'All SCF computations finished. Generating PSI4 input file...'
 else
  write(6,'(/,A)') 'All SCF computations finished. Generating GAMESS input fil&
                   &e...'
 end if
end subroutine frag_guess_wfn

! read the parameter eda_type from a given .gjf file
subroutine read_eda_type_from_gjf(gjfname, eda_type, stab_chk)
 implicit none
 integer :: i, j, k, len_buf, fid
 integer, intent(out) :: eda_type
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 character(len=69), parameter :: error_str = 'ERROR in subroutine read_eda_type&
  &_from_gjf: wrong syntax in gjf file.'
 logical, intent(out) :: stab_chk ! True/False for 'stable=opt' or not

 eda_type = 0; stab_chk = .true.
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while

 read(fid,'(A)') buf ! Title Card line
 close(fid)
 len_buf = LEN_TRIM(buf)

 i = index(buf, '{')
 j = index(buf, '}')
 if(i==0 .or. j==0) then
  write(6,'(/,A)') error_str
  write(6,'(A)') "No '{' or '}' found in file "//TRIM(gjfname)
  stop
 end if

 buf = ADJUSTL(buf(i+1:j-1))
 k = index(buf,',')
 if(k > 0) j = k
 if(j < i) then
  write(6,'(/,A)') error_str
  write(6,'(A)') 'buf='//TRIM(buf)
  stop
 end if

 call lower(buf)

 select case(buf(1:j-1))
 case('frag')
 case('morokuma')
  eda_type = 1
 case('lmo')
  eda_type = 2
 case('gks')
  eda_type = 3
 case('sapt')
  eda_type = 4
 case default
  write(6,'(/,A)') error_str
  write(6,'(A)') 'You are supposed to write one of {morokuma}, {lmo}, {gks} or&
                 & {sapt,bronze} in'
  write(6,'(A)') 'the Title Card line.'
  stop
 end select

 k = index(buf,',')

 if(k > 0) then
  buf = buf(k+1:)
  select case(TRIM(buf))
  case('bronze') ! SAPT,bronze
  case('nostab')
   stab_chk = .false.
  case default
   write(6,'(/,A)') error_str
   write(6,'(A)') "Invalid keyword '"//TRIM(buf)//"'"
   stop
  end select
 end if

end subroutine read_eda_type_from_gjf

! read memory and nprocshared from a given .gjf file
subroutine read_mem_and_nproc_from_gjf(gjfname, mem, np)
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: mem, np
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 ! default settings
 mem = 1000 ! 1000 MB
 np = 1     ! 1 core

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call lower(buf)
  if(buf(1:1) == '#') exit

  j = LEN_TRIM(buf)
  k = index(buf,'=')
  if(buf(1:4) == '%mem') then
   read(buf(k+1:j-2),*) mem
   select case(buf(j-1:j))
   case('gb')
    mem = 1000*mem ! You like 1024? I prefer 1000
   case('mb')
   case('gw')
    mem = 1000*8*mem
   case('mw')
    mem = 8*mem
   case default
    write(6,'(/,A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: memory&
                    & unit cannot be recognized.'
    write(6,'(A)') "Only 'GB', 'MB', 'GW', and 'MW' are accepted."
    write(6,'(A)') 'unit = '//TRIM(buf(j-1:j))
    stop
   end select
  else if(buf(1:6) == '%nproc') then
   read(buf(k+1:),*) np
  end if
 end do ! for while

 close(fid)

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: incomplete&
                & file '//TRIM(gjfname)
  stop
 end if
end subroutine read_mem_and_nproc_from_gjf

! read the method and basis set from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_method_and_basis_from_buf(buf, method, basis, wfn_type)
 implicit none
 integer :: i, j
 integer, intent(out) :: wfn_type ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 character(len=1200), intent(in) :: buf
 character(len=9), intent(out) :: method
 character(len=21), intent(out) :: basis

 j = index(buf, '/')
 if(j == 0) then
  write(6,'(/,A)') "ERROR in subroutine read_method_and_basis_from_buf: no '/'&
                   & symbol found in"
  write(6,'(A)') "the '#' line. Failed to identify the method name. This utilit&
                 &y does not support"
  write(6,'(A)') "support syntax like 'M062X cc-pVDZ'. You should use the '/' &
                 &symbol like 'M062X/cc-pVDZ'."
  stop
 end if

 i = index(buf(1:j-1), ' ', back=.true.)
 method = ' '
 method = buf(i+1:j-1)

 if(method(1:1) == 'u') then
  wfn_type = 3 ! UHF
  method = method(2:)
 else if(method(1:2) == 'ro') then
  wfn_type = 2 ! ROHF
  method = method(3:)
 else if(method(1:1) == 'r') then
  wfn_type = 1 ! RHF
  method = method(2:)
 else
  wfn_type = 0 ! undetermined
 end if

 i = index(buf(j+1:), ' ')
 basis = ' '
 basis = buf(j+1:j+i-1)
end subroutine read_method_and_basis_from_buf

! read nfrag (the number of fragments) from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_nfrag_from_buf(buf, nfrag)
 implicit none
 integer :: i, j, k
 integer, intent(out) :: nfrag
 character(len=1200), intent(in) :: buf

 i = index(buf, 'guess')
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nfrag_from_buf: keyword 'guess'&
                   & not found in buf."
  write(6,'(A)') 'buf = '//TRIM(buf)
  write(6,'(/,A)') "You should write 'guess(fragment=N)' in route section&
                     & to specify the number of fragments."
  stop
 end if

 i = i + 4
 nfrag = 0
 j = i+1 + index(buf(i+2:), '=')

 select case(buf(i+1:i+1))
 case('=')
  k = i+1 + index(buf(i+2:), ' ')
 case('(')
  k = i+1 + index(buf(i+2:), ')')
 case default
  write(6,'(A)') 'ERROR in subroutine read_nfrag_from_buf: syntax error.'
  write(6,'(A)') 'buf = '//TRIM(buf)
  stop
 end select
 read(buf(j+1:k-1),*) nfrag

 if(nfrag == 1) then
  write(6,'(A)') 'ERROR in subroutine read_nfrag_from_buf: nfrag = 1.'
  write(6,'(A)') 'Only one fragment. No need to use fragment-guess.'
  stop
 else if(nfrag < 1) then
  write(6,'(A,I0)') 'ERROR in subroutine read_nfrag_from_buf: invalid&
                    & nfrag = ', nfrag
  stop
 end if

end subroutine read_nfrag_from_buf

! generate a .gjf file from a type frag
subroutine gen_gjf_from_type_frag(frag0, guess_read, stab_chk, basname)
 use frag_info, only: frag
 use theory_level, only: mem, nproc, method, basis, scrf, sph, eda_type, disp_type
 implicit none
 integer :: i, k(3),fid
 character(len=9) :: method0, method1
 character(len=240), intent(in) :: basname
 type(frag), intent(in) :: frag0
 logical, intent(in) :: guess_read, stab_chk

 k = [index(method,'o3lyp'),index(method,'ohse2pbe'),index(method,'ohse1pbe')]

 if(ANY(k > 0)) then
  write(6,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: O3LYP, OHSE2PBE&
                   &, and OHSE1PBE'
  write(6,'(A)') 'functionals are not supported since these functional names ar&
                 &e confusing.'
  stop
 end if

 if(method(1:7)=='ccsd(t)' .and. frag0%wfn_type==2) then
  write(6,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: ROCCSD(T) not &
                    supported.'
  write(6,'(A)') 'You can use ROCCSD instead.'
  stop
 end if

 ! If MP2, CCSD or CCSD(T) is required in LMO-EDA, only R(O)HF orbitals are wr-
 ! itten into the generated .inp file. So we call Gaussian to compute R(O)HF
 if(method(1:3)=='mp2' .or. method(1:4)=='ccsd') then
  if(eda_type /= 2) then
   write(6,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: MP2 or CC meth&
                    &ods are only available'
   write(6,'(A)') 'for LMO-EDA, but you are not using LMO-EDA.'
   stop
  end if
  if(frag0%wfn_type > 2) then
   write(6,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: only R(O)HF ba&
                    &sed MP2 or CC methods'
   write(6,'(A)') 'are supported in LMO-EDA. It seems that you are trying to us&
                  &e UHF-based methds.'
   stop
  end if
  method1 = 'hf'
 else
  method1 = method
 end if

 method0 = ' '
 select case(frag0%wfn_type)
 case(1)
  method0 = 'R'//TRIM(method1)
 case(2)
  method0 = 'RO'//TRIM(method1)
 case(3)
  method0 = 'U'//TRIM(method1)
 case default
  write(6,'(A,I0)') 'ERROR in subroutine gen_gjf_from_type_frag: invalid&
                    & wfn_type=', frag0%wfn_type
  stop
 end select

 if(eda_type==1 .and. TRIM(method0)/='Rhf') then
  write(6,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: Morokuma-EDA&
                   & can only be applied'
  write(6,'(A)') 'to RHF. But found method='//TRIM(method0)
  stop
 end if

 open(newunit=fid,file=TRIM(frag0%fname),status='replace')
 i = index(frag0%fname, '.gjf', back=.true.)
 write(fid,'(A)') '%chk='//frag0%fname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'MB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p '//TRIM(method0)//'/'//TRIM(basis)//' nosymm'

 call lower(method0)
 select case(TRIM(method0))
 case('rhf','rohf','uhf')
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 case default
  write(fid,'(A)',advance='no') ' int(nobasistransform,ultrafine)'
 end select

 select case(frag0%wfn_type)
 case(1) ! RHF
  write(fid,'(A)',advance='no') ' scf(xqc,maxcycle=300)'
 case(2) ! ROHF
  write(fid,'(A)',advance='no') ' scf(maxcycle=300,NoIncFock,NoVarAcc)'
 case(3) ! UHF
  write(fid,'(A)',advance='no') ' scf(xqc,maxcycle=300,NoIncFock,NoVarAcc,conver=7)'
 case default
  write(6,'(A)') 'ERROR in subroutine gen_gjf_from_type_frag: wfn_type&
                 & out of range.'
  write(6,'(A,I0)') 'frag0%wfn_type=', frag0%wfn_type
  stop
 end select

 if(LEN_TRIM(scrf) > 0) write(fid,'(A)',advance='no') ' '//TRIM(scrf)

 if(sph) then
  write(fid,'(A)',advance='no') ' 5D 7F'
 else
  write(fid,'(A)',advance='no') ' 6D 10F'
 end if

 if(disp_type == 1) then
  write(fid,'(A)',advance='no') ' em=GD3'
 else if(disp_type == 2) then
  write(fid,'(A)',advance='no') ' em=GD3BJ'
 end if

 if(guess_read) then
  write(fid,'(A)') ' guess=read'
 else
  if(frag0%noiter) then
   write(fid,'(A)') ' guess(only,save)'
  else
   if(frag0%wfn_type==3 .and. frag0%mult==1) then
    write(fid,'(A)') ' guess=mix'
   else
    write(fid,'(/)',advance='no')
   end if
  end if
 end if

 write(fid,'(/,A,/)') 'generated by utility frag_guess_wfn of MOKIT'
 write(fid,'(I0,1X,I0)') frag0%charge, frag0%mult

 do i = 1, frag0%natom, 1
  write(fid,'(A)',advance='no') TRIM(frag0%elem(i))
  if(frag0%ghost(i)) write(fid,'(A3)',advance='no') '-Bq'
  write(fid,'(3(1X,F18.10))') frag0%coor(1:3,i)
 end do ! for i

 write(fid,'(/)',advance='no')
 close(fid)

 call copy_gen_basis_bas2gjf(basname, frag0%fname)

 if(ANY(frag0%ghost .eqv. .true.) .and. index(basis,'genecp')>0) &
  call del_ecp_of_ghost_in_gjf(frag0%fname)

 if((.not.guess_read) .and. (.not.frag0%noiter) .and. frag0%wfn_type==3) then
  open(newunit=fid,file=TRIM(frag0%fname),status='old',position='append')
  write(fid,'(/,A)') '--Link1--'
  i = index(frag0%fname, '.gjf', back=.true.)
  write(fid,'(A)') '%chk='//frag0%fname(1:i-1)//'.chk'
  write(fid,'(A,I0,A)') '%mem=', mem, 'MB'
  write(fid,'(A,I0)') '%nprocshared=', nproc
  write(fid,'(A)',advance='no') '#p '//TRIM(method0)//' chkbasis guess=read&
   & geom=allcheck nosymm scf(xqc,maxcycle=300,NoIncFock,NoVarAcc)'
  if(stab_chk) write(fid,'(A)',advance='no') ' stable=opt'

  if(LEN_TRIM(scrf) > 0) write(fid,'(A)',advance='no') ' '//TRIM(scrf)
  if(disp_type == 1) then
   write(fid,'(A)',advance='no') ' em=GD3'
  else if(disp_type == 2) then
   write(fid,'(A)',advance='no') ' em=GD3BJ'
  end if

  select case(TRIM(method0))
  case('rhf','rohf','uhf')
   write(fid,'(A,/)') ' int=nobasistransform'
  case default
   write(fid,'(A,/)') ' int(nobasistransform,ultrafine)'
  end select
  close(fid)
 end if
end subroutine gen_gjf_from_type_frag

! generate extended basis set gjf files for a fragment
subroutine gen_extend_bas_frag(frag_i, frag_n, frag_k)
 use frag_info, only: frag
 implicit none
 integer :: i, natom
 type(frag), intent(in) :: frag_i, frag_n
 type(frag), intent(out) :: frag_k

 frag_k%charge   = frag_i%charge
 frag_k%mult     = frag_i%mult
 frag_k%wfn_type = frag_i%wfn_type
 frag_k%pos      = frag_i%pos
 frag_k%coor     = frag_n%coor
 frag_k%elem     = frag_n%elem

 natom = frag_n%natom
 frag_k%natom = natom
 allocate(frag_k%ghost(natom))
 frag_k%ghost = .true.

 do i = 1, natom, 1
  if( ANY(frag_i%atm_map == i) ) frag_k%ghost(i) = .false.
 end do ! for i
end subroutine gen_extend_bas_frag

! generate .fch files from .chk files, and create GAMESS .inp file (if EDA requested)
subroutine gen_inp_of_frags()
 use frag_info, only: nfrag0, nfrag, frags
 use util_wrapper, only: fch2psi_wrap, fch2inp_wrap
 use theory_level, only: nproc, wfn_type, eda_type, scrf, method
 use periodic_table, only: read_elem_from_gjf, elem2vdw_radii
 implicit none
 integer :: i, k
 integer :: natom  ! No. of atoms
 real(kind=8), allocatable :: radii(:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: buf, fchname, fileA, fileB
 character(len=240) :: inpname1, inpname2, outname
 logical :: r2u
 logical, allocatable :: ghost(:)

 if(eda_type == 0) return
 natom = frags(nfrag)%natom

 k = index(frags(1)%fname, '.gjf', back=.true.)
 ! do not change k below in this subroutine
 inpname1 = ' '; inpname2 = ' '; outname = ' '

 do i = 1, nfrag, 1
  fchname = frags(i)%fname(1:k-1)//'.fch'
  if(eda_type == 4) then ! SAPT in PSI4
   call fch2psi_wrap(fchname)
   inpname1 = frags(i)%fname(1:k-1)//'_psi.inp'
   if(i < nfrag) call delete_files(2, [fchname, inpname1])
  else
   call fch2inp_wrap(fchname, .false., 0, 0)
   if(eda_type==1 .and. i==nfrag) call delete_file(fchname)
  end if
 end do ! for i

 i = index(frags(nfrag)%fname, '-')
 if(eda_type == 4) then ! SAPT in PSI4
  inpname1 = frags(nfrag)%fname(1:k-1)//'_psi.inp'
 else
  inpname1 = frags(nfrag)%fname(1:k-1)//'.inp'
 end if
 inpname2 = frags(nfrag)%fname(i+1:k-1)//'.inp'

 allocate(radii(natom), source=0d0)

 if(LEN_TRIM(scrf) > 0) then
  if(index(scrf,'smd') == 0) then
   buf = frags(nfrag)%fname(i+1:k+3)
   call read_natom_from_gjf(buf, natom)
   allocate(elem(natom), ghost(natom))
   call read_elem_from_gjf(buf, natom, elem, ghost)
   deallocate(ghost)
   forall(i = 1:natom) radii(i) = elem2vdw_radii(elem(i))
   deallocate(elem)
  end if
 end if

 if(eda_type == 4) then ! SAPT in PSI4
  do i = 1, nfrag0, 1
   if(wfn_type==3 .and. frags(i)%wfn_type/=3) then ! RHF copied to UHF
    fileA = frags(i)%fname(1:k-1)//'.A'
    fileB = frags(i)%fname(1:k-1)//'.B'
    call copy_file(fileA, fileB, .false.)
   end if
  end do ! for i

  inpname1 = frags(nfrag)%fname(1:k-1)//'_psi.inp'
  call copy_and_modify_psi4_sapt_file(inpname1, inpname2)
 else ! not SAPT

  call copy_and_modify_gms_eda_file(natom, radii, inpname1, inpname2)
  select case(eda_type)
  case(1)
   call delete_file(inpname1)
  case(2,3) ! LMO/GKS-EDA
   call copy_vec_to_append_another_inp(inpname1, inpname2, 0, .false., .true.,&
                                       .true., .false.)
  end select

  do i = 1, nfrag0, 1
   r2u = .false.
   if(wfn_type==3 .and. frags(i)%wfn_type/=3) r2u = .true.
   inpname1 = frags(i)%fname(1:k-1)//'.inp'
   call copy_vec_to_append_another_inp(inpname1, inpname2, i, .false., .true.,&
                                       .true., r2u)

   if(eda_type==2 .or. eda_type==3) then ! LMO/GKS-EDA
    inpname1 = frags(i+nfrag0)%fname(1:k-1)//'.inp'
    call copy_vec_to_append_another_inp(inpname1, inpname2, i, .true., .true.,&
                                        .true., r2u)
   end if
  end do ! for i
 end if

 deallocate(radii, frags)
 write(6,'(/,A)',advance='no') 'Done. Now you can submit file '//TRIM(inpname2)&
                              //' to '

 i = index(inpname2, '.inp',back=.true.)

 if(eda_type == 4) then ! SAPT in PSI4
  outname = inpname2(1:i-1)//'.out'
  write(6,'(A,/,A)') 'PSI4 (DO NOT delete *.A and', '*.B files) like'
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A,I0,A)') 'psi4 '//TRIM(inpname2)//' '//TRIM(outname)//' -n ',&
                       nproc," &"
 else
  outname = inpname2(1:i-1)//'.gms'
  write(6,'(A)') 'GAMESS or XEDA like'
  if(wfn_type==2 .and. method(1:4)=='ccsd') then ! ROCCSD
   i = 1
  else
   i = nproc
  end if
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A,I0,A)') 'xeda '//TRIM(inpname2)//' 00 ',i,' >'//TRIM(outname)&
                      //" 2>&1 &"
 end if

 write(6,'(A)') REPEAT('-',79)
end subroutine gen_inp_of_frags

! copy content of a provided .inp file and modify it to be SAPT job
subroutine copy_and_modify_psi4_sapt_file(inpname1, inpname2)
 use frag_info, only: nfrag0, nfrag, frags
 use theory_level, only: mem, basis, wfn_type
 use periodic_table, only: elem2nuc
 implicit none
 integer :: i, j, natom, fid1, fid2
 character(len=2) :: str
 character(len=240) :: buf, fileA, fileB
 character(len=240), intent(in) :: inpname1, inpname2

 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')
 write(fid2,'(A)') '# generated by utility frag_guess_wfn of MOKIT'
 write(fid2,'(A,I0,A)') 'memory ',mem,' MB'
 write(fid2,'(/,A)') 'molecule mymol {'
 write(fid2,'(A)') 'symmetry c1'
 write(fid2,'(A)') 'no_reorient'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:11) == 'no_reorient') exit
 end do ! for while
 read(fid1,'(A)') buf

 do i = 1, nfrag0, 1
  if(i > 1) write(fid2,'(A)') '--'
  write(fid2,'(I0,1X,I0)') frags(i)%charge, frags(i)%mult
  natom = COUNT(frags(i)%ghost .eqv. .false.)
  do j = 1, natom, 1
   read(fid1,'(A)') buf
   write(fid2,'(A)') TRIM(buf)
  end do ! for j
 end do ! for i
 read(fid1,'(A)') buf
 write(fid2,'(A)') '}'

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:1) == '}') exit
 end do ! for while
 close(fid1, status='delete')

 write(fid2,'(/,A)') 'dimer = psi4.get_active_molecule()'
 write(fid2,'(/,A)') 'set {'
 write(fid2,'(A)') ' scf_type df' ! pk is slow and disk-consuming
 write(fid2,'(A)') ' s_tolerance 1e-6'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'
 i = frags(3)%wfn_type
 select case(i)
 case(1)
  write(fid2,'(A)') ' reference rhf'
 case(2)
  write(fid2,'(A)') ' reference rohf'
 case(3)
  write(fid2,'(A)') ' reference uhf'
 case default
  write(6,'(A)') 'ERROR in subroutine copy_and_modify_psi4_sapt_file: wfn_ty&
                 &pe out of range.'
  write(6,'(A,I0)') 'wfn_type=', i
  stop
 end select
 write(fid2,'(A)') '}'
 write(fid2,'(A)') 'df_basis_sapt {'
 if(basis(1:4) == 'def2') basis = 'def2-'//TRIM(basis(5:))
 write(fid2,'(A)') ' assign '//TRIM(basis)//'-ri'
 if(TRIM(basis) == 'jun-cc-pvdz') then
  natom = frags(3)%natom
  do i = 1, natom, 1
   str = frags(3)%elem(i); j = elem2nuc(str)
   if(i > 2) then
    if(ANY(frags(3)%elem(1:i-1) == str)) cycle
   end if
   if(j==3 .or. j==11 .or. (j>18 .and. j<31) .or. j>36) then
    write(fid2,'(A)') ' assign '//TRIM(str)//' def2-qzvpp-ri'
   end if
  end do ! for i
 end if
 write(fid2,'(A)') '}'

 write(fid2,'(/,A)') 'set df_ints_io save'
 write(fid2,'(A)') "psi4.IO.set_default_namespace('dimer')"
 write(fid2,'(A)') "Edim, wfn_dimer = energy('scf', molecule=dimer, return_wfn=True)"
 write(fid2,'(A)') 'set df_ints_io load'
 write(fid2,'(A)') '# the above scf makes every array allocated'
 i = index(inpname2, '.inp', back=.true.)
 write(fileA,'(I3.3,A)') nfrag, '-'//inpname2(1:i-1)//'.A'
 write(fileB,'(I3.3,A)') nfrag, '-'//inpname2(1:i-1)//'.B'
 write(fid2,'(A)') "wfn_dimer.Ca().load('"//TRIM(fileA)//"')"
 if(wfn_type == 3) write(fid2,'(A)') "wfn_dimer.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'wfn_dimer.to_file(wfn_dimer.get_scratch_filename(180))'
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "Edim, wfn_dimer = energy('scf', molecule=dimer, return_wfn=True)"

 write(fid2,'(/,A)') 'monomerA = dimer.extract_subsets(1,2)'
 write(fid2,'(A)') "psi4.IO.change_file_namespace(97, 'dimer', 'monomerA')"
 write(fid2,'(A)') "psi4.IO.set_default_namespace('monomerA')"
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonA, wfn_monA = energy('scf', molecule=monomerA, return_wfn=True)"
 write(fid2,'(A)') '# the above scf makes every array allocated'
 i = index(inpname2, '.inp', back=.true.)
 write(fileA,'(A)') '001-'//inpname2(1:i-1)//'.A'
 write(fileB,'(A)') '001-'//inpname2(1:i-1)//'.B'
 write(fid2,'(A)') "wfn_monA.Ca().load('"//TRIM(fileA)//"')"
 if(wfn_type == 3) write(fid2,'(A)') "wfn_monA.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'wfn_monA.to_file(wfn_monA.get_scratch_filename(180))'
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonA, wfn_monA = energy('scf', molecule=monomerA, return_wfn=True)"

 write(fid2,'(/,A)') 'monomerB = dimer.extract_subsets(2,1)'
 write(fid2,'(A)') "psi4.IO.change_file_namespace(97, 'monomerA', 'monomerB')"
 write(fid2,'(A)') "psi4.IO.set_default_namespace('monomerB')"
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonB, wfn_monB = energy('scf', molecule=monomerB, return_wfn=True)"
 write(fid2,'(A)') '# the above scf makes every array allocated'
 i = index(inpname2, '.inp', back=.true.)
 write(fileA,'(A)') '002-'//inpname2(1:i-1)//'.A'
 write(fileB,'(A)') '002-'//inpname2(1:i-1)//'.B'
 write(fid2,'(A)') "wfn_monB.Ca().load('"//TRIM(fileA)//"')"
 if(wfn_type == 3) write(fid2,'(A)') "wfn_monB.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'wfn_monB.to_file(wfn_monB.get_scratch_filename(180))'
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonB, wfn_monB = energy('scf', molecule=monomerB, return_wfn=True)"

 write(fid2,'(/,A)') "psi4.IO.change_file_namespace(97, 'monomerB', 'dimer')"
 write(fid2,'(A)') "psi4.IO.set_default_namespace('dimer')"

 write(fid2,'(/,A)') "aux_basis = psi4.core.BasisSet.build(wfn_dimer.molecule(),&
                   & ""DF_BASIS_SAPT"","
 write(fid2,'(A)') " psi4.core.get_global_option(""DF_BASIS_SAPT""),"
 write(fid2,'(A)') " ""RIFIT"", psi4.core.get_global_option(""BASIS""))"
 write(fid2,'(A)') "wfn_dimer.set_basisset(""DF_BASIS_SAPT"", aux_basis)"
 write(fid2,'(A)') "wfn_dimer.set_basisset(""DF_BASIS_ELST"", aux_basis)"

 write(fid2,'(/,A)') 'psi4.sapt(wfn_dimer, wfn_monA, wfn_monB)'
 close(fid2)
end subroutine copy_and_modify_psi4_sapt_file

! copy content of a provided .inp file and modify it to be EDA job
subroutine copy_and_modify_gms_eda_file(natom, radii, inpname1, inpname2)
 use frag_info, only: nfrag0, frags
 use theory_level, only: mem, nproc, method, eda_type, scrf, solvent, &
  solvent_gau, disp_type
 implicit none
 integer :: i, j, m, fid1, fid2
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: radii(natom)
 character(len=9) :: dft_in_gms
 character(len=240) :: buf1, buf2
 character(len=240), intent(in) :: inpname1, inpname2

 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')

 read(fid1,'(A)') buf1
 i = index(buf1, 'RUNTYP=ENERGY')
 buf2 = buf1(1:i+6)//'EDA'//TRIM(buf1(i+13:))
 write(fid2,'(A)') TRIM(buf2)

 call convert_dft_name_gau2gms(method, dft_in_gms)
 read(fid1,'(A)') buf1

 i = index(buf1, '$END')
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine gen_inp_of_frags: no '$END' found&
                 & in the 2nd"
  write(6,'(A)') 'line of file '//TRIM(inpname1)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)',advance='no') buf1(1:i-1)//'DFTTYP='//TRIM(dft_in_gms)
 if(method(1:3) == 'mp2') then
  write(fid2,'(A)',advance='no') ' MPLEVL=2'
 else if(TRIM(method) == 'ccsd') then
  write(fid2,'(A)',advance='no') ' CCTYP=CCSD'
 else if(TRIM(method) == 'ccsd(t)') then
  write(fid2,'(A)',advance='no') ' CCTYP=CCSD(T)'
 end if
 write(fid2,'(A)') ' $END'

 read(fid1,'(A)') buf1 ! this line should be $SYSTEM
 i = FLOOR(DBLE(mem)*0.95d0/(8d0*DBLE(nproc)))
 write(fid2,'(A,I0)',advance='no') ' $SYSTEM MWORDS=',i

 if(method(1:3)=='mp2' .or. method(1:4)=='ccsd') then
  i = FLOOR(DBLE(mem)*0.05d0/8d0)
  write(fid2,'(A,I0)',advance='no') ' MEMDDI=',i
 end if
 write(fid2,'(A)') ' $END'

 if(TRIM(dft_in_gms) /= 'NONE') then
  write(fid2,'(A)',advance='no') ' $DFT NRAD0=99 NLEB0=590 NRAD=99 NLEB=590'
  if(disp_type == 1) then
   write(fid2,'(A)',advance='no') ' IDCVER=3' ! GD3
! Some dispersion parameters are not built-in in GAMESS, and some built-in
! dispersion parameters are wrong. So I have to explicitly add them
   select case(TRIM(dft_in_gms))
   case('B3PW91')
    write(fid2,'(A)',advance='no') ' DCSR=1.176 DCS8=1.775'
   case('BMK')
    write(fid2,'(A)',advance='no') ' DCSR=1.931 DCS8=2.168'
   case('BPBE')
    write(fid2,'(A)',advance='no') ' DCSR=1.087 DCS8=2.033'
   case('CAMB3LYP')
    write(fid2,'(A)',advance='no') ' DCSR=1.378 DCS8=1.217'
   case('M05')
    write(fid2,'(A)',advance='no') ' DCSR=1.373 DCS8=0.595'
   case('M05-2X')
    write(fid2,'(A)',advance='no') ' DCSR=1.417 DCS8=1d-9'
   case('M06')
    write(fid2,'(A)',advance='no') ' DCSR=1.325 DCS8=1d-9'
   case('M06L')
    write(fid2,'(A)',advance='no') ' DCSR=1.581 DCS8=1d-9'
   case('M06-HF')
    write(fid2,'(A)',advance='no') ' DCSR=1.446 DCS8=1d-9'
   case('M06-2X')
    write(fid2,'(A)',advance='no') ' DCSR=1.619 DCS8=1d-9'
   case('M11L')
    write(fid2,'(A)',advance='no') ' DCSR=2.3933 DCS8=1.1129'
   case('MN15L')
    write(fid2,'(A)',advance='no') ' DCSR=3.3388 DCS8=1d-9'
   case('PBE0')
    write(fid2,'(A)',advance='no') ' DCSR=1.287 DCS8=0.928'
   case('revTPSS')
    write(fid2,'(A)',advance='no') ' DCSR=1.3491 DCS8=1.3666'
   end select
  else if(disp_type == 2) then
   write(fid2,'(A)',advance='no') ' IDCVER=4' ! GD3BJ
  end if
  write(fid2,'(A)') ' $END'
 end if

 do while(.true.)
  read(fid1,'(A)') buf1
  if(buf1(2:7) == '$GUESS') exit
  if(buf1(2:5) == '$SCF') then
   if(eda_type == 1) then
    buf1 = ' $SCF DIRSCF=.F. $END'
   else
    if(TRIM(dft_in_gms) /= 'NONE') buf1 = ' $SCF DIRSCF=.T. DIIS=.F. SOSCF=.T. $END'
   end if
  end if
  write(fid2,'(A)') TRIM(buf1)
 end do ! for while

 select case(eda_type)
 case(1) ! MOROKUMA-EDA
  write(fid2,'(A)',advance='no') ' $MOROKM IATM(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%natom
  end do ! for i

  write(fid2,'(A)',advance='no') ' ICHM(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%charge
  end do ! for i

  write(fid2,'(A)') ' $END'

 case(2,3) ! 2/3 for LMO-EDA/GKS-EDA, respectively
  write(fid2,'(A)',advance='no') ' $LMOEDA MATOM(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%natom
  end do ! for i

  write(fid2,'(A)',advance='no') ' MCHARG(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%charge
  end do ! for i

  write(fid2,'(A)',advance='no') ' MMULT(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   if(frags(i)%pos) then
    write(fid2,'(I0)',advance='no') frags(i)%mult
   else
    write(fid2,'(I0)',advance='no') -frags(i)%mult
   end if
  end do ! for i

  if(eda_type == 2) then
   write(fid2,'(/,A)',advance='no') '  EDATYP=NONE'
  else
   write(fid2,'(/,A)',advance='no') '  EDATYP=GKS'
  end if
  write(fid2,'(A,/,A)') ' RDVECM=.T. $END',' $GUESS GUESS=HCORE $END'

 case default
  write(6,'(A,I0)') 'ERROR in subroutine gen_inp_of_frags: invalid eda_type=',&
                     eda_type
  close(fid1)
  close(fid2)
  stop
 end select

 if(LEN_TRIM(scrf) > 0) then
  write(fid2,'(A)',advance='no') ' $PCM IEF=-3 SOLVNT='//TRIM(solvent)

  if(index(scrf,'smd') > 0) then
   write(fid2,'(A)') ' SMD=.T. $END'
   write(6,'(/,A)') 'Warning from subroutine gen_inp_of_frags: scrf=SMD detected.'
   write(6,'(A)') 'If you encounter SCF convergence problems later in GAMESS,&
                  & you can change'
   write(6,'(A)') 'to scrf=PCM or remove scrf=SMD (i.e. gas phase computation).'
  else ! PCM
   if(TRIM(solvent) == 'INPUT') then
    select case(TRIM(solvent_gau))
    case('n,n-dimethylacetamide')
     write(fid2,'(A)',advance='no') ' RSOLV=1.935 EPS=37.781'
    case('n,n-dimethylformamide')
     write(fid2,'(A)',advance='no') ' RSOLV=1.885 EPS=37.219'
    end select
   end if
   write(fid2,'(A,/,A)') ' $END',' $PCMCAV ALPHA(1)=1.1'
   j = natom/8
   do i = 1, j, 1
    write(fid2,'(A,I0,A,F6.4,7(A1,F6.4))') '  RIN(',8*i-7,')=', radii(8*i-7),&
                                                  (',',radii(m),m=8*i-6,8*i)
   end do ! for i
   if(natom-8*j > 0) then
    write(fid2,'(A,I0,A,F6.4,7(A1,F6.4))') '  RIN(',8*j+1,')=', radii(8*j+1),&
                                                (',',radii(m),m=8*j+2,natom)
   end if
   write(fid2,'(A)') ' $END'
  end if
 end if

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf1
  if(i /= 0) exit
  if(buf1(2:5) == '$VEC') exit
  write(fid2,'(A)') TRIM(buf1)
 end do ! for while

 close(fid1)
 close(fid2)

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine copy_and_modify_gms_eda_file: no '$VEC&
                &' found in file "//TRIM(inpname1)
  stop
 end if
end subroutine copy_and_modify_gms_eda_file

! copy $VEC from file inpname1 into inpname2, by appending
subroutine copy_vec_to_append_another_inp(inpname1, inpname2, ivec, extended, &
                                          deleted, occ, r2u)
 implicit none
 integer :: i, nline, fid1, fid2
 integer :: na, nb, nif, nbf, nbf1
 ! Note: nbf1>=nbf always holds. nbf1 is always the number of basis functions
 !  of Cartesian functions, but nbf depends on the keywords(ISPHER=1) used
 integer, intent(in) :: ivec ! ivec>=-1
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname1, inpname2
 logical :: uhf
 logical, intent(in) :: extended ! $VEC or $VEA
 logical, intent(in) :: deleted  ! whether to delete inpname1
 logical, intent(in) :: occ      ! whether to copy only occupied orbitals or all orbitals
 logical, intent(in) :: r2u      ! whether to copy nb orbitals from Alpha MOs

 call check_uhf_in_gms_inp(inpname1, uhf)
 if(uhf .and. r2u) then
  write(6,'(A)') "ERROR in subroutine copy_vec_to_append_another_inp: both&
                 & logical variables 'uhf' and"
  write(6,'(A)') "'r2u' are .True. This is not allowed."
  write(6,'(A)') 'inpname1='//TRIM(inpname1)//', inpname2='//TRIM(inpname2)
  stop
 end if

 if(occ) then
  call read_na_nb_nif_nbf_from_gms_inp(inpname1, na, nb, nif, nbf)
  call read_nbf_from_dat(inpname1, nbf1)
!  write(*,'(5(A,I0))') 'na=',na,',nb=',nb,',nif=',nif,',nbf=',nbf,',nbf1=',nbf1
 end if

 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == '$VEC') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine copy_vec_to_append_another_inp:'
  write(6,'(A)') "No '$VEC' found in file "//TRIM(inpname1)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(inpname2),status='old',position='append')
 if(extended) then
  if(ivec > 0) then
   write(fid2,'(A,I0)') ' $VEA',ivec
  else
   write(6,'(A)') 'ERROR in subroutine copy_vec_to_append_another_inp:'
   write(6,'(A)') 'ivec<=0 when extended=.True. This is not allowed.'
   close(fid1)
   close(fid2)
   stop
  end if
 else ! not extended basis set
  if(ivec > 0) then
   write(fid2,'(A,I0)') ' $VEC',ivec
  else if(ivec == 0) then
   write(fid2,'(A)') ' $VEC0'
  else
   write(fid2,'(A)') ' $VEC'
  end if
 end if

 if(occ) then ! only copy occupied MOs of a fragment
  nline = nbf1/5
  if(nbf1-5*nline > 0) nline = nline + 1

  do i = 1, nline*na, 1
   read(fid1,'(A)') buf
   write(fid2,'(A)') TRIM(buf)
  end do ! for i

  if(uhf) then
   do i = 1, nline*(nif-na), 1
    read(fid1,'(A)') buf
   end do ! for i
   do i = 1, nline*nb, 1
    read(fid1,'(A)') buf
    write(fid2,'(A)') TRIM(buf)
   end do ! for i
  end if

  if(r2u) then ! copy beta MOs from alpha MOs
   do while(.true.)
    BACKSPACE(fid1)
    BACKSPACE(fid1)
    read(fid1,'(A)') buf
    if(buf(2:5) == '$VEC') exit
   end do ! for while
   do i = 1, nline*nb, 1
    read(fid1,'(A)') buf
    write(fid2,'(A)') TRIM(buf)
   end do ! for i
  end if

 else ! copy all MOs of a fragment

  do while(.true.)
   read(fid1,'(A)') buf
   if(buf(2:5) == '$END') exit
   write(fid2,'(A)') TRIM(buf)
  end do ! for while

  if(r2u) then ! copy beta MOs from alpha MOs
   do while(.true.)
    BACKSPACE(fid1)
    BACKSPACE(fid1)
    read(fid1,'(A)') buf
    if(buf(2:5) == '$VEC') exit
   end do ! for while
   do while(.true.)
    read(fid1,'(A)') buf
    if(buf(2:5) == '$END') exit
    write(fid2,'(A)') TRIM(buf)
   end do ! for while
  end if
 end if

 write(fid2,'(A)') ' $END'
 if(deleted) then
  close(fid1,status='delete')
 else
  close(fid1)
 end if

 close(fid2)
end subroutine copy_vec_to_append_another_inp

! read the scrf string from a longbuf
subroutine read_scrf_string_from_buf(longbuf, scrf)
 implicit none
 integer :: i, j
 character(len=1200), intent(in) :: longbuf
 character(len=60), intent(out) :: scrf

 scrf = ' '
 i = index(longbuf, 'scrf')
 if(i == 0) return

 j = index(longbuf(i:), ' ')
 if(j == 0) then
  j = i + 3
 else
  j = j + i - 2
 end if

 scrf = longbuf(i:j)
end subroutine read_scrf_string_from_buf

! determine implicit solvent model from the scrf string
! Note: the input scrf must be in lower case
subroutine determine_solvent_from_gau2gms(scrf, solvent)
 use theory_level, only: solvent_gau
 implicit none
 integer :: i, j, k(3)
 character(len=60), intent(in) :: scrf
 character(len=40), intent(out) :: solvent

 solvent = ' '
 if(LEN_TRIM(scrf) == 0) return

 i = index(scrf,'read')
 if(i > 0) then
  write(6,'(/,A)') "ERROR in subroutine determine_solvent_from_gau2gms: 'read'&
                  & in 'scrf()' is not"
  write(6,'(A)') 'supported for frag_guess_wfn.'
  stop
 end if

 i = index(scrf,'solvent')
 if(i == 0) then
  solvent = 'H2O'
  return
 end if

 k = [index(scrf(i+10:),')'),index(scrf(i+10:),','),index(scrf(i+10:),' ')]
 forall(j=1:3, k(j)==0) k(j) = LEN_TRIM(scrf)
 j = MINVAL(k)
 solvent_gau = ' '
 solvent_gau = scrf(i+8:i+8+j)

 select case(TRIM(solvent_gau))
 case('ch3ch2oh')
  solvent = 'C2H5OH'
 case('chloroform')
  solvent = 'CHCl3'
 case('carbontetrachloride')
  solvent = 'CCl4'
 case('methylenechloride')
  solvent = 'CH2Cl2'
 case('dichloroethane','1,2-dichloroethane','ch2clch2cl')
  solvent = 'C2H4Cl2'
 case('heptane','n-heptane')
  solvent = 'C7H16'
 case('acetone','ch3coch3')
  solvent = 'CH3COCH3'
 case('tetrahydrofuran')
  solvent = 'THF'
 case('dimethylsulfoxide','ch3soch3')
  solvent = 'DMSO'
 case('n-hexane')
  solvent = 'HEXANE'
 case('1,4-dioxane')
  solvent = 'dioxane'
 case('n,n-dimethylacetamide','n,n-dimethylformamide')
  solvent = 'INPUT'
 case('water','h2o','methanol','ch3oh','ethanol','hexane','acetonitrile','dmso',&
      'thf','benzene','c6h6','nitromethane','ch3no2','aniline','c6h5nh2',&
      'cyclohexane','c6h12','ccl4','dichloromethane','ch2cl2','toluene','c6h5ch3',&
      'chlorobenzene','c6h5cl')
  solvent = TRIM(solvent_gau)
 case default
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning in subroutine determine_solvent_from_gau2gms: implic&
                 &it solvent name in'
  write(6,'(A)') 'Gaussian .gjf file cannot be identified. You should open GAM&
                 &ESS .inp file and'
  write(6,'(A)') 'add it by yourself. Solvent in .gjf file: '//TRIM(solvent_gau)
  write(6,'(A)') REPEAT('-',79)
  solvent = 'INPUT'
 end select
end subroutine determine_solvent_from_gau2gms

! convert DFT name of Gaussian to that of GAMESS
subroutine convert_dft_name_gau2gms(method, dft_in_gms)
 implicit none
 character(len=9), intent(in) :: method
 character(len=9), intent(out) :: dft_in_gms

 dft_in_gms = ' '

 select case(TRIM(method))
 case('b3lyp','blyp','x3lyp','b3p86','b3pw91','wB97','wB97X','b98','m05','m06',&
      'm11','bmk','bp86','tpssh','revtpss')
  dft_in_gms = method
 case('m052x')
  dft_in_gms = 'M05-2X'
 case('m06hf')
  dft_in_gms = 'M06-HF'
 case('m062x')
  dft_in_gms = 'M06-2X'
 case('m06l')
  dft_in_gms = 'M06-L'
 case('m08hx')
  dft_in_gms = 'M08-HX'
 case('m11l')
  dft_in_gms = 'M11-L'
 case('wb97xd')
  dft_in_gms = 'wB97X-D'
 case('b971')
  dft_in_gms = 'B97-1'
 case('b972')
  dft_in_gms = 'B97-2'
 case('b97d')
  dft_in_gms = 'B97-D'
 case('pbepbe')
  dft_in_gms = 'PBE'
 case('pbe1pbe')
  dft_in_gms = 'PBE0'
 case('cam-b3lyp')
  dft_in_gms = 'CAMB3LYP'
 case('tpsstpss')
  dft_in_gms = 'TPSS'
 case('hf','rhf','rohf','uhf','mp2','ccsd','ccsd(t)')
  dft_in_gms = 'NONE'
 case default
  write(6,'(A)') 'Warning in subroutine convert_dft_name_gau2gms: functional na&
                 &me cannot be'
  write(6,'(A)') 'recognized. DFTTYP=NONE will be set in .inp file. You can mod&
                 &ify it by yourself.'
  dft_in_gms = 'NONE'
 end select

end subroutine convert_dft_name_gau2gms

! delete ECP/PP of ghost atoms in a given .gjf file
subroutine del_ecp_of_ghost_in_gjf(gjfname)
 use periodic_table, only: read_elem_from_gjf
 implicit none
 integer :: i, j, nbat1, nbat2, fid, fid1, RENAME
 integer :: nblank, natom
 character(len=2), allocatable :: elem(:)
 character(len=3) :: str
 character(len=6) :: str1
 character(len=240) :: buf, gjfname1
 character(len=240), intent(in) :: gjfname
 logical :: skipped
 logical, allocatable :: ghost(:) ! .True./.False. for ghost atom or not

 call read_natom_from_gjf(gjfname, natom)
 allocate(elem(natom), ghost(natom))
 call read_elem_from_gjf(gjfname, natom, elem, ghost)

 gjfname1 = TRIM(gjfname)//'.t'
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(gjfname1),status='replace')
 nblank = 0

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 4) exit
 end do ! for while

 read(fid,'(A)') buf

 if(LEN_TRIM(buf) > 0) then
  do while(.true.)
   read(buf,*) str
   str = ADJUSTL(str)
   if(str(1:1) == '-') str = str(2:)

   j = IACHAR(str(1:1))
   if(j>96 .and. j<123) str(1:1) = ACHAR(j-32)
   j = IACHAR(str(2:2))
   if(j>64 .and. j<91) str(2:2) = ACHAR(j+32)

   skipped = .false.
   do i = 1, natom, 1
    if(elem(i) == TRIM(str)) then
     if(ghost(i)) then
      skipped = .true.
      exit
     end if
    end if
   end do ! for i

   if(.not. skipped) write(fid1,'(A)') TRIM(buf)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str1, nbat1

   if(skipped) then
    if(i == 0) then ! ECP/PP data, not name
     nbat1 = nbat1 + 1
     do i = 1, nbat1, 1
      read(fid,'(A)') buf
      read(fid,*) nbat2
      do j = 1, nbat2, 1
       read(fid,'(A)') buf
      end do ! for j
     end do ! for i
    end if
   else ! not skipped
    write(fid1,'(A)') TRIM(buf)

    if(i == 0) then ! ECP/PP data, not name
     nbat1 = nbat1 + 1
     do i = 1, nbat1, 1
      read(fid,'(A)') buf
      write(fid1,'(A)') TRIM(buf)
      read(fid,'(A)') buf
      write(fid1,'(A)') TRIM(buf)
      read(buf,*) nbat2
      do j = 1, nbat2, 1
       read(fid,'(A)') buf
       write(fid1,'(A)') TRIM(buf)
      end do ! for j
     end do ! for i
    end if
   end if

   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do ! for while
 end if

 write(fid1,'(/)',advance='no')
 close(fid,status='delete')
 close(fid1)
 deallocate(elem, ghost)
 i = RENAME(TRIM(gjfname1), TRIM(gjfname))
end subroutine del_ecp_of_ghost_in_gjf

! construct supermolecule occ MOs from direct sum of fragment occ MOs
! construct supermolecule vir MOs by PAO construction
subroutine direct_sum_frag_mo2super_mo(n, gjfname0, wfn_type0, pos, gjfname, &
                                       wfn_type, stab_chk)
 use util_wrapper, only: formchk, unfchk
 implicit none
 integer :: i, j, fid0, fid, k1, k2, k5, RENAME
 integer :: nbf0, nif0, nbf, nif
 integer :: na0, nb0, na, nb   ! No. of alpha/beta electrons
 integer, intent(in) :: n
 integer, intent(in) :: wfn_type0(n), wfn_type ! 1/2/3 for RHF/ROHF/UHF
 character(len=240) :: buf, fchname0, fchname, chkname
 character(len=240), intent(in) :: gjfname0(n), gjfname
 real(kind=8), allocatable :: mo_a0(:,:), mo_b0(:,:) ! fragment MOs
 real(kind=8), allocatable :: mo_a(:,:), mo_b(:,:)   ! supermolecule MOs
 real(kind=8), allocatable :: ovlp(:,:) ! AO overlap matrix
 logical :: alive
 logical, intent(in) :: pos(n), stab_chk

 i = index(gjfname, '.gjf', back=.true.)
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine direct_sum_frag_mo2super_mo: no '.gjf'&
                & string found in file "//TRIM(gjfname)
  stop
 end if

 if(ANY(wfn_type0 == 0) .or. wfn_type==0) then
  write(6,'(A)') 'ERROR in subroutine direct_sum_frag_mo2super_mo: there exists&
                & at least one file'
  write(6,'(A)') 'whose wfn_type is 0. Not allowed.'

  do i = 1, n, 1
   write(6,'(2(A,I0))') 'i=',i,','//TRIM(gjfname0(i))//',wfn_type0(i)=',wfn_type0(i)
  end do ! for i
  write(6,'(A,I0)') 'wfn_type=', wfn_type
  stop
 end if

 ! modify 'guess(only,save)' in gjf file into 'guess=read geom=allcheck'
 fchname = TRIM(gjfname)//'.t'
 open(newunit=fid0,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid,file=TRIM(fchname),status='replace')
 do while(.true.)
  read(fid0,'(A)') buf
  i = index(buf,'guess')
  if(i > 0) then
   buf = buf(1:i+4)//'=read geom=allcheck'
   i = index(buf,'/')
   j = index(buf(i+1:),' ')
   buf = buf(1:i)//'chkbasis'//buf(i+j:)

   i = index(buf,',conver=')
   if(i > 0) then
    j = index(buf(i+1:),')')
    buf = buf(1:i-1)//')'//buf(i+j+1:)
   end if

   if(wfn_type==3 .and. stab_chk) buf = TRIM(buf)//' stable=opt'
  end if
  write(fid,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while
 close(fid0,status='delete')
 close(fid)
 i = RENAME(TRIM(fchname), TRIM(gjfname))
 ! modify done

 i = index(gjfname, '.gjf', back=.true.)
 fchname = gjfname(1:i-1)//'.fch'
 chkname = gjfname(1:i-1)//'.chk'
 inquire(file=TRIM(fchname), exist=alive)
 if(.not. alive) call formchk(chkname)

 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 allocate(mo_a(nbf,nif), source=0d0)
 if(wfn_type == 3) allocate(mo_b(nbf,nif), source=0d0)
 k1 = 0; k2 = 0; k5 = 0

 do i = 1, n, 1
  j = index(gjfname0(i), '.gjf', back=.true.)
  fchname0 = gjfname0(i)(1:j-1)//'.fch'
  call read_nbf_and_nif_from_fch(fchname0, nbf0, nif0)
  allocate(mo_a0(nbf0,nif0))
  if(pos(i)) then
   call read_na_and_nb_from_fch(fchname0, na0, nb0)
   call read_mo_from_fch(fchname0, nbf0, nif0, 'a', mo_a0)
  else ! negative spin
   call read_na_and_nb_from_fch(fchname0, nb0, na0)
   if(wfn_type0(i) == 3) then
    call read_mo_from_fch(fchname0, nbf0, nif0, 'b', mo_a0)
   else ! ROHF
    call read_mo_from_fch(fchname0, nbf0, nif0, 'a', mo_a0)
   end if
  end if

  if(k1+nbf0 > nbf) then
   write(6,'(A)') 'ERROR in subroutine direct_sum_frag_mo2super_mo: the 1st&
                 & dimension of array mo_a out of range!'
   write(6,'(5(A,I0))') 'k1=',k1,',nbf0=',nbf0,',nbf=',nbf,',i=',i,',n=',n
   stop
  end if
  mo_a(k1+1:k1+nbf0, k2+1:k2+na0) = mo_a0(:,1:na0)

  if(wfn_type == 3) then ! UHF
   if(wfn_type0(i) == 3) then ! UHF fragment
    allocate(mo_b0(nbf0,nif0))
    if(pos(i)) then
     call read_mo_from_fch(fchname0, nbf0, nif0, 'b', mo_b0)
    else ! negative spin
     call read_mo_from_fch(fchname0, nbf0, nif0, 'a', mo_b0)
    end if
    mo_b(k1+1:k1+nbf0, k5+1:k5+nb0) = mo_b0(:,1:nb0)
    deallocate(mo_b0)
   else ! RHF/ROHF fragment
    mo_b(k1+1:k1+nbf0, k5+1:k5+nb0) = mo_a0(:,1:nb0)
   end if
   k5 = k5 + nb0
  end if

  deallocate(mo_a0)
  k1 = k1 + nbf0
  k2 = k2 + na0
 end do ! for i

 ! In previous version, virtual MOs of the supermolecule were constructed from
 ! virtual MOs of fragments. But it did not work for linear dependence case,
 ! so I change it to construction of virtual MOs using PAO

 allocate(ovlp(nbf,nbf), mo_a0(nbf,nif))
 call get_ao_ovlp_using_fch(fchname, nbf, ovlp)
 ! orthonormalize occupied MOs
 call orthonormalize_orb(nbf, k2, ovlp, mo_a(:,1:k2), mo_a0(:,1:k2))
 ! generate virtual MOs using PAO
 call construct_vir(nbf, nif, k2+1, mo_a0, ovlp, mo_a)
 call write_mo_into_fch(fchname, nbf, nif, 'a', mo_a)
 deallocate(mo_a, mo_a0)

 if(wfn_type == 3) then
  allocate(mo_b0(nbf,nif))
  call orthonormalize_orb(nbf, k5, ovlp, mo_b(:,1:k5), mo_b0(:,1:k5))
  call construct_vir(nbf, nif, k5+1, mo_b0, ovlp, mo_b)
  call write_mo_into_fch(fchname, nbf, nif, 'b', mo_b)
  deallocate(mo_b, mo_b0)
 end if

 deallocate(ovlp)
 call unfchk(fchname, chkname)
end subroutine direct_sum_frag_mo2super_mo

! sum fragment Total SCF Densities and print the 
subroutine sum_frag_density_and_prt_into_fch(n, fname0, pos, fname)
 use util_wrapper, only: formchk
 implicit none
 integer :: i, j, nif, nbf, nbf1
 integer, intent(in) :: n
 character(len=240) :: chkname
 character(len=240), intent(in) :: fname0(n), fname
 character(len=240), allocatable :: fchname(:)
 real(kind=8), allocatable :: dm(:,:), dm1(:,:)
 logical :: alive
 logical, allocatable :: has_spin_density(:)
 logical, intent(in) :: pos(n) ! negative spin means to switch alpha/beta

 allocate(fchname(n+1))

 ! if fname0(n) and fname are .gjf files, convert filenames into .fch
 do i = 1, n, 1
  j = index(fname0(i), '.fch', back=.true.)
  if(j == 0) then
   j = index(fname0(i), '.gjf', back=.true.)
   fchname(i) = fname0(i)(1:j-1)//'.fch'
  else ! j > 0
   fchname(i) = fname0(i)
  end if
 end do ! for i

 j = index(fname, '.fch', back=.true.)
 if(j == 0) then
  j = index(fname, '.gjf', back=.true.)
  fchname(n+1) = fname(1:j-1)//'.fch'
 else ! j > 0
  fchname(n+1) = fname
  j = index(fname, '.fch', back=.true.)
 end if
 chkname = fname(1:j-1)//'.chk'
 inquire(file=TRIM(fchname(n+1)), exist=alive)
 if(.not. alive) call formchk(chkname, fchname(n+1))

 call read_nbf_and_nif_from_fch(fchname(n+1), nbf, nif)
 allocate(dm(nbf,nbf), dm1(nbf,nbf))
 dm = 0d0

 do i = 1, n, 1
  call read_nbf_and_nif_from_fch(fchname(i), nbf1, nif)
  if(nbf1 /= nbf) then
   write(6,'(A)') 'ERROR in subroutine sum_frag_density_and_prt_into_fch:&
                  & nbf1 /= nbf.'
   write(6,'(A,I0,A)') 'Inconsistent nbf between fragment ',i,' and the&
                       & total system.'
   write(6,'(2(A,I0))') 'nbf1=', nbf1, ', nbf=', nbf
   stop
  end if
  call read_density_from_fch(fchname(i), 1, nbf, dm1)
  dm = dm + dm1
 end do ! for i
 call write_density_into_fch(fchname(n+1), nbf, .true., dm)

 allocate(has_spin_density(n))
 has_spin_density = .false.
 dm = 0d0

 do i = 1, n, 1
  call detect_spin_scf_density_in_fch(fchname(i), has_spin_density(i))
  if(.not. has_spin_density(i)) cycle
  call read_density_from_fch(fchname(i), 2, nbf, dm1)
  if(pos(i)) then
   dm = dm + dm1
  else
   dm = dm - dm1 ! interchange alpha/beta spin, which is equal to -dm1
  end if
 end do ! for i

 call detect_spin_scf_density_in_fch(fchname(n+1), alive)

 if(ANY(has_spin_density .eqv. .true.)) then
  if(alive) then
   call write_density_into_fch(fchname(n+1), nbf, .false., dm)
  else ! no Spin SCF Density in the total system .fch file
   write(6,'(A)') 'ERROR in subroutine sum_frag_density_and_prt_into_fch: so&
                  &me fragment has UHF-type wave'
   write(6,'(A)') 'function, but the total system has RHF-type wave function.&
                  & Inconsistency detected.'
   stop
  end if
 else ! all fragments has RHF-type wave function
  if(alive) call write_density_into_fch(fchname(n+1), nbf, .false., dm)
 end if

 deallocate(fchname, dm, dm1, has_spin_density)
end subroutine sum_frag_density_and_prt_into_fch

! modify 'guess(only,save)' in a given .gjf file
subroutine modify_guess_only_in_gjf(gjfname)
 implicit none
 integer :: i, j, fid1, fid2, RENAME
 character(len=240) :: buf, fname
 character(len=240), intent(in) :: gjfname

 fname = TRIM(gjfname)//'.t'
 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fname),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:1) == '#') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')

 i = index(buf, '/')
 j = index(buf(i+1:),' ') + i - 1
 buf = buf(1:i-1)//' chkbasis'//buf(j+1:)

 i = index(buf, 'guess')
 buf(i:) = 'guess=read geom=allcheck'

 ! assuming '#p R...' or '#p U...'
 if(buf(4:4)=='u' .or. buf(4:4)=='U') then
  buf = TRIM(buf)//' stable=opt'
  i = index(buf, ',conver=6')
  if(i == 0) i = index(buf, ',conver=7')
  if(i > 0) buf = buf(1:i-1)//buf(i+9:)
 end if

 write(fid2,'(A,/)') TRIM(buf)
 close(fid2)
 i = RENAME(TRIM(fname), TRIM(gjfname))
end subroutine modify_guess_only_in_gjf

