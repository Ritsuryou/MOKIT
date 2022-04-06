! written by jxzou at 20210612: transform basis sets in GAMESS format to Dalton format
!  This file is initially copied from file bas_gms2molcas. Modifications are made

! Note: Currently isotopes are not tested.
program main
 use pg, only: iout
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format
 logical :: spherical

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A)') 'Example1: bas_gms2dal a.inp (generate an a.input file)'
  write(iout,'(A,/)') "Example2: bas_gms2dal a.inp -sph (without 'Cartesian')"
  stop
 end if

 str = ' '
 fname = ' '
 spherical = .false.
 call getarg(1, fname)
 call require_file_exist(fname)

 if(i == 2) then
  call getarg(2,str)

  if(str == '-sph') then
   spherical = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine bas_gms2dal: wrong command line arguments!'
   write(iout,'(A)') "The 2nd argument can only be '-sph'. But got '"//str//"'"
   stop
  end if
 end if

 call bas_gms2dal(fname, spherical)
 stop
end program main

! Transform the basis sets in GAMESS format to those in Dalton format
subroutine bas_gms2dal(fort7, spherical)
 use pg, only: iout, natom, ram, ntimes, elem, coor, highest, all_ecp, ecp_exist
 implicit none
 integer :: i, nline, rc, charge, mult, fid1, fid2
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, dalfile, molfile
 character(len=1) :: stype
 logical :: uhf
 logical, intent(in) :: spherical

 buf = ' '   ! initialization
 i = index(fort7, '.', back=.true.)
 dalfile = fort7(1:i-1)//'.dal'
 molfile = fort7(1:i-1)//'.mol'

 call read_natom_from_gms_inp(fort7, natom)
 allocate(ram(natom), elem(natom), coor(3,natom), ntimes(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, ram, coor)
 ! ram cannot be deallocated here since subroutine prt_prim_gau will use it

 call calc_ntimes(natom, elem, ntimes)
 call read_charge_and_mult_from_gms_inp(fort7, charge, mult, uhf, ecp_exist)
 if(uhf) then
  write(iout,'(/,A)') 'WARNING in subroutine bas_gms2dal: Dalton does not support UHF.'
  write(iout,'(A)') 'Basis set data will still be written.'
 end if
 if(ecp_exist) call create_dir('ecp_data')

 open(newunit=fid1,file=TRIM(dalfile),status='replace')
 write(fid1,'(A)') '**DALTON INPUT'
 write(fid1,'(A)') '.RUN WAVE FUNCTIONS'
 write(fid1,'(A)') '**WAVE FUNCTIONS'
 write(fid1,'(A)') '.HF'
 if(mult > 1) then
  write(fid1,'(A)') '*SCF INPUT'
  write(fid1,'(A)') '.DOUBLY OCCUPIED'
  write(fid1,'(I0)') (SUM(ram)-charge-mult+1)/2
  write(fid1,'(A)') '.SINGLY OCCUPIED'
  write(fid1,'(I0)') mult-1
 end if
 write(fid1,'(A)') '*ORBITAL INPUT'
 write(fid1,'(A)') '.PUNCHOUTPUTORBITALS'
 write(fid1,'(A)') '**END OF INPUT'
 close(fid1)

 call read_all_ecp_from_gms_inp(fort7)

 ! find the $DATA section
 open(newunit=fid1,file=TRIM(fort7),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)',iostat=rc) buf
  if(rc /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:6))
   if(buf(3:6) == 'DATA') exit
  end if
 end do ! for while

 if(rc /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: No $DATA section found&
                   & in file '//TRIM(fort7)//'.'
  close(fid1)
  stop
 end if

 ! skip 2 lines: the Title line and the Point Group line
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()

 open(newunit=fid2,file=TRIM(molfile),status='replace')
 write(fid2,'(A)') 'ATOMBASIS'
 write(fid2,'(A)') 'generated by utility bas_gms2dal in MOKIT'
 write(fid2,'(A)') 'Basis set specified with gen'
 write(fid2,'(2(A,I0),A)',advance='no') 'AtomTypes=', natom, &
  ' Integrals=1.0D-14 Charge=',charge,' NoSymmetry Angstrom'
 if(spherical) then
  write(fid2,'(/)',advance='no')
 else
  write(fid2,'(A)') ' Cartesian'
 end if

 do i = 1, natom, 1
  read(fid1,'(A)',iostat=rc) buf
  if(rc /= 0) exit
  ! 'buf' contains the element, ram and coordinates
  if(elem(i) == 'Bq') then
   write(fid2,'(A)') 'Charge=0. Atoms=1 Basis=INTGRL Ghost'
   write(fid2,'(A2,3(1X,F15.8))') elem(i), coor(1:3,i)
   cycle
  end if

  ! deal with primitive gaussians
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
   call read_prim_gau(stype, nline, fid1)
  end do ! for while

  call get_highest_am()
  write(fid2,'(3(A,I0))',advance='no') 'Charge=',ram(i),'. Atoms=1 Basis=INTGRL&
   & Blocks=',highest+1
  write(fid2,'(A)',advance='no') REPEAT(' 1', highest+1)
  if(all_ecp(i)%ecp) then
   write(fid2,'(A,I0,A)') ' ECP=',i,'_ecp'
  else
   write(fid2,'(/)',advance='no')
  end if
  write(fid2,'(A2,3(1X,F15.8))') elem(i), coor(1:3,i)

  ! print basis sets and ECP/PP (if any) of this atom in Dalton format
  call prt_prim_gau_dalton(i, fid2)

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 close(fid1)
 ! now ram can be deallocated
 deallocate(ram, elem, ntimes, coor, all_ecp)

 if(rc /= 0) then
  write(iout,'(A)') "ERROR in subroutine bas_gms2dal: it seems the '$DATA'&
                   & has no corresponding '$END'."
  write(iout,'(A)') 'Incomplete file '//TRIM(fort7)
  close(fid2,status='delete')
  stop
 end if

 close(fid2)
 return
end subroutine bas_gms2dal

! print primitive gaussians in Dalton format
subroutine prt_prim_gau_dalton(iatom, fid)
 use pg, only: prim_gau, all_ecp, ecp_exist, ram
 implicit none
 integer :: i, j, k, m, n, nline, ncol, ecpid
 integer, intent(in) :: iatom, fid
 character(len=240) :: ecpname

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(A1,I4,I5)') 'H', nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(F20.10)',advance='no') prim_gau(i)%coeff(j,1)
   do k = 2, ncol
    write(fid,'(F20.10)',advance='no') prim_gau(i)%coeff(j,k)
    if(MOD(k-1,3)==0 .and. k<ncol) write(fid,'(/,20X)',advance='no')
   end do ! for k
   write(fid,'(/)',advance='no')
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(iatom)%ecp) then
  write(ecpname,'(A,I0,A)') 'ecp_data/',iatom,'_ecp'
  open(newunit=ecpid,file=TRIM(ecpname),status='replace')
  write(ecpid,'(A,/,A,I4,/,A)') '$','a', ram(iatom), '$'
  m = all_ecp(iatom)%highest
  write(ecpid,'(2I4,I0)') m, all_ecp(iatom)%core_e

  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   write(ecpid,'(I12)') n
   do j = 1, n, 1
    write(ecpid,'(I2,2F20.12)') all_ecp(iatom)%potential(i)%col2(j),&
     all_ecp(iatom)%potential(i)%col3(j), all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
  write(ecpid,'(A)') '$'
 end if

 return
end subroutine prt_prim_gau_dalton

