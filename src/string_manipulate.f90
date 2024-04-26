! written by jxzou at 20201208: move string manipulation subroutines into this file

module phys_cons ! physics constants
 implicit none
 real(kind=8), parameter :: au2ev = 27.211396d0
 real(kind=8), parameter :: au2kcal = 627.51d0 ! a.u. to kcal/mol
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), parameter :: GPa2eVA3 = 6.241509074460764d-3
end module phys_cons

! transform a string into upper case
subroutine upper(buf)
 implicit none
 integer :: i, k
 character(len=*), intent(inout) :: buf

 do i = 1, LEN(buf), 1
  k = IACHAR(buf(i:i))
  if(k>=97 .and. k<=122) buf(i:i) = ACHAR(k-32)
 end do
end subroutine upper

! transform a string into lower case
subroutine lower(buf)
 implicit none
 integer :: i, j, k
 character(len=*), intent(inout) :: buf

 k = LEN_TRIM(buf)
 do i = 1, k, 1
  j = IACHAR(buf(i:i))
  if(j>=65 .and. j<=90) buf(i:i) = ACHAR(j+32)
 end do ! for i
end subroutine lower

! convert a (character) stype to (integer) itype
subroutine stype2itype(stype, itype)
 implicit none
 integer, intent(out) :: itype
 character(len=1), intent(in) :: stype

 ! 'S', 'P', 'D', 'F', 'G', 'H', 'I'
 !  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7
 select case(stype)
 case('S')
  itype = 1
 case('P')
  itype = 2
 case('D')
  itype = 3
 case('F')
  itype = 4
 case('G')
  itype = 5
 case('H')
  itype = 6
 case('I')
  itype = 7
 case('L') ! 'L' is 'SP'
  itype = 0
 case default
  write(6,'(A)') 'ERROR in subroutine stype2itype: stype out of range.'
  write(6,'(A)') 'stype= '//TRIM(stype)
  stop
 end select
end subroutine stype2itype

! check whether there exists DKH keywords in a given GAMESS .inp file
subroutine check_DKH_in_gms_inp(inpname, order)
 implicit none
 integer :: i, k, fid
 integer, intent(out) :: order
! -2: no DKH
! -1: RESC
!  0: DKH 0th-order
!  2: DKH2
!  4: DKH4 with SO
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 character(len=1200) :: longbuf

 longbuf = ' '
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  longbuf = TRIM(longbuf)//TRIM(buf)
  call upper(buf)
  if(INDEX(buf,'$END') /= 0) exit
 end do ! for while
 close(fid)

 call upper(longbuf)

 order = -2
 if(INDEX(longbuf,'RELWFN') == 0) then
  return
 else
  if(INDEX(longbuf,'RELWFN=DK') == 0) then
   write(6,'(A)') 'Warning in subroutine check_DKH_in_gms_inp: unsupported&
                    & relativistic method detected.'
   write(6,'(A)') '(Open)Molcas does not support RELWFN=LUT-IOTC, IOTC,&
                    & RESC, or NESC in GAMESS. Only RELWFN=DK is supported.'
   write(6,'(A)') 'The MO transferring will still be proceeded. But the result&
                    & may be non-sense.'
  end if
 end if

 order = 2
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf)
  if(INDEX(buf,'$RELWFN') /= 0) exit
 end do ! for while
 close(fid)

 if(i == 0) then
  k = INDEX(buf,'NORDER=')
  if(k /= 0) then
   read(buf(k+7:),*) order
  end if
 end if
end subroutine check_DKH_in_gms_inp

! check whether X2C appears in a given GAMESS .inp file
! Note: GAMESS does not support X2C, this is just for the utility bas_gms2molcas
!  to recognize X2C and pass it into (Open)Molcas .input file
subroutine check_X2C_in_gms_inp(inpname, X2C)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical, intent(out) :: X2C

 X2C = .false. ! default

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do i = 1, 5
  read(fid,'(A)') buf
  if(INDEX(buf,'X2C') /= 0) X2C = .true.
  if(INDEX(buf,'$END') /= 0) exit
 end do ! for i

 close(fid)
end subroutine check_X2C_in_gms_inp

subroutine check_sph_in_gjf(gjfname, sph)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=1200) :: longbuf
 character(len=240), intent(in) :: gjfname
 logical, intent(out) :: sph

 sph = .true.
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '#') exit
 end do ! for while

 longbuf = buf
 do i = 1, 5
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  longbuf = TRIM(longbuf)//TRIM(buf)
 end do ! for i

 close(fid)
 if(INDEX(longbuf,'6D')>0 .or. INDEX(longbuf,'6d')>0) sph = .false.
end subroutine check_sph_in_gjf

! convert a filename into which molpro requires, i.e. in lowercase
subroutine convert2molpro_fname(fname, suffix)
 implicit none
 integer :: i, len1, len2
 character(len=240), intent(inout) :: fname
 character(len=2), intent(in) :: suffix

 if(LEN_TRIM(fname) == 0) then
  write(6,'(/,A)') 'ERROR in subroutine convert2molpro_fname: input fname is NU&
                   &LL.'
  stop
 end if

 if(fname(1:1) == ' ') fname = ADJUSTL(fname)

 i = INDEX(fname, '.', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine convert2molpro_fname: '.' character not&
                   & found in"
  write(6,'(A)') 'filename '//TRIM(fname)
  stop
 end if

 len1 = INDEX(fname, '.', back=.true.) - 1
 if(len1 == -1) len1 = LEN_TRIM(fname)
 len2 = LEN(suffix)

 if(len1+len2 > 32) then
  fname = fname(1:32-len2)//suffix
 else
  fname = fname(1:len1)//suffix
 end if

 call lower(fname)
end subroutine convert2molpro_fname

! add DKH2 related keywords into a given GAMESS .inp file,
! and switch the default SOSCF into DIIS
subroutine add_DKH2_into_gms_inp(inpname)
 implicit none
 integer :: i, k, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do i = 1, 3
  read(fid1,'(A)') buf

  if(INDEX(buf, 'RELWFN=DK') /= 0) then
   close(fid1)
   close(fid2,status='delete')
   return
  end if

  k = INDEX(buf,'$END')
  if(k /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for i

 write(fid2,'(A)') buf(1:k-1)//' RELWFN=DK $END'

 do while(.true.)
  read(fid1,'(A)') buf
  k = INDEX(buf,'$END')
  if(k /= 0) exit
  if(INDEX(buf,'$DATA') /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(k == 0) then
  write(fid2,'(A)') '$SCF DIRSCF=.TRUE. DIIS=.T. SOSCF=.F. $END'
 else
  write(fid2,'(A)') buf(1:k-1)//' DIIS=.T. SOSCF=.F. $END'
 end if

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_DKH2_into_gms_inp

! add DKH2 keyword into a given Gaussian .fch(k) file
subroutine add_DKH2_into_fch(fchname)
 implicit none
 integer :: i, j, k, nline, nterm, fid, fid1, RENAME
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 character(len=1200) :: longbuf, longbuf1
 logical :: no_route, alive(3)

 buf = ' '; longbuf = ' '; nterm = 0
 i = INDEX(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')
 no_route = .false.

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Charge') then
   no_route = .true.
   exit
  end if

  if(buf(1:5) == 'Route') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(no_route) then
  write(fid1,'(A5,38X,A)') 'Route','C   N=           3'
  write(fid1,'(A)') '#p int(nobasistransform,DKH2) nosymm'
 else
  if(i /= 0) then
   write(6,'(A)') 'ERROR in subroutine add_DKH2_into_fch: incomplete .fch file.'
   write(6,'(A)') "Neither 'Route' nor 'Charge' is detected in file "//TRIM(fchname)
   close(fid)
   close(fid1,status='delete')
   stop
  else
   k = INDEX(buf, '='); read(buf(k+1:),*) nterm
   do while(.true.)
    read(fid,'(A)',iostat=i) buf
    if(i /= 0) exit
    if(buf(1:6) == 'Charge') exit
    longbuf = TRIM(longbuf)//TRIM(buf)
   end do ! for while

   if(i /= 0) then
    write(6,'(A)') 'ERROR in subroutine add_DKH2_into_fch: incomplete .fch file.'
    write(6,'(A)') "No 'Charge' is detected in file "//TRIM(fchname)
    close(fid)
    close(fid1,status='delete')
    stop
   else
    longbuf1 = longbuf
    call upper(longbuf1)
    alive = [(INDEX(longbuf1,'DKH2')/=0),(INDEX(longbuf1,'DOUGLASKROLLHESS')/=0),&
             (INDEX(longbuf1,'DKH')/=0 .and. INDEX(longbuf1,'DKH4')==0 .and. &
              INDEX(longbuf1,'NODKH')==0 .and. INDEX(longbuf1,'DKHSO')==0)]
    if(ALL(alive .eqv. .false.)) then
     nterm = nterm + 1
     longbuf = TRIM(longbuf)//' int=DKH2'
    end if
    write(fid1,'(A5,38X,A,I2)') 'Route','C   N=          ', nterm
    k = LEN_TRIM(longbuf)
    nline = k/60
    if(k-60*nline > 0) nline = nline + 1
    do i = 1, nline, 1
     j = min(60*i,k)
     write(fid1,'(A)') longbuf(60*i-59:j)
    end do ! for i
   end if
  end if
 end if

 ! copy remaining content
 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine add_DKH2_into_fch

! add X2C keyword into a given Gaussian .fch(k) file
! Note: Obviously, Gaussian cannot use X2C. This is just for other utilities to
!  recognize the X2C keyword, so that other utilities can add related keywords
!  when generating input files of other programs
subroutine add_X2C_into_fch(fchname)
 implicit none
 integer :: i, j, k, nline, nterm, fid, fid1, RENAME
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 character(len=1200) :: longbuf, longbuf1
 logical :: no_route

 buf = ' '; longbuf = ' '; nterm = 0
 i = INDEX(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')
 no_route = .false.

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Charge') then
   no_route = .true.
   exit
  end if

  if(buf(1:5) == 'Route') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(no_route) then
  write(fid1,'(A5,38X,A)') 'Route','C   N=           3'
  write(fid1,'(A)') '#p int(nobasistransform,X2C) nosymm'
 else
  if(i /= 0) then
   write(6,'(A)') 'ERROR in subroutine add_X2C_into_fch: incomplete .fch file.'
   write(6,'(A)') "Neither 'Route' nor 'Charge' is detected in file "//TRIM(fchname)
   close(fid)
   close(fid1,status='delete')
   stop
  else
   k = INDEX(buf, '='); read(buf(k+1:),*) nterm
   do while(.true.)
    read(fid,'(A)',iostat=i) buf
    if(i /= 0) exit
    if(buf(1:6) == 'Charge') exit
    longbuf = TRIM(longbuf)//TRIM(buf)
   end do ! for while

   if(i /= 0) then
    write(6,'(A)') 'ERROR in subroutine add_X2C_into_fch: incomplete .fch file.'
    write(6,'(A)') "No 'Charge' is detected in file "//TRIM(fchname)
    close(fid)
    close(fid1,status='delete')
    stop
   else
    longbuf1 = longbuf
    call upper(longbuf1)
    j = INDEX(longbuf1, 'DKH'); k = INDEX(longbuf1, 'NODKH')
    if(j/=0 .and. k==0) then
     longbuf(j:j+2) = 'X2C'
    else
     nterm = nterm + 1
     longbuf = TRIM(longbuf)//' int=X2C'
    end if
    write(fid1,'(A5,38X,A,I2)') 'Route','C   N=          ', nterm
    k = LEN_TRIM(longbuf)
    nline = k/60
    if(k-60*nline > 0) nline = nline + 1
    do i = 1, nline, 1
     j = min(60*i,k)
     write(fid1,'(A)') longbuf(60*i-59:j)
    end do ! for i
   end if
  end if
 end if

 ! copy remaining content
 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine add_X2C_into_fch

! detect the number of columns of data in a string buf
function detect_ncol_in_buf(buf) result(ncol)
 implicit none
 integer :: i, ncol
 character(len=24), allocatable :: sbuf(:)
 character(len=240), intent(in) :: buf

 if(LEN_TRIM(buf) == 0) then
  ncol = 0
  return
 end if

 ncol = 1
 do while(.true.)
  allocate(sbuf(ncol))
  read(buf,*,iostat=i) sbuf(1:ncol)
  deallocate(sbuf)
  if(i /= 0) exit
  ncol = ncol + 1
 end do ! for while

 ncol = ncol - 1
end function detect_ncol_in_buf

! modify the memory in a given .inp file
subroutine modify_memory_in_gms_inp(inpname, mem, nproc)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: mem, nproc
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'MWORDS') /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine modify_memory_in_gms_inp: no 'MWORDS' fou&
                 &nd in file "//TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',FLOOR(DBLE(mem)*1d3/(8d0*DBLE(nproc))),' $END'

 ! copy the remaining content
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_memory_in_gms_inp

! modify memory in a given PSI4 input file
! Note: input mem is in unit GB
subroutine modify_memory_in_psi4_inp(inpname, mem)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: mem
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf)
 read(fid,'(A)') buf
 write(fid1,'(A,I0,A)') 'memory ', mem, ' GB'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_memory_in_psi4_inp

! modify memory in a given Q-Chem input file
! Note: input mem is in unit GB
subroutine modify_memory_in_qchem_inp(mem, inpname)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: mem
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 i = INDEX(inpname, '.in', back=.true.)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == 'mem_tot') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A,I0)') 'mem_total ',mem*1000 ! in MB

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_memory_in_qchem_inp

! add given/specified RIJK basis set into a PSI4 input file
subroutine add_RIJK_bas_into_psi4_inp(inpname, RIJK_bas)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: inpname1
 character(len=21), intent(in) :: RIJK_bas
 character(len=240), intent(in) :: inpname

 if(LEN_TRIM(RIJK_bas) == 0) then
  write(6,'(A)') 'ERROR in subroutine add_RIJK_bas_into_psi4_inp:'
  stop
 end if

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_RIJK_bas_into_psi4_inp

! add given/specified RIJK basis set into an ORCA input file
subroutine add_RIJK_bas_into_orca_inp(inpname, RIJK_bas)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: inpname1
 character(len=21), intent(in) :: RIJK_bas
 character(len=240), intent(in) :: inpname

 if(LEN_TRIM(RIJK_bas) == 0) then
  write(6,'(A)') 'ERROR in subroutine add_RIJK_bas_into_orca_inp: input RI&
                   & basis set is null string.'
  stop
 end if
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_RIJK_bas_into_orca_inp

! detect whether there exists the charge keyword in a given .gjf file
function detect_charge_key_in_gjf(gjfname) result(has_charge)
 implicit none
 integer :: i, j, nblank, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 logical :: has_charge

 has_charge = .true.; nblank = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '#') then
   call lower(buf)
   if(INDEX(buf,'charge') > 0) then
    close(fid)
    return
   end if
  end if
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 1) exit
 end do ! for while

 read(fid,'(A)') buf
 close(fid)
 i = INDEX(buf,'{')
 j = INDEX(buf,'}')
 if(i>0 .and. j>0) then
  call lower(buf(i+1:j-1))
  if(INDEX(buf(i+1:j-1),'charge') > 0) return
 end if
 has_charge = .false.
end function detect_charge_key_in_gjf

! copy mixed/user-defined basis set in a given .gjf file to a .bas file
subroutine record_gen_basis_in_gjf(gjfname, basname, add_path)
 implicit none
 integer :: i, nblank0, nblank, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 character(len=240), intent(out) :: basname
 logical, intent(in) :: add_path
 logical, external :: detect_charge_key_in_gjf
 logical :: nobasis

 if(detect_charge_key_in_gjf(gjfname)) then
  nblank0 = 4
 else
  nblank0 = 3
 end if
 i = INDEX(gjfname, '.gjf', back=.true.)
 basname = gjfname(1:i-1)//'.bas'

 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == nblank0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine record_gen_basis_in_gjf: incomplete fil&
                   &e: '//TRIM(gjfname)
  close(fid1)
  stop
 end if

 read(fid1,'(A)',iostat=i) buf
 nobasis = .false.
 if(i /= 0) nobasis = .true.
 if((.not.nobasis) .and. LEN_TRIM(buf)==0) nobasis = .true.

 if(nobasis) then
  write(6,'(/,A)') 'ERROR in subroutine record_gen_basis_in_gjf: no mixed/user-&
                   &defined basis'
  write(6,'(A)') 'set detected in file '//TRIM(gjfname)
  close(fid1)
  stop
 end if

 buf = ADJUSTL(buf)
 if(buf(1:1) == '-') then
  i = IACHAR(buf(2:2))
 else
  i = IACHAR(buf(1:1))
 end if

 if(.not. ((i>96 .and. i<123) .or. (i>64 .and. i<91))) then
  write(6,'(A)') 'ERROR in subroutine record_gen_basis_in_gjf: the first charac&
                 &ter in mixed/user-defined'
  write(6,'(A)') 'basis set is neither a-z, nor A-Z. This is not an element sym&
                 &bol. This format of basis'
  write(6,'(A)') 'set cannot be recognized by automr. Problematic file: '//&
                  TRIM(gjfname)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(basname),status='replace')
 write(fid2,'(A)') TRIM(buf)

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  ! we cannot use LEN_TRIM(buf)==0 to judge since there may exist ECP
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(/)')
 close(fid1)
 close(fid2)

 call add_hyphen_for_elem_in_basfile(basname)
 if(add_path) call add_mokit_path_to_genbas(basname)
end subroutine record_gen_basis_in_gjf

! add '-' symbol before elements, in a .bas file
subroutine add_hyphen_for_elem_in_basfile(basname)
 implicit none
 integer :: i, j, nbat1, nbat2, fid, fid1, RENAME
 character(len=7) :: str
 character(len=240) :: buf0, buf, basname1
 character(len=240), intent(in) :: basname

 basname1 = TRIM(basname)//'.t'
 open(newunit=fid,file=TRIM(basname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(basname1),status='replace')
 buf0 = '****'

 ! deal with the basis set data
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '!') then
   write(fid1,'(A)') TRIM(buf)
   cycle
  end if

  if(buf0(1:4) == '****') then
   if(LEN_TRIM(buf) == 0) exit
   call add_hyphen_for_elem_in_buf(buf)
   write(fid1,'(A)') TRIM(buf)
  else
   write(fid1,'(A)') TRIM(buf)
  end if

  buf0 = buf
 end do ! for i

 write(fid1,'(/)',advance='no')

 ! deal with the pseudo potential data
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit

  call add_hyphen_for_elem_in_buf(buf)
  write(fid1,'(A)') TRIM(buf)

  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  read(buf,*,iostat=i) str, nbat1

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
 end do ! for while

 write(fid1,'(/)')
 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(basname1), TRIM(basname))
end subroutine add_hyphen_for_elem_in_basfile

! delete '-' symbol before elements, in a .bas file
subroutine del_hyphen_for_elem_in_basfile(basname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: basname1
 character(len=240), intent(in) :: basname
 character(len=300) :: buf

 basname1 = TRIM(basname)//'.t'
 open(newunit=fid,file=TRIM(basname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(basname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  select case(buf(1:1))
  case('!')
   write(fid1,'(A)') TRIM(buf)
   cycle
  case('-')
   i = IACHAR(buf(2:2))
   if((i>64 .and. i<91) .or. (i>96 .and. i<123)) buf = TRIM(buf(2:))
  end select

  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(basname1), TRIM(basname))
end subroutine del_hyphen_for_elem_in_basfile

function get_mokit_root() result(mokit_root)
 implicit none
 integer :: i, fid
 character(len=240) :: home, mokit_root !, buf

 mokit_root = ' '
 call getenv('MOKIT_ROOT', mokit_root)

 if (len_trim(mokit_root) < 1) then
  call getenv('HOME', home)
  open(newunit=fid,file=TRIM(home)//'/.mokitrc',status='old',position='rewind')
  read(fid,'(A)',iostat=i) mokit_root
  if (len_trim(mokit_root) < 1) then
   write(6,'(/,A)') 'ERROR in subroutine get_mokit_root: invalid MOKIT_ROOT'
   stop
  end if
  close(fid)
 end if
end function get_mokit_root

! add MOKIT_ROOT path into basis sets like ANO-RCC-VDZP, DKH-def2-SVP in file
! basname because MOKIT has these basis sets in $MOKIT_ROOT/mokit/basis/
subroutine add_mokit_path_to_genbas(basname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=11) :: sbuf = ' '
 character(len=240) :: mokit_root, basname1
 character(len=240), external :: get_mokit_root
 character(len=240), intent(in) :: basname
 character(len=480) :: buf = ' '
 logical :: alive(7)

 sbuf = ' '
 !mokit_root = ' '
 !call getenv('MOKIT_ROOT', mokit_root)
 mokit_root = get_mokit_root()
 basname1 = TRIM(basname)//'.t'

 open(newunit=fid,file=TRIM(basname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(basname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  sbuf = TRIM(buf(1:11))
  
  call upper(sbuf)
  alive = [(sbuf(1:3)=='X2C'), (sbuf(1:6)=='PCSSEG'), (sbuf(1:5)=='ANO-R'), &
           (sbuf(1:8)=='DKH-DEF2'), (sbuf(9:11)=='X2C'), (sbuf=='MA-DKH-DEF2'),&
           (sbuf(8:11)=='-F12')]
  if(ANY(alive .eqv. .true.)) then
   buf = '@'//TRIM(mokit_root)//'/mokit/basis/'//TRIM(buf)//'/N'
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(basname1), TRIM(basname))
end subroutine add_mokit_path_to_genbas

subroutine create_basfile(basfile, basis)
 implicit none
 integer :: fid
 character(len=240), intent(in) :: basfile
 character(len=*), intent(in) :: basis

 open(newunit=fid,file=TRIM(basfile),status='replace')
 write(fid,'(A)') TRIM(basis)
 write(fid,'(/)')
 close(fid)
 call add_mokit_path_to_genbas(basfile)
end subroutine create_basfile

! add '-' symbol for each element, in a buf
subroutine add_hyphen_for_elem_in_buf(buf)
 implicit none
 integer :: i, j
 integer, parameter :: max_nelem = 21
 character(len=3) :: str(max_nelem)
 character(len=240), intent(inout) :: buf

 do i = 1, max_nelem, 1
  read(buf,*,iostat=j) str(1:i)
  if(j /= 0) exit
 end do ! for i

 i = i - 1
 if(TRIM(str(i)) == '0') i = i - 1
 forall(j=1:i, str(j)(1:1)/='-') str(j) = '-'//TRIM(str(j))

 buf = TRIM(str(1))
 do j = 2, i, 1
  buf = TRIM(buf)//' '//TRIM(str(j))
 end do ! for i

 buf = TRIM(buf)//' 0'
end subroutine add_hyphen_for_elem_in_buf

! copy mixed/user-defined basis set content from file basname to gjfname
subroutine copy_gen_basis_bas2gjf(basname, gjfname)
 implicit none
 integer :: i, nblank, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: basname, gjfname

 if(LEN_TRIM(basname) == 0) return

 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 3) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine copy_gen_basis_bas2gjf: incomplete file&
                   & '//TRIM(gjfname)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(basname),status='old',position='rewind')
 do while(.true.)
  read(fid2,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid1)
 close(fid2)
end subroutine copy_gen_basis_bas2gjf

! read the version of dispersion correction from a .gjf file
subroutine read_disp_ver_from_gjf(gjfname, itype)
 implicit none
 integer :: i, fid
 integer, intent(out) :: itype
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 itype = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_disp_ver_from_gjf: no '#' symbol&
                   & found in file "//TRIM(gjfname)
  stop
 end if

 call lower(buf)
 if(INDEX(buf,'em=gd3bj')>0 .or. INDEX(buf,'empiricaldispersion=gd3bj')>0) then
  itype = 2
 else if(INDEX(buf,'em=gd3')>0 .or. INDEX(buf,'empiricaldispersion=gd3')>0) then
  itype = 1
 end if
end subroutine read_disp_ver_from_gjf

! print Fock operator coupling coefficients for ROGVB when nopen>=3
subroutine prt_gvb_couple_coeff(fid, ndb, nopen)
 implicit none
 integer :: i, ia
 integer, intent(in) :: fid, ndb, nopen
 character(len=3), allocatable :: f(:), alpha(:)
 character(len=4), allocatable :: beta(:)

 if(nopen < 3) then
  write(6,'(/,A)') 'ERROR in subroutine prt_gvb_couple_coeff: nopen<3.'
  close(fid)
  stop
 end if

 allocate(f(nopen))
 f = '0.5'
 if(ndb > 0) then
  ia = nopen*(nopen+3)/2
 else
  ia = nopen*(nopen+1)/2
 end if
 allocate(alpha(ia))
 alpha = '0.5'
 if(ndb > 0) then
  forall(i = 1:nopen) alpha(i*(i+1)/2) = '1.0'
 end if
 allocate(beta(ia))
 beta = '-0.5'

 if(ndb > 0) then
  write(fid,'(A,17(A1,A3))') '  F(1)=1.0', (',',f(i),i=1,nopen)
  write(fid,'(A,15(A1,A3))') '  ALPHA(1)=2.0', (',',alpha(i),i=1,ia)
  write(fid,'(A,12(A1,A4))') '  BETA(1)=-1.0', (',',beta(i),i=1,ia)
 else
  write(fid,'(A,17(A3,A1))') '  F(1)=', (f(i),',',i=1,nopen)
  write(fid,'(A,16(A3,A1))') '  ALPHA(1)=', (alpha(i),',',i=1,ia)
  write(fid,'(A,13(A4,A1))') '  BETA(1)=', (beta(i),',',i=1,ia)
 end if

 deallocate(f, alpha, beta)
end subroutine prt_gvb_couple_coeff

! copy GVB CI coefficients (or called pair coefficients) from a .dat file into
! another one 
subroutine copy_and_add_pair_coeff(addH_dat, datname, nopen)
 implicit none
 integer :: i, j, npair, fid1, fid2, fid3, RENAME
 integer, intent(in) :: nopen
 character(len=240) :: buf, new_dat
 character(len=240), intent(in) :: addH_dat, datname

 i = INDEX(addH_dat, '.dat', back=.true.)
 new_dat = addH_dat(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(addH_dat),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:6) == '$DATA') exit
 end do ! for while

 open(newunit=fid2,file=TRIM(new_dat),status='replace')
 write(fid2,'(A)') ' $DATA'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:5) == '$VEC') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 open(newunit=fid3,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid3,'(A)') buf
  if(buf(2:5) == '$SCF') exit
 end do ! for while

 i = INDEX(buf, 'CICOEF')
 if(i > 0) then
  write(fid2,'(A)') TRIM(buf)
  npair = 1
 else
  npair = 0
 end if

 do while(.true.)
  read(fid3,'(A)') buf
  i = INDEX(buf, 'CICOEF')
  if(i > 0) then
   j = INDEX(buf, '$END')
   if(j > 0) then
    buf(j:j+3) = '    '
    write(fid2,'(A)') TRIM(buf)
    exit
   else
    write(fid2,'(A)') TRIM(buf)
   end if
   npair = npair + 1
  else
   exit
  end if
 end do ! for while

 close(fid3)
 do i = 1, nopen, 1
  j = 2*(npair+i) - 1
  write(fid2,'(3X,A,I3,A)') 'CICOEF(',j,')= 0.7071067811865476,-0.7071067811865476'
 end do ! for i
 write(fid2,'(A,/,A)') ' $END', ' $VEC'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(new_dat), TRIM(addH_dat))
end subroutine copy_and_add_pair_coeff

! add the force keyword into a PySCF input file
subroutine add_force_key2py_script(mem, pyname, ccsd_t)
 implicit none
 integer :: fid
 integer, intent(in) :: mem ! GB
 character(len=240), intent(in) :: pyname
 logical, intent(in) :: ccsd_t

 open(newunit=fid,file=TRIM(pyname),status='old',position='append')
 if(ccsd_t) then
  write(fid,'(A)') 'from pyscf.grad import ccsd_t as ccsd_t_grad'
  write(fid,'(A)') 'mcg = ccsd_t_grad.Gradients(mc)'
 else
  write(fid,'(A)') 'from pyscf import grad'
  write(fid,'(A)') 'mcg = mc.Gradients()'
 end if

 write(fid,'(A,I0,A)') 'mcg.max_memory = ',mem*1000,' # MB'
 write(fid,'(A)') 'mcg.kernel()'
 close(fid)
end subroutine add_force_key2py_script

! standardize a set of elements, e.g. he -> He
subroutine standardize_elem(natom, elem)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 character(len=2), intent(inout) :: elem(natom)

 do i = 1, natom, 1
  elem(i) = ADJUSTL(elem(i))

  j = IACHAR(elem(i)(1:1))
  if(j>96 .and. j<123) elem(i)(1:1) = ACHAR(j-32)

  j = IACHAR(elem(i)(2:2))
  if(j>64 .and. j<91) elem(i)(2:2) = ACHAR(j+32)
 end do ! for i
end subroutine standardize_elem

! add the force keyword into a GAMESS input file
subroutine add_force_key2gms_inp(inpname)
 implicit none
 integer :: i, j, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do i = 1, 2
  read(fid,'(A)') buf
  j = INDEX(buf,'RUNTYP=ENERGY')
  if(j > 0) buf = buf(1:j+6)//'GRADIENT'//TRIM(buf(j+13:))
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_force_key2gms_inp

! add the force keyword into an ORCA input file
subroutine add_force_key2orca_inp(inpname)
 implicit none
 integer :: i, j, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do i = 1, 5
  read(fid,'(A)') buf
  j = INDEX(buf,'TightSCF')
  if(j > 0) buf = buf(1:j+7)//' EnGrad'//TRIM(buf(j+8:))
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_force_key2orca_inp

! read the basis set name from the 2nd line of a specified .fch file
subroutine read_basis_name_from_fch(fchname, basname)
 implicit none
 integer :: i, k, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 character(len=240), intent(out) :: basname

 basname = ' '
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)

 buf = TRIM(buf)
 k = LEN_TRIM(buf)
 i = INDEX(buf(1:k), ' ', back=.true.)
 basname = TRIM(buf(i+1:))
end subroutine read_basis_name_from_fch

subroutine find_specified_suffix(fname, suffix, i)
 implicit none
 integer, intent(out) :: i
 character(len=*), intent(in) :: suffix
 character(len=240), intent(in) :: fname

 i = INDEX(fname, suffix, back=.true.)

 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine find_specified_suffix: suffix '"//&
                   suffix//"' not found"
  write(6,'(A)') 'in filename '//TRIM(fname)
  stop
 end if
end subroutine find_specified_suffix

subroutine strip_ip_ea_eom(method)
 implicit none
 character(len=15), intent(inout) :: method

 call upper(method)
 select case(method(1:7))
 case('IP-EOM-','EOM-IP-','EA-EOM-','EOM-EA-')
  method = method(8:)
  return
 end select

 if(method(1:4) == 'EOM-') then
  method = method(5:)
 else if(method(1:3) == 'EOM') then
  method = method(4:)
 end if
end subroutine strip_ip_ea_eom

! calculate the integer array shell_coor ('Coordinates of each shell' in fch)
subroutine calc_shell_coor(ncontr, natom, shl2atm, coor, shell_coor)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, j
 integer, intent(in) :: ncontr, natom
 integer, intent(in) :: shl2atm(ncontr)
 real(kind=8), intent(in) :: coor(3,natom)
 real(kind=8), intent(out) :: shell_coor(3*ncontr)

 shell_coor = 0d0
 do i = 1, ncontr, 1
  j = shl2atm(i)
  shell_coor(3*i-2:3*i) = coor(:,j)/Bohr_const
 end do ! for i
end subroutine calc_shell_coor

! read the array size of shell_type and shell_to_atom_map from a given .fch(k) file
subroutine read_ncontr_from_fch(fchname, ncontr)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ncontr
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 ncontr = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:18) == 'Number of contract') exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_ncontr_from_fch: missing&
                 & 'Number of contract' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, ncontr
 close(fid)
end subroutine read_ncontr_from_fch

! read shell_type and shell_to_atom_map from a given .fch(k) file
subroutine read_shltyp_and_shl2atm_from_fch(fchname, k, shltyp, shl2atm)
 implicit none
 integer :: i, fid
 integer, intent(in) :: k
 integer, intent(out) :: shltyp(k), shl2atm(k)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find and read Shell types
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Shell types') exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_shltyp_and_shl2atm_from_fch:&
                 & missing 'Shell types' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 shltyp = 0
 read(fid,'(6(6X,I6))') (shltyp(i),i=1,k)
 ! read Shell types done

 ! find and read Shell to atom map
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Shell to atom') exit
 end do
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_shltyp_and_shl2atm_from_fch:&
                 & missing 'Shell to atom map' section in file "//TRIM(fchname)
  close(fid)
  return
 end if

 shl2atm = 0
 read(fid,'(6(6X,I6))') (shl2atm(i),i=1,k)
 close(fid)
end subroutine read_shltyp_and_shl2atm_from_fch

! replace Cartesian coordinates in fch
subroutine replace_coor_in_fch(fchname, natom, coor)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, k, ncontr, nline, fid, fid1, RENAME
 integer, intent(in) :: natom
 integer, allocatable :: shltyp(:), shl2atm(:)
 real(kind=8), intent(in) :: coor(3,natom) ! in Angstrom
 real(kind=8), allocatable :: shell_coor(:)
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shltyp(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shltyp, shl2atm)
 deallocate(shltyp)
 allocate(shell_coor(3*ncontr))
 call calc_shell_coor(ncontr, natom, shl2atm, coor, shell_coor)
 deallocate(shl2atm)

 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:12) == 'Current cart') exit
 end do ! for while
 write(fid1,'(5(1X,ES15.8))') coor/Bohr_const

 ! skip 'Current cartesian coordinates' in the old file
 k = 3*natom
 nline = k/5
 if(k - 5*nline > 0) nline = nline + 1
 do i = 1, nline, 1
  read(fid,'(A)') buf
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:12) == 'Coordinates ') exit
 end do ! for while
 write(fid1,'(5(1X,ES15.8))') shell_coor
 deallocate(shell_coor)

 ! skip 'Coordinates of each shell' in the old file
 k = 3*ncontr
 nline = k/5
 if(k - 5*nline > 0) nline = nline + 1
 do i = 1, nline, 1
  read(fid,'(A)') buf
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine replace_coor_in_fch

subroutine add_RI_kywd_into_molcas_inp(inpname, ricd)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: ricd ! T/F for RICD/Cholesky

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 if(ricd) then
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:4) == 'noCD') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while
  write(fid1,'(A)') 'RICD'
 else
  do while(.true.)
   read(fid,'(A)') buf
   write(fid1,'(A)') TRIM(buf)
   if(buf(1:7) == "&SEWARD") exit 
  end do ! for while
  write(fid1,'(A,/,A,/,A)') 'CHOLESKY','Threshold = 1d-14','Cutoff = 1d-16'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_RI_kywd_into_molcas_inp

! read memory and nprocshared from a given .gjf file
! mem is in MB
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
  k = INDEX(buf,'=')
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
  write(6,'(/,A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: incomplete&
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
 character(len=11), intent(out) :: method
 character(len=21), intent(out) :: basis

 j = INDEX(buf, '/')
 if(j == 0) then
  write(6,'(/,A)') "ERROR in subroutine read_method_and_basis_from_buf: no '/'&
                   & symbol found in"
  write(6,'(A)') "the '#' line. Failed to identify the method name. This utilit&
                 &y does not support"
  write(6,'(A)') "support syntax like 'M062X cc-pVDZ'. You should use the '/' &
                 &symbol like 'M062X/cc-pVDZ'."
  stop
 end if

 i = INDEX(buf(1:j-1), ' ', back=.true.)
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

 i = INDEX(buf(j+1:), ' ')
 basis = ' '
 basis = buf(j+1:j+i-1)
end subroutine read_method_and_basis_from_buf

subroutine read_title_card_from_gjf(gjfname, title)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 character(len=240), intent(out) :: title

 buf = ' '; title = ' '
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while

 read(fid,'(A)') title
 close(fid)
end subroutine read_title_card_from_gjf

