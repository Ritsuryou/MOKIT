! written by jxzou at 20210831: do rigid/relaxed scan

! currently only rigid scan is supported
subroutine do_pes_scan()
 use mol, only: ndb, nacto, scan_itype
 use mr_keyword, only: rigid_scan, relaxed_scan, gjfname, scan_nstep, scan_val,&
  casnofch
 use util_wrapper, only: unfchk
 implicit none
 integer :: i, k, fid
 real(kind=8), external :: calc_an_int_coor
 character(len=24) :: data_string = ' '
 character(len=33), parameter :: error_warn='ERROR in subroutine do_pes_scan: '
 character(len=240) :: filelist, new_gjf, new_chk, new_fch, old_fch
 logical :: alive

 if(.not. (rigid_scan .or. relaxed_scan)) return
 write(6,'(//,A)') 'Enter subroutine do_pes_scan...'

 call read_scan_var_from_gjf()
 if(rigid_scan) then
  write(6,'(A)') 'Rigid scan values:'
 else
  write(6,'(A)') 'Relaxed scan values:'
  write(6,'(/,A)') error_warn//'relaxed scan is unsupported currently.'
  stop
 end if

 select case(scan_itype)
 !case(1,2,3)
 ! write(6,'(10F6.3)') (scan_val(i),i=1,scan_nstep)
 case(4)
  write(6,'(A)') 'using provided geometries from filelist.'
  call find_specified_suffix(gjfname, '.gjf', i)
  filelist = gjfname(1:i-1)//'.filelist'
  call find_nfile_in_filelist(filelist, scan_nstep)
 case default
  write(6,'(/,A)') error_warn//'currently only scan_itype=4 is supported.'
  write(6,'(A,I0)') 'But got scan_itype=', scan_itype
  stop
 end select

 write(6,'(2(A,I0))') 'MO projection using nmo = ndb + nacto. ndb=', ndb, &
                      ', nacto=', nacto
 open(newunit=fid,file=TRIM(filelist),status='old',position='rewind')
 ! Currently we use the Gaussian program to do MO projection, I do not know
 ! what happens to make_orb_resemble() and whether there exits any bug in that
 ! function.
 old_fch = casnofch

 do i = 1, scan_nstep, 1
  read(fid,'(A)') new_gjf
  inquire(file=TRIM(new_gjf),exist=alive)
  if(.not. alive) then
   write(6,'(/,A)') error_warn//TRIM(new_gjf)//' not found!'
   close(fid)
   stop
  end if
  call find_specified_suffix(new_gjf, '.gjf', k)
  new_chk = new_gjf(1:k-1)//'.chk'
  new_fch = new_gjf(1:k-1)//'.fch'
  call unfchk(old_fch, new_chk)
  call gen_fch_from_gjf(new_gjf, new_fch, .true., .true.)
  call add_automr_kywd2gjf(new_gjf)
  call submit_automr_job(new_gjf)
  old_fch = new_fch
 end do ! for i

 close(fid)
 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_pes_scan at '//TRIM(data_string)
end subroutine do_pes_scan

! read scan variables/coordinates from gjf
subroutine read_scan_var_from_gjf()
 use mol, only: scan_itype, scan_atoms
 use mr_keyword, only: gjfname, hf_fch, skiphf, scan_nstep, scan_val
 implicit none
 integer :: i, j, k, m, nblank0, nblank, fid
 integer, external :: detect_ncol_in_buf
 real(kind=8) :: rtmp0, stepsize
 real(kind=8), external :: calc_an_int_coor
 real(kind=8), allocatable :: coor(:,:), rtmp(:,:)
 character(len=240) :: buf
 character(len=44), parameter :: error_warn = 'ERROR in subroutine read_scan_va&
                                              &r_from_gjf: '
 scan_atoms = 0; nblank = 0; buf = ' '
 if(skiphf) then
  nblank0 = 2
 else
  nblank0 = 3
 end if

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == nblank0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') error_warn//'wrong format of scan coordinate.'
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 close(fid)

 call upper(buf)
 if(TRIM(buf) == 'FILELIST') then
  scan_itype = 4 ! use geometries from filelist
  return
 end if

 select case(buf(1:1))
 case('B')
  i = 5
 case('A')
  i = 6
 case('D')
  i = 7
 case default
  write(6,'(A)') error_warn//'invalid scan variable: '//buf(1:1)
  stop
 end select

 j = detect_ncol_in_buf(buf)
 if(j < i) then
  write(6,'(A)') error_warn//'invalid rigid scan syntax.'
  write(6,'(A)') 'buf='//TRIM(buf)
  stop
 end if

 select case(buf(1:1))
 case('B') ! bond
  scan_itype = 1
  read(buf(3:),*) scan_atoms(1:2), scan_nstep
 case('A') ! angle
  scan_itype = 2
  read(buf(3:),*) scan_atoms(1:3), scan_nstep
 case('D') ! dihedral
  scan_itype = 3
  read(buf(3:),*) scan_atoms(1:4), scan_nstep
 case default
  write(6,'(A)') error_warn//'invalid scan variable "'//buf(1:1)//'"'
  stop
 end select

 if(scan_nstep < 1) then
  write(6,'(A,I0)') error_warn//'invalid scan_nstep=', scan_nstep
  stop
 end if

 allocate(scan_val(scan_nstep), source=0d0)
 j = LEN_TRIM(buf)

 if(buf(j:j) == '}') then ! given a set of values {}
  i = INDEX(buf, '{')
  read(buf(i+1:j-1),fmt=*,iostat=m) (scan_val(k),k=1,scan_nstep)
  if(m /= 0) then
   write(6,'(/,A)') error_warn//'wrong scan syntax.'
   write(6,'(A)') 'buf='//TRIM(buf)
   stop
  end if

 else ! given the step length/interval
  ! calculate the current value of the target scan variable
  call read_natom_from_fch(hf_fch, k)
  allocate(coor(3,k))
  call read_coor_from_fch(hf_fch, k, coor)
  k = scan_itype + 1
  allocate(rtmp(3,k))
  forall(i = 1:k) rtmp(:,i) = coor(:,scan_atoms(i))
  deallocate(coor)
  rtmp0 = calc_an_int_coor(k, rtmp(:,1:k))
  deallocate(rtmp)
  i = INDEX(buf(1:j), ' ', back=.true.)
  read(buf(i+1:j),*) stepsize
  forall(i = 1:scan_nstep) scan_val(i) = rtmp0 + DBLE(i)*stepsize
 end if

 call check_scan_val(scan_itype, scan_nstep, scan_val)
end subroutine read_scan_var_from_gjf

! check whether the array scan_val is valid/reasonable
subroutine check_scan_val(scan_itype, n, scan_val)
 implicit none
 integer :: i, j
 integer, intent(in) :: scan_itype, n
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: scan_val(n)
 character(len=36), parameter :: error_warn = 'ERROR in subroutine check_scan_v&
                                              &ar: '
 do i = 1, n-1, 1
  rtmp = scan_val(i) - scan_val(i+1)

  select case(scan_itype)
  case(1) ! scan a bond
   if(rtmp < 0d0) then
    write(6,'(/,A)') error_warn//'scan values must be in descending order.'
    write(6,'(A,10F6.3)') 'scan_val=',(scan_val(j),j=1,n)
    stop
   end if
   if(rtmp<1d-3 .or. rtmp>5d0) then
    write(6,'(/,A)') error_warn//'intervals of the scanned bond are'
    write(6,'(A)') 'too large or too small. It should be >=0.001 and <=5.0 Angs&
                   &trom.'
    stop
   end if
  case(2,3) ! scan an angle or a dihedral
   rtmp = DABS(rtmp)
   if(rtmp<1d0 .or. rtmp>90d0) then
    write(6,'(/,A)') error_warn//'intervals of the scanned angle/dihedral'
    write(6,'(A)') 'are too large or too small. It should be >=1.0 and <=90.0 d&
                   &egree.'
    stop
   end if
  case default
   write(6,'(/,A,I0)') error_warn//'invalid scan_itype=',scan_itype
   stop
  end select
 end do ! for i
end subroutine check_scan_val

! add automr keywords into a specified .gjf
subroutine add_automr_kywd2gjf(gjfname)
 use mr_keyword, only: iroot, nstate, nmr, icss, polar, hardwfn, crazywfn, &
  casscf_prog, polar_prog
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, gjfname1, fchname
 character(len=240), intent(in) :: gjfname

 call find_specified_suffix(gjfname, '.gjf', i)
 gjfname1 = gjfname(1:i-1)//'.t'
 fchname = gjfname(1:i-1)//'.fch'

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(gjfname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while

 write(fid1,'(A)',advance='no') 'mokit{ist=5,readno='''//TRIM(fchname)//''''
 if(TRIM(casscf_prog) /= 'pyscf') then
  write(fid1,'(A)',advance='no') ',CASSCF_prog='//TRIM(casscf_prog)
 end if

 if(iroot > 0) write(fid1,'(A)',advance='no') ',Root=', iroot
 if(nstate > 0) write(fid1,'(A)',advance='no') ',Nstates=', nstate

 if(polar) then
  write(fid1,'(A)',advance='no') ',Polar,Polar_prog='//TRIM(polar_prog)
 end if

 if(nmr) write(fid1,'(A)',advance='no') ',NMR'
 if(icss) write(fid1,'(A)',advance='no') ',ICSS'
 if(hardwfn) write(fid1,'(A)',advance='no') ',HardWFN'
 if(crazywfn) write(fid1,'(A)',advance='no') ',CrazyWFN'
 write(fid1,'(A)') '}'

 read(fid,'(A)') buf ! skip one line
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(gjfname1), TRIM(gjfname))
end subroutine add_automr_kywd2gjf

