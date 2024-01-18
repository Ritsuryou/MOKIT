! written by jxzou at 20210207: originally copied from bas_gau2molcas.f90

! Note: Currently isotopes are not tested.
program main
 implicit none
 integer :: i
 character(len=240) :: fname
 ! fname: file contains only basis sets data in Gaussian format

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine bas_gau2molcas: wrong command line&
                      & argument!'
  write(6,'(A,/)') ' Example: bas_gau2molcas cc-pVTZ.gbs (generate CC-PVTZ)'
  stop
 end if

 fname = ' '
 call getarg(1, fname)
 call require_file_exist(fname)
 call bas_gau2molcas(fname)
end program main

! Transform the basis sets in Gaussian format to those in (Open)Molcas format
subroutine bas_gau2molcas(inpname)
 use pg, only: prim_gau, natom, nuc, elem
 use fch_content, only: elem2nuc
 implicit none
 integer :: i, nline, fid1, fid2, RENAME
 character(len=240), intent(in) :: inpname
 character(len=240) :: buf, outname, gmslike
 character(len=1) :: stype
 character(len=21) :: str1, str2

 ! initialization
 buf = ' '
 outname = ' '

 i = INDEX(inpname, '.gjf')
 if(i == 0) i = INDEX(inpname, '.com')
 if(i > 0) then
  write(6,'(/,A)') 'ERROR in subroutine bas_gau2molcas: .gjf/.com file not supp&
                   &orted currently.'
  stop
 end if

 i = INDEX(inpname, '.gbs', back=.true.)
 if(i == 0) then
  i = LEN_TRIM(inpname)
 else
  i = i - 1
 end if
 gmslike = inpname(1:i)//'.gmsbs'
 call bas_gau2gmslike(inpname, gmslike)
 ! transform to GAMESS-like basis set data format, so we can call subroutines
 ! read_prim_gau1 and read_prim_gau2

 outname = inpname(1:i)
 call upper(outname)
 if(outname == inpname) i = RENAME(TRIM(inpname), TRIM(inpname)//'.bak')

 natom = 1
 allocate(nuc(1), source=0) ! atomic number
 allocate(elem(1))
 elem = ' '

 open(newunit=fid1,file=TRIM(gmslike),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(outname),status='replace')
 write(fid2,'(A)') '* generated by utility bas_gau2molcas in MOKIT'

 call clear_prim_gau() ! initialization: clear all primitive gaussians

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  read(buf,*) elem(1)
  nuc(1) = elem2nuc(elem(1))

  ! deal with primitive gaussians
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
   call read_prim_gau(stype, nline, fid1)
  end do ! for while

  call gen_contracted_string(prim_gau(:)%nline,prim_gau(:)%ncol,str1,str2)
  write(fid2,'(A)') '/'//TRIM(elem(1))//'.'//TRIM(outname)//'..'//TRIM(str1)//&
                  & '.'//TRIM(str2)//'.'
  write(fid2,'(A,/,A)') 'Comment1', 'Comment2'

  ! print basis sets of this atom in Molcas format
  call prt_prim_gau_molcas2(fid2)
  call clear_prim_gau() ! clear/deallocate for next cycle
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 deallocate(nuc, elem)
end subroutine bas_gau2molcas

! print primitive gaussians
subroutine prt_prim_gau_molcas2(fid)
 use pg, only: prim_gau, nuc, highest
 implicit none
 integer :: i, j, k, nline, ncol
 integer, intent(in) :: fid

 call get_highest_am()
 write(fid,'(5X,I0,A1,3X,I1)') nuc(1), '.', highest

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)') '* '//prim_gau(i)%stype//'-type functions'
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(2(1X,I4))') nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(3X,ES16.9)') prim_gau(i)%coeff(j,1)
  end do ! for j
  do j = 1, nline, 1
   write(fid,'(10(ES16.9,2X))') (prim_gau(i)%coeff(j,k), k=2,ncol)
  end do ! for j
 end do ! for i

end subroutine prt_prim_gau_molcas2

! transform to GAMESS-like basis set data format, so we can call subroutines
! read_prim_gau1 and read_prim_gau2
subroutine bas_gau2gmslike(inpname, gmslike)
 implicit none
 integer :: i, j, k, nline, fid1, fid2
 real(kind=8) :: rtmp(3)
 character(len=2) :: stype
 character(len=3) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname, gmslike

 elem = ' '
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(gmslike),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:1) /= '!') exit
 end do ! for while
 BACKSPACE(fid1)

 do while(.true.) ! while 1
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit

  read(buf,*) elem
  if(elem(1:1) == '-') elem = elem(2:3)//' '
  write(fid2,'(A)') TRIM(elem)

  do while(.true.) ! while 2
   read(fid1,'(A)') buf
   if(buf(1:4) == '****') exit

   read(buf,*) stype, nline
   if(stype == 'SP') then
    stype = 'L'
    k = 3
   else
    k = 2
   end if
   write(fid2,'(3X,A1,2X,I0)') TRIM(stype), nline
 
   rtmp = 0d0
   do j = 1, nline, 1
    read(fid1,*) rtmp(1:k)
    write(fid2,'(2X,I2,3(2X,ES16.9))') j, rtmp(1:3)
   end do ! for j
  end do ! for while 2

  write(fid2,'(/)',advance='no')
 end do ! for while 1

 close(fid1)
 close(fid2)
end subroutine bas_gau2gmslike

