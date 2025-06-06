! written by jxzou at 20190411: read and adjust the orders of d,f,g, etc functions
!  in Gaussian .chk file, to the orders in PySCF

! Step1:
!  transform the .chk file to _chk.txt file using Gaussian utility chkchk,
!  and read the MOs from _chk.txt
! Step2:
!  transform the .chk file to _tmp.fchk file using Gaussian utility formchk,
!  and read shell types from _tmp.fchk
!  (The _chk.txt and _tmp.fchk file will be deleted after reading)
! Step3:
!  adjust its d,f,g, etc. functions order of Gaussian to that of PySCF
subroutine chk2py(chkname, nbf, nif, ab, coeff)
 implicit none
 integer :: i, k, length, chkid
 integer :: ncoeff, nbf, nif
!f2py intent(in) :: nbf, nif
 integer :: system
 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell_to_atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)

 real(kind=8) :: coeff(nbf,nif)
!f2py depend(nbf,nif) :: coeff
!f2py intent(out) :: coeff

 character(len=1) :: ab
!f2py intent(in) :: ab
 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'
 character(len=240) :: chkname, txtname, buffer
!f2py intent(in) :: chkname
 logical :: alive

 inquire(file=TRIM(chkname),exist=alive)
 if(.not. alive) then
  write(6,'(A)') 'ERROR in subroutine chk2py: file does not exist!'
  write(6,'(A)') 'Filename='//TRIM(chkname)
  stop
 end if

 i = INDEX(chkname,'.chk')
 txtname = chkname(1:i-1)//'_chk.txt'

 ! Step1: transform .chk to _chk.txt using Gaussian utility chkchk
 i = SYSTEM('chkchk -p '//TRIM(chkname)//' > '//TRIM(txtname))
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine chk2py: fail to transform the .chk file t&
                 &o *_chk.txt file'
  write(6,'(A)') 'using Gaussian utility chkchk. File: '//TRIM(chkname)
  write(6,'(A)') "You can use 'which chkchk' to check if this command exists."
  stop
 end if

 ! read Alpha MO or Beta MO from _chk.txt
 buffer = ' '
 ncoeff = 0
 key = key1
 if(ab/='a' .and. ab/='A') then
  key = key2//' '
 end if

 open(newunit=chkid,file=TRIM(txtname),status='old',position='rewind')
 do while(.true.)
  read(chkid,'(A)') buffer
  if(buffer(7:14) == key) exit
 end do
 BACKSPACE(chkid)

 coeff = 0d0
 do i = 1, nif, 1
  read(chkid,'(A)') buffer ! skip one line
  read(chkid,'(4D20.12)') (coeff(k,i), k=1,nbf)
 end do
 ! MOs in _chk.txt read done
 close(chkid,status='delete')


 ! Step2: transform .chk to _tmp.fchk using Gaussian utility formchk
 i = INDEX(chkname,'.chk')
 txtname = chkname(1:i-1)//'_tmp.fchk'
 i = SYSTEM('formchk '//TRIM(chkname)//' '//TRIM(txtname)//" >junk_tmp 2>&1")
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine chk2py: fail to transform the .chk file&
                   & to _tmp.fchk'
  write(6,'(A)') 'file, using Gaussian utility formchk.'
  write(6,'(A)') TRIM(chkname)
  write(6,'(A)') "You can use 'which formchk' to check if this command exists."
  stop
 end if

 ! Normally, the content in file junk_tmp is:
 ! ' Read checkpoint file xxx.chk'
 ! ' Write formatted file xxx.fchk'
 ! delete the file junk_tmp
 open(newunit=chkid,file='junk_tmp',status='old')
 close(chkid,status='delete')

 open(newunit=chkid,file=TRIM(txtname),status='old',position='rewind')
 do while(.true.)
  read(chkid,'(A)') buffer
  if(buffer(1:17) == 'Number of basis f') exit
 end do
 BACKSPACE(chkid)

 read(chkid,'(A49,2X,I10)') buffer, i
 if(i /= nbf) then
  write(6,'(A)') 'ERROR in subroutine chk2py: inconsistent number of basis&
                   & functions in .chk file and in PySCF script.'
  write(6,'(A)') TRIM(chkname)
  close(chkid)
  stop
 end if

 read(chkid,'(A49,2X,I10)') buffer, i
 if(i /= nif) then
  write(6,'(A)') 'ERROR in subroutine chk2py: inconsistent number of MOs in&
                   & .chk file and in PySCF script.'
  write(6,'(A)') TRIM(chkname)
  close(chkid)
  stop
 end if

 ! find and read Shell types
 do while(.true.)
  read(chkid,'(A)',iostat=i) buffer
  if(buffer(1:11) == 'Shell types') exit
 end do
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine chk2py: missing the 'Shell types'&
                   & section in .fchk file!"
  write(6,'(A)') TRIM(chkname)
  close(chkid)
  stop
 end if

 BACKSPACE(chkid)
 read(chkid,'(A49,2X,I10)') buffer, k
 allocate(shell_type(2*k))
 shell_type = 0
 read(chkid,'(6(6X,I6))') (shell_type(i),i=1,k)
 ! read Shell types done

 ! find and read Shell to atom map
 do while(.true.)
  read(chkid,'(A)',iostat=i) buffer
  if(buffer(1:13) == 'Shell to atom') exit
 end do
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine chk2py: missing the 'Shell to atom map'&
                   & section in .fchk file!"
  write(6,'(A)') TRIM(chkname)
  close(chkid)
  stop
 end if
 allocate(shell_to_atom_map(2*k))
 shell_to_atom_map = 0
 read(chkid,'(6(6X,I6))') (shell_to_atom_map(i),i=1,k)
 ! read Shell to atom map done

 ! all information in .fchk file read done
 close(chkid,status='delete')


! Step3:
! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell_to_atom_map, length)

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell_to_atom_map, nbf, nif, coeff)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n6dmark = 0
 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 n5dmark = 0
 n7fmark = 0
 n9gmark = 0
 n11hmark = 0
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))
 d_mark = 0
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do i = 1,k,1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  end select
 end do
 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions
 do i = 1,n5dmark,1
  call chk2py_permute_5d(nif,coeff(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1,n6dmark,1
  call chk2py_permute_6d(nif,coeff(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1,n7fmark,1
  call chk2py_permute_7f(nif,coeff(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1,n10fmark,1
  call chk2py_permute_10f(nif,coeff(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1,n9gmark,1
  call chk2py_permute_9g(nif,coeff(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1,n15gmark,1
  call chk2py_permute_15g(nif,coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1,n11hmark,1
  call chk2py_permute_11h(nif,coeff(h_mark(i):h_mark(i)+10,:))
 end do
 do i = 1,n21hmark,1
  call chk2py_permute_21h(nif,coeff(h_mark(i):h_mark(i)+20,:))
 end do
! adjustment finished

 deallocate(d_mark, f_mark, g_mark, h_mark)
end subroutine chk2py

subroutine chk2py_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in PySCF
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 allocate(coeff2(5,nif), source=0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine chk2py_permute_5d

subroutine chk2py_permute_6d(nif,coeff)
 use Sdiag_parameter, only: Sdiag_d
 implicit none
 integer :: i
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in PySCF
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 allocate(coeff2(6,nif), source=coeff)
 forall(i = 1:6) coeff(i,:) = coeff2(order(i),:)/Sdiag_d(i)
 deallocate(coeff2)
end subroutine chk2py_permute_6d

subroutine chk2py_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in PySCF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 allocate(coeff2(7,nif), source=0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine chk2py_permute_7f

subroutine chk2py_permute_10f(nif,coeff)
 use Sdiag_parameter, only: Sdiag_f
 implicit none
 integer :: i
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in PySCF
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 allocate(coeff2(10,nif), source=coeff)
 forall(i = 1:10) coeff(i,:) = coeff2(order(i),:)/Sdiag_f(i)
 deallocate(coeff2)
end subroutine chk2py_permute_10f

subroutine chk2py_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in PySCF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 allocate(coeff2(9,nif), source=0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine chk2py_permute_9g

subroutine chk2py_permute_15g(nif,coeff)
 use Sdiag_parameter, only: Sdiag_g
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 allocate(coeff2(15,nif), source=coeff)
 forall(i = 1:15) coeff(i,:) = coeff2(16-i,:)/Sdiag_g(i)
 deallocate(coeff2)
end subroutine chk2py_permute_15g

subroutine chk2py_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 allocate(coeff2(11,nif), source=0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine chk2py_permute_11h

subroutine chk2py_permute_21h(nif,coeff)
 use Sdiag_parameter, only: Sdiag_h
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in PySCF
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 allocate(coeff2(21,nif), source=coeff)
 forall(i = 1:21) coeff(i,:) = coeff2(22-i,:)/Sdiag_h(i)
 deallocate(coeff2)
end subroutine chk2py_permute_21h

