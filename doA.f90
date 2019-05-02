!--------------
 program doA
!--------------
! computes A matrix or inverse from a ped file
! tocompile:
!   gfortran doA.f90 -o doa
! usage:
!   cat pedfile | ./doa -nind nind [-inv] > afile
 implicit none
 integer   :: i, j, k, l, ii, is, id, nind=0, ios, n, p(3), par_stat
 character :: xc(20)*30, cmd*100
 real      :: x(20), w(3) = (/1., -.5, -.5/), &
                     res(4)  = (/2., 1.333333333, 1., 0./)
 real(SELECTED_REAL_KIND( 15, 307 )), allocatable  :: a(:,:)
 logical :: inverse=.false.

 call get_command (cmd)
 call nums2(cmd, n=n, x=x, xc=xc)
 do i=2, n
   select case (xc(i))
      case ('-n','-nind')
        nind=x(i+1)
      case ('-inv')
        inverse = .true.
   end select
 enddo
 if (nind==0) STOP('USAGE: ')
 allocate(A(nind,nind))
 A=0.

 !--> computes inverse of numerator relationship matrix
 if (inverse) then
    do i=1, nind
        read(*,*) p(1:3) 
        par_stat = 1
        if(p(2)==0) par_stat = par_stat+1
        if(p(3)==0) par_stat = par_stat+1
        do k=1,3
           do l=1,3
              if (p(k) /= 0 .and. p(l) /=0) then
                 A(p(k),p(l)) = A(p(k),p(l)) + w(k) * w(l) * res(par_stat)
              endif
           enddo
        enddo
    enddo

 !--> computes numerator relationship matrix
 else
     forall(i=1:nind) A(i,i)=1.
     do i=1, nind
        read(*,*) ii,is,id
        if(is>0 .and. id>0) a(ii,ii) = a(ii,ii) + 0.5*a(is,id)
        do j=1, i-1
           if(is>0) a(j,ii) = a(j,ii) + 0.5*a(j,is)
           if(id>0) a(j,ii) = a(j,ii) + 0.5*a(j,id)
           a(i,j) = a(j,i)
        enddo
     enddo
 endif

 do i=1, nind
    print*, a(i,:)
 enddo

CONTAINS

!----------------
 subroutine nums2 (a, n, x, xc)
!----------------
! separates array a into items delimited by blanks. character elements are
! put into optional character vector xc, decoded numeric values
! into optional real vector x, and n contains the number of items. The
! dimension of x and xc can be lower than n.
! format changed slightly to read real numbers in scientific notation (mpe)
! 2/23/2000 IMisztal

 character (*)          :: a
 character (*),optional :: xc(:)
 real,optional          :: x(:)
 integer :: n, curr, first, last, lena, stat, i

 curr=1
 lena=len(a)
 n=0

 do
   !--> search for first nonspace
   first=0
   do i=curr,lena
     if (a(i:i) /= ' ') then
        first=i
        exit
     endif
   enddo
   if (first == 0) exit

   !--> search for first space
   curr=first+1
   last=0
   do i=curr,lena
      if (a(i:i) == ' ') then
        last=i
        exit
      endif
   enddo

   if (last == 0) last=lena

   n=n+1
   if (present(xc)) then
      if (size(xc) >= n) then
         xc(n)=a(first:last)
      else
         print*, 'NUMS2: Increase size of XC'
      endif
   endif

   if (present(x)) then
      if (size(x) >= n) then
        read(a(first:last),*,iostat=stat) x(n)    !NEW (mpe)
        if (stat .ne. 0) x(n)=0
      else
         print*, 'NUMS2: Increase size of X'
      endif
   endif

   curr=last+1
 enddo
!--------------
 end subroutine
!--------------

end program


