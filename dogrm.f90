! gfortran dogrm.f90 -lblas -O4 -o dogrm
! -hap option removed, polyploids allowed
! cat genotype_file | ./dogrm -nind Nind [-ploidy ploidy] [-dom] [-maf maf] [-maxmiss maxmiss]'
! if -dom computes Vitezika's 2013 dominance variance 
! if -hap genotypes instead of alleles are assumed to be read, eg, 1 2 instead of 1 0 1 1
program dogrm
implicit none
integer, parameter  :: nsnp_block=10000

integer   :: i, j, n, ios, isnp, nind=0, nmiss, ploidy=2, hap=2
real(SELECTED_REAL_KIND( 15, 307 )), allocatable :: grm(:,:), g(:,:), ig(:)
real(SELECTED_REAL_KIND( 15, 307 ))  :: one=1.d0, freq, maf=1.d-6, grm_d
character :: cmd*200, xc(10)*20
real      :: x(size(xc)), maxmiss=0.
logical   :: dom = .false.

call get_command (cmd)
call nums2(cmd, n=n, x=x, xc=xc)
do i=2, n
   select case (xc(i))
      case ('-dom')
        dom = .true.
      !case ('-hap')
      !  hap = 1
      case ('-nind')
        nind=x(i+1)
      case ('-maf')
        maf=x(i+1)
      case ('-maxmiss')
        maxmiss = x(i+1)
      case ('-ploidy')
        ploidy = x(i+1)
   end select
enddo
hap = ploidy
if (nind==0) &
   STOP 'Usage: cat genotype.txt | ./dogrm -nind nind [-maf maf] [-ploidy ploidy] [-dom] [-hap] > matrix.grm '
if (dom .and. ploidy>2) STOP 'Dominance GRM defined only for diploids'
allocate(grm(nind,nind), g(nind,nsnp_block), ig(nind*hap))
grm=0
grm_d=0
isnp=0

!--> usual GRM
if (.not. dom) then
   do
      read(*,*,iostat=ios) ig(:)
      if(ios.ne.0) EXIT
      !--> n missing, skip if larger than allowed
      nmiss = count(ig>ploidy)
      if (real(nmiss)/nind > maxmiss) CYCLE
      freq = sum(pack(ig, ig<ploidy)) / ((nind-nmiss)*ploidy)
      if(min(freq,1.-freq) < maf) CYCLE
      !--> replace miss with mean
      where(ig>ploidy) ig=freq/ploidy
      isnp=isnp+1
      g(:,isnp) = ig
      !--> sum in ploidy chunks or in genotypes (hap=1)
      g(:,isnp) = (/(sum(ig(j:(j+hap-1))), j=1,nind*hap,hap)/)
      g(:,isnp) = g(:,isnp) - ploidy*freq
      !--> denominator
      grm_d = grm_d + ploidy*freq*(1.-freq)

      if (mod(isnp,nsnp_block)==0) then
         !--> computes XX' and adds to GRM
         call dgemm('n','t',nind,nind,isnp,one, g(:,1:isnp), nind, g(:,1:isnp),nind,one, GRM,nind)
         isnp=0
      endif
   enddo

!--> dominance variance GRM
else
   do
      read(*,*,iostat=ios) ig(:)
      if(ios.ne.0) EXIT
      nmiss = count(ig>ploidy)
      if (real(nmiss)/nind > 0) CYCLE
      freq = sum(ig(:)) / (nind*ploidy)
      if(min(freq,1.-freq) < maf) CYCLE
      isnp=isnp+1
      g(:,isnp) = ig
      !--> sum in ploidy chunks or in genotypes (hap=1)
      g(:,isnp) = (/(sum(ig(j:(j+hap-1))), j=1,nind*hap,hap)/)

      !--> freq is frequency of allele '1'
      where (g(:,isnp)==0) g(:,isnp) = -2*freq**2
      where (g(:,isnp)==1) g(:,isnp) = 2.*freq*(1.-freq)
      where (g(:,isnp)==2) g(:,isnp) = -2*(1.-freq)**2

      grm_d = grm_d + (ploidy*freq*(1.-freq))**2

      if (mod(isnp,nsnp_block)==0) then
         !--> computes XX' and adds to GRM
         call dgemm('n','t',nind,nind,isnp,one, g(:,1:isnp), nind, g(:,1:isnp),nind,one, GRM,nind)
         isnp=0
      endif
   enddo

endif

!--> remaining snps
if (mod(isnp,nsnp_block)/=0) then
   call dgemm('n','t',nind,nind,isnp,one, g(:,1:isnp), nind, g(:,1:isnp),nind,one, GRM,nind)
endif

grm = grm / grm_d

do i=1, nind
   write(*,*) grm(i,:)
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


