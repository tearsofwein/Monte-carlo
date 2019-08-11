MODULE s_wrtmat
CONTAINS

! Subroutine of wrtmat
! #########################################################################
  SUBROUTINE wrtmat(a)

! -------------------------------------------------------------------------
! Purpose  :  Routine to display matrix a
!
! Author   :  Feng Zhou
!
! Created  :  09-Nov-2013
!
! Last mod.:  
!
! Changes  :  09-Nov-2013 ZF: Create prototype of wrtmat
!
! SR called:  
!
! Remarks  :  
!
! Reference:  
! -------------------------------------------------------------------------

! List of Parameters
! ---------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT(in) :: a      ! Input matrix a (irow x icol)

! Local Variables
! ---------------------------------
    INTEGER(KIND=4), PARAMETER :: LFNERR = 6
    INTEGER(KIND=4) :: i, j, irow, icol
    
    irow = SIZE(a,1)
    icol = SIZE(a,2)
    WRITE(LFNERR,'(a,i3,a,i3,a)') ' Matrix has ', irow, ' rows and ', icol, 'columns'
    DO i = 1, irow
       DO j = 1, icol
          WRITE(LFNERR,'(1x,f8.4,$)') a(i,j)
       END DO
       WRITE(LFNERR,*) ' '
    END DO
    

  END SUBROUTINE wrtmat
  
  END MODULE s_wrtmat




! Before programing input 'ulimit -s unlimited'
! Before programing input 'ulimit -c unlimited'
! The above two steps are required to release the limitation of the stack size
! E=0
! Ex_Rnd = 0.7
! Ey_Rnd = 0.7
! One option: Let Polarizition have value, whose range is randomly from 0 to 1
! , not guassian distribution. Another option is guassian distribution.


!***************************************************
!***************************************************
module global


! parallel threads
integer :: threads,threads1
parameter( threads = 12 )

! natural parameters
real(kind=8),parameter :: pi = 4.0 * atan ( 1.0 )

! parameters  for lattice
integer :: L,R,num_loc
parameter( L = 63 , R = 8 , num_loc = 0 )


! Numbers of loop steps
! Judgestep * t : steps for reaching the equilibrium
! Estep: electric steps
! 't' is a very important parameter, which should be at least 3 times greater than
! threads. Because t/threads is the times of one nodes being visited.
integer,parameter :: Judgestep = L*L , Estep = 1800,t = 36


real(kind=8) :: DemonSum( Judgestep * t / threads , Estep+1  )
! Coordinates of the system


common Demonsum

end module
!***************************************************
!***************************************************



Program read_test
use global
USE s_wrtmat

real(kind=8) :: DemonSum1( Judgestep * t / threads , Estep+1  )
integer(kind=8) :: N2
real(kind=8) :: deltaTmean
real(kind=8) :: Temperature
CHARACTER(LEN=80) :: FMT   ! definition for format
!***************************************************
!***************************************************

! Read file from DemonSum1.bin
open( unit = 10, file ='DemonSum1.bin' ,status = 'old' ,  &
recl =  (Estep+1) * Judgestep * t / threads  *8 ,access = 'direct' , &
form='unformatted')
read(unit=10,rec=1)DemonSum1
close(unit=10)



! Read file from Kt1.bin
open( unit = 10, file ='KT1.bin' ,status = 'old', &
recl = (Estep+1)*8, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)Temperature
close(unit=10)


! ! Calculate the average polarization
!k=start-1
write(*,*)size(DemonSum1,dim=2)
N2 = 1600
deltaTmean = sum( DemonSum1(1:(Judgestep * t / threads) , N2:(Estep+1) ) - Temperature ) &
/ (Estep+2-N2) / (Judgestep * t / threads)



write(*,*) deltaTmean

!!write file to  pshist.txt
FMT   = "(F10.5,F10.5)"
open( unit = 10, file ='temperature.txt',status = 'replace')
do i=1,1
write(10,FMT) Temperature,deltaTmean
enddo
close(unit=10)


! ! END OF DOMAIN CONFIGURATION FIGURE
!***************************************************
!***************************************************


End 
