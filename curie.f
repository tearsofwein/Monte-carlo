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
parameter( threads = 12)

! natural parameters
real(kind=8),parameter :: pi = 4.0 * atan ( 1.0 )

! parameters  for lattice
integer :: L,R,num_loc
parameter( L = 63 , R = 8 , num_loc = 0 )
! Size of the cell
real(kind=8) :: det_x,det_y
parameter( det_x = 1.0 , det_y = 1.0 )

! Parameters for Landau free energy
real(kind=8) :: a,b,c
parameter( a = -6.7284  ,b = 15.007 )

! Parameters for gradient energy
real(kind=8) :: g1,g2,g22
parameter( g1 = 1.2267 , g2 =  1.2267 , g22 = 1.2267 )

! Prefactors
! SJ: prefactor for gradient energy, SE:electric field, SA:dipole-dipole interaction
real(kind=8) :: SJ,SE,SA
parameter( SJ = 2.0 , SE = 1.0 , SA = 7.35 )


! Numbers of loop steps
! Judgestep * t : steps for reaching the equilibrium
! Estep: electric steps
! 't' is a very important parameter, which should be at least 3 times greater than
! threads. Because t/threads is the times of one nodes being visited.
integer,parameter :: Judgestep = L*L , Estep = 1800,t = 36


! frequency
real(kind=8),parameter :: f = 0.001

! Parameters for Temperature
real(kind=8) :: kT,kT_temp
parameter(kT = 8.0 )

! Electric field
real(kind=8) :: Ex

! Magnitude of polarization
real(kind=8) :: Px( L+2*R , L+2*R ),Py( L+2*R , L+2*R ),Pxs( L+2*R , L+2*R ),Pys( L+2*R , L+2*R )
! Demon of energy
integer(kind=8),parameter :: df=1
real(kind=8) :: Demon( L+2*R , L+2*R ,df )
real(kind=8) :: DemonProcess( L+2*R , L+2*R , Estep+1)

real(kind=8) :: Demon_temp( L+2*R , L+2*R ),DemonProcess_temp( L+2*R , L+2*R , Estep+1)

real(kind=8) :: ave_demon , ave_demon_process(Estep+1)

real(kind=8) :: DemonSum( Judgestep * t / threads , Estep+1  )
! Coordinates of the system
real :: xcoord( L+2*R , L+2*R ), ycoord( L+2*R , L+2*R )

! common loop expression
integer :: I,J,K,loop_m

! invariable for range setting
integer :: rangetotal,RR,RRR,RL,RLL

! coordinates for EnergyofSelectnoeds
real(kind=8) ::  dists(2*R+1,2*R+1),dist2s(2*R+1,2*R+1),  dxcoords(2*R+1,2*R+1), dycoords(2*R+1,2*R+1)

! Px and Py without boundary
real(kind=8) :: PxTotal( 1:L , 1:L ) , PyTotal( 1:L , 1:L )

! local electric field influence parameter
real(kind=8) :: Ex_Rnd( L+2*R , L+2*R ), Ey_Rnd( L+2*R , L+2*R )

! the energy profile at each site
real(kind=8) :: Energy_profile(2*R+L , 2*R+L )

common Ex_Rnd,Ey_Rnd,Ex,Px,Py,Energy_profile,kT_temp,dists,dxcoords,dycoords,PxTotal,Pytotal
common Demon,rangetotal,RR,RRR,RL,RLL

end module
!***************************************************
!***************************************************


Program read_test
use global
USE s_wrtmat
integer :: n,n1,n2,n3,n4
real(kind=8) :: E(Estep+1),ii,Px_ave(L,L),Py_ave(L,L)
real(kind=8) :: PxD2(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PyD2(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PxD2_1(L , L , Estep+1 )
real(kind=8) :: PyD2_1(L , L , Estep+1 )
real(kind=8) :: temp(4)
real(kind=8) :: P_ave_site

real(kind=8) :: Temperature
real(kind=8),ALLOCATABLE :: qu1(:,:),qu2(:,:),qu3(:,:),qu4(:,:)
real(8) :: x
real(8), dimension(1000,2) :: curve1
real(8), dimension(1000,2) :: curve2
real(8), dimension(1000,2) :: xydata1
real(8), dimension(1000,2) :: xydata2
integer :: plot_type    
integer :: ret,start
character :: xlabel,ylabel,title1,title2
character(LEN=80) :: plotinput
CHARACTER(LEN=80) :: FMT
real(kind=8),ALLOCATABLE :: Px_his(:,:,:),Py_his(:,:,:),Px_ave_his(:,:,:),Py_ave_his(:,:,:)
real(kind=8),ALLOCATABLE :: px_site_his(:),py_site_his(:), p_site(:),p_site_his(:)
real(kind=8),allocatable :: Px_2(:), Py_2(:)




!***************************************************
!***************************************************

! Read file from PxProcess1.bin
open( unit = 10, file ='PxProcess1.bin' ,status = 'old', &
recl =  (L+2*R) * (L+2*R) * (Estep+1) *8, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)PxD2
close(unit=10)


! Read file from PyProcess1.bin
open( unit = 10, file ='PyProcess1.bin' ,status = 'old', &
recl =  (L+2*R) * (L+2*R) * (Estep+1) *8, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)PyD2
close(unit=10)




! Read file from Kt1.bin
open( unit = 10, file ='KT1.bin' ,status = 'old', &
recl = (Estep+1)*8, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)Temperature
close(unit=10)


! ! Calculate the average polarization
!k=start-1
N2 = 100
N3 = 200
allocate(Px_his(L,L,N3-N2+1))
allocate(Py_his(L,L,N3-N2+1))
allocate(Px_ave_his(L,L,N3-N2+1))
allocate(Py_ave_his(L,L,N3-N2+1))
allocate(px_site_his(N3-N2+1))
allocate(py_site_his(N3-N2+1))
allocate(p_site(N3-N2+1))
allocate(p_site_his(N3-N2+1))
allocate(Px_2(N3-N2+1))
allocate(Py_2(N3-N2+1))
Px_his = 0.0
Py_his = 0.0
start = 2
Px_his(:,:,1) =  PxD2( (R+1):(R+L), (R+1):(R+L) ,N2)
Py_his(:,:,1) =  PyD2( (R+1):(R+L), (R+1):(R+L) ,N2)

PxD2_1 = PxD2 ( (R+1) : (L+R) , (R+1) : (L+R) , 1:N ) 
PyD2_1=  PyD2 ( (R+1) : (L+R) , (R+1) : (L+R) , 1:N ) 
! %%%%%%%%%%%%%%%%%




! %%%%%%%%%%%%%%%%%
! % concept 2013 12 04
do i=1,L
do j=1,L
Px_2 = PxD2_1(i,j,N2:N3)
Px_ave(i,j) = sum( Px_2) / (N3-N2+1)
enddo
enddo	




do i=1,L
do j=1,L
Py_2 = PyD2_1(i,j,N2:N3)
Py_ave(i,j) = sum( Py_2 ) / (N3-N2+1)
enddo
enddo

P_ave_site = sum(sqrt(Px_ave**2+Py_ave**2 ) )  / (L*L)
write(*,*)P_ave_site


!!write file to  pshist.txt
FMT   = "(F10.5, F10.5)"
open( unit = 10, file ='curie.txt',status = 'replace')
do i=1,1
write(10,FMT) Temperature , P_ave_site
enddo
close(unit=10)




! ! END OF DOMAIN CONFIGURATION FIGURE
!***************************************************
!***************************************************


End 
