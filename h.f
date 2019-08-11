
!***************************************************
!***************************************************
module global


! parallel threads
integer :: threads,threads1
parameter( threads = 18)

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
parameter(kT = 3.8 )

! Electric field
real(kind=8) :: Ex

! Magnitude of polarization
real(kind=8) :: Px( L+2*R , L+2*R ),Py( L+2*R , L+2*R ),Pxs( L+2*R , L+2*R ),Pys( L+2*R , L+2*R )
! Demon of energy
integer(kind=8),parameter :: df=5
real(kind=8) :: Demon( L+2*R , L+2*R ,df )
real(kind=8) :: DemonProcess( L+2*R , L+2*R , Estep+1)

real(kind=8) :: Demon_temp( L+2*R , L+2*R ),DemonProcess_temp( L+2*R , L+2*R , Estep+1)

real(kind=8) :: ave_demon , ave_demon_process(Estep+1)

real(kind=8) :: DemonSum( Judgestep * t / threads , Estep+1  )
! Coordinates of the system
real :: xcoord( L+2*R , L+2*R ), ycoord( L+2*R , L+2*R )

! common loop expression
integer :: I,J,K,loop_m,id

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

integer :: n,n1,n2,n3,n4
real(kind=8) :: E(Estep+1),ii
real(kind=8) :: PxD2(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PyD2(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: temp(4)
real(kind=8),ALLOCATABLE :: qu1(:,:),qu2(:,:),qu3(:,:),qu4(:,:)
real(8) :: x
real(8), dimension(1000,2) :: curve1
real(8), dimension(1000,2) :: curve2
real(8), dimension(1000,2) :: xydata1
real(8), dimension(1000,2) :: xydata2
integer :: plot_type    
integer :: ret,start
character :: xlabel,ylabel,title1,title2 
CHARACTER(LEN=80) :: FMT,FMT2
real(kind=8),ALLOCATABLE :: Px_his(:,:,:),Py_his(:,:,:),Px_ave_his(:,:,:),Py_ave_his(:,:,:)
real(kind=8),ALLOCATABLE :: px_site_his(:),py_site_his(:), p_site(:)
! energy
real(kind=8) :: HeProcess( Estep+1 ),HfldProcess( Estep+1 ),HdipProcess( Estep+1 )
real(kind=8) :: HgrProcess( Estep+1 ) ,  HtotalProcess( Estep+1 ) ,deltaT(Estep+1)
real(kind=8) :: ht( Estep+1 ) 




!***************************************************
!***************************************************
! READ NECESSARY DATA
! Read file from E1.bin
open( unit = 10, file ='E1.bin' ,status = 'old', &
recl = (Estep+1)*8, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)E
close(unit=10)


open( unit = 12, file ='HeProcess1.bin', status = 'old' , access = 'direct' & 
, recl =  (Estep+1) *8 , form='unformatted')
read(unit=12,rec=1)HeProcess
close(unit=12)
HeProcess = HeProcess/L/L


! Read file for landau
open( unit = 12, file ='HfldProcess1.bin' ,status = 'old' , access = 'direct' & 
, recl =  (Estep+1) *8 , form='unformatted')
read(unit=12,rec=1)HfldProcess
close(unit=12)
HfldProcess = HfldProcess/L/L



open( unit = 12, file ='HdipProcess1.bin', status = 'old' , access = 'direct' & 
, recl =  (Estep+1) *8 , form='unformatted')
read(unit=12,rec=1)HdipProcess
close(unit=12)
HdipProcess = HdipProcess/L/L


open( unit = 12, file ='HgrProcess1.bin', status = 'old' , access = 'direct' & 
, recl =  (Estep+1) *8 , form='unformatted')
read(unit=12,rec=1)HgrProcess
close(unit=12)
HgrProcess = HgrProcess/L/L



open( unit = 12, file ='HtotalProcess1.bin', status = 'old' , access = 'direct' & 
, recl =  (Estep+1) *8 , form='unformatted')
read(unit=12,rec=1)HtotalProcess
close(unit=12)
HtotalProcess = HtotalProcess/L/L



open( unit = 12, file ='DemonSum1.bin' ,status = 'old' , access = 'direct' & 
, recl =  (Estep+1) * Judgestep * t / threads  *8 , form='unformatted')
read(unit=12,rec=1)DemonSum
close(unit=12)


open( unit = 12, file ='DemonProcess1.bin' ,status = 'old' , access = 'direct'  & 
, recl =  (L+2*R) * (L+2*R) * (Estep+1) *8 , form='unformatted' )
read(unit=12,rec=1)DemonProcess
close(unit=12)

deltaT = SUM ( SUM ( DemonProcess( (R+1):(L+R) , (R+1):(L+R) , :  ), DIM=1) ,DIM=1 )

deltaT = deltaT /L/L

! write(*,*)shape( sum(Demonsum,dim=1) )
! ! deltaT = sum(Demonsum,dim=1) / (Judgestep * t / threads) * df
write(*,*)shape(delta),shape(HtotalProcess),shape(ht)



Ht =  deltaT + HtotalProcess


! format for energy
! deltaT = deltaT / df




FMT2 = "(F6.1,F6.2,F6.2,F6.2,F6.2,F6.2,F6.2,F6.2)"
! !write file to energy.txt
open( unit=11,file='energy.txt',status = 'replace')
do i=1,Estep+1
write(11,FMT2)float(i),HeProcess(i),HfldProcess(i),HdipProcess(i),HgrProcess(i),HtotalProcess(i),deltaT(i),ht(i)
enddo
close(unit=11)

! ! END OF DOMAIN CONFIGURATION FIGURE
!***************************************************
!***************************************************


!create gnuplot command file
OPEN(10,ACCESS='SEQUENTIAL',FILE='gp.gnu')
! write(10,*) 'set terminal png size 1800,1800 enhanced font "Helvetica,20" crop'
write(10,*) ' set terminal postscript eps enhanced color solid 20'
!***************************************************


write(10,*) 'set output "energy.eps"'

!***************************************************
! plot hysteresis
! write(10,*) 'set xtics 0,100,900'
! write(10,*) 'set ytics -0.4,0.2,0.4'
! write(10,*) 'set xrange [0:900]'
! ! write(10,*) 'set yrange [-0.4:0.4]'
write(10,*) 'set key font ",10"'
write(10,*) 'set xlabel "steps','"offset 0,0.5'
write(10,*) 'set ylabel "Energy" offset 1.5,0'
write(10,*) 'plot &
"energy.txt" using 1:2 title "H_e",&
            "energy.txt" using 1:3 title "H_D", &
            "energy.txt" using 1:4 title "H_{dip}", &
            "energy.txt" using 1:5 title "H_{gr}", &
            "energy.txt" using 1:6 title "H_{total}" , &       
            "energy.txt" using 1:7 title "Demon", &
            "energy.txt" using 1:8 title "sum"'
CLOSE(10,STATUS='KEEP')

! ! END OF DOMAIN CONFIGURATION FIGURE
!***************************************************
!***************************************************
ret=SYSTEM('gnuplot gp.gnu')

!***************************************************
!***************************************************
! History of Ps history


End 

