

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
parameter( L = 63 , R = 8 , num_loc = 3493 )
! Size of the cell
real(kind=8) :: det_x,det_y
parameter( det_x = 1.0 , det_y = 1.0 )

! Parameters for Landau free energy
real(kind=8) :: a,b,a1,b1
parameter( a = -6.7284  ,b = 15.007 )
parameter( a1 = 2.9989  ,b1 = 0.75034 )
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
parameter(kT = 1.0 )

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
integer :: I,J,K,loop_m

! invariable for range setting
integer :: rangetotal,RR,RRR,RL,RLL

! coordinates for EnergyofSelectnoeds
real(kind=8) ::  dists(2*R+1,2*R+1),dist2s(2*R+1,2*R+1),  dxcoords(2*R+1,2*R+1), dycoords(2*R+1,2*R+1)

! Px and Py without boundary
real(kind=8) :: PxTotal( L , L ) , PyTotal( L , L ) , rndmap( L , L )

! local electric field influence parameter
real(kind=8) :: Ex_Rnd( L+2*R , L+2*R ), Ey_Rnd( L+2*R , L+2*R )

! the energy profile at each site
real(kind=8) :: Energy_profile(2*R+L , 2*R+L )

common Ex_Rnd,Ey_Rnd,Ex,Px,Py,Energy_profile,kT_temp,dists,dxcoords,dycoords,PxTotal,Pytotal
common Demon,rndmap,rangetotal,RR,RRR,RL,RLL

end module
!***************************************************
!***************************************************


Program read_test
use global
integer :: n,n1,n2,n3,n4,n5,ndomain,ii,ibegin=400,iend=400
real(kind=8) :: E(Estep+1)
real(kind=8) :: PxD2(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PyD2(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PxD3(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PyD3(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: temp(4)
real(kind=8),ALLOCATABLE :: qu1(:,:),qu2(:,:),qu3(:,:),qu4(:,:),qu5(:,:)
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






!***************************************************
!***************************************************
! READ NECESSARY DATA
! Read file from E1.bin
open( unit = 10, file ='E1.bin' ,status = 'old', &
recl = (Estep+1)*8, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)E
close(unit=10)
E =  E * 46.875 

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


! Read file from xaxis1.bin
open( unit = 10, file ='xaxis1.bin' ,status = 'old', &
recl = ( L+2*R ) * ( L+2*R ) * 4, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)xcoord
close(unit=10)


! Read file from yaxis1.bin
open( unit = 10, file ='yaxis1.bin' ,status = 'old', &
recl = ( L+2*R ) * ( L+2*R ) * 4, access="DIRECT" &
, form='unformatted' )
read(unit=10,rec=1)ycoord
close(unit=10)
! END OF READING NECESSARY DATA
!***************************************************
!***************************************************

! format for domain
FMT   = "(F10.5, F10.5 , F10.5 , F10.5)"
! format for hysteresis
FMT2 = "(F10.5,F10.5)"

! original value
PxD3 = 0.4 * PxD2 
PyD3 = 0.4 * PyD2

 
do ii=ibegin,iend
ndomain=ii

!***************************************************
!***************************************************
! ! DOMAIN CONFIGURATION FIGURE
n1 = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0
do n =ndomain,ndomain
     do i=R+1,(L+R+1)
     do j=R+1,(L+R+1)
     	
     if (PxD2(i,j,n)<0.0 .AND. PxD2(i,j,n)>-1.0) then     
     
     n1 = n1 +1 
     
     elseif (PxD2(i,j,n)>0.0 .AND. PxD2(i,j,n)<1.0) then
     
     n2 = n2 +1 

     elseif (PyD2(i,j,n)<0.0) then
     
     n3 = n3 +1 
     
     elseif (PyD2(i,j,n)>0.0) then  
     
     n4 = n4 +1      
     
     else
     
      n5 = n5 +1      
     
     endif
     enddo
     enddo
enddo


  allocate(qu1(n1,4))
  allocate(qu2(n2,4))
  allocate(qu3(n3,4))
  allocate(qu4(n4,4))
  allocate(qu5(n5,4))

! 
n1 = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0

do n =ndomain,ndomain
     do i=R+1,(L+R+1)
     do j=R+1,(L+R+1)
     if (PxD2(i,j,n)<0.0 .AND. PxD2(i,j,n)>-1.0) then     
      n1 = n1 +1 
      qu1(n1,:) = (/DBLE(xcoord(i,j)),DBLE(ycoord(i,j)),PxD2(i,j,n),DBLE(0.0)/)

     elseif (PxD2(i,j,n)>0.0 .and. PxD2(i,j,n)<1.0) then
      n2 = n2 +1 
      qu2(n2,:) = (/DBLE(xcoord(i,j)),DBLE(ycoord(i,j)),PxD2(i,j,n),DBLE(0.0)/)     

     elseif (PyD2(i,j,n)<0.0) then     
      n3 = n3 +1 
      qu3(n3,:) = (/DBLE(xcoord(i,j)),DBLE(ycoord(i,j)),DBLE(0.0),PyD2(i,j,n)/)

     elseif (PyD2(i,j,n)>0.0) then
      n4 = n4 +1 
      qu4(n4,:) = (/DBLE(xcoord(i,j)),DBLE(ycoord(i,j)),DBLE(0.0),PyD2(i,j,n)/)

     else
     n5 = n5 +1 
     qu5(n5,:) = (/DBLE(xcoord(i,j)),DBLE(ycoord(i,j)),PxD2(i,j,n),DBLE(0.0)/)

     endif
     enddo
     enddo
enddo


   
! !write file to  qui2.TXT
if (n1==0) then
n1 = 0
else
open( unit = 10, file ='qui1.txt',status = 'replace'  )
do i=1,n1
write(10,FMT)qu1(i,:)
enddo
close(unit=10)
endif
! 
! !write file to  qui2.TXT
if (n2==0) then
n2 = 0
else
open( unit = 10, file ='qui2.txt',status = 'replace'  )
do i=1,n2
write(10,FMT)qu2(i,:)
enddo
close(unit=10)
endif
! 
! !write file to  qui3.TXT
if (n3==0) then
n3 = 0
else
open( unit = 10, file ='qui3.txt',status = 'replace')
do i=1,n3
write(10,FMT)qu3(i,:)
enddo
close(unit=10)
endif

! 
! !write file to  qui4.TXT
if (n4==0) then
n4 = 0
else
open( unit = 11, file ='qui4.txt',status = 'replace' )
do i=1,n4
write(11,FMT)qu4(i,:)
enddo
close(unit=11)
endif

! !write file to  qui5.TXT
if (n5==0) then
n5 = 0
else
open( unit = 11, file ='qui5.txt',status = 'replace' )
do i=1,n5
write(11,FMT)qu5(i,:)
enddo
close(unit=11)
endif

! !write file to  hyster.TXT
open( unit=11,file='hyster.txt',status = 'replace')
do i=1,ndomain
write(11,FMT2) E(i),sum(PxD3((R+1):(L+R),(R+1):(L+R),i))/L/L
enddo
close(unit=11)



!create gnuplot command file
OPEN(10,ACCESS='SEQUENTIAL',FILE='gp3.gnu')
! write(10,*) 'set terminal png size 1800,1800 enhanced font "Helvetica,20" crop'
write(10,*) 'set term epslatex size 3,3'
!***************************************************
!***************************************************
!***************************************************
! begin multiplot
if (ndomain<10) then
write(10,"(A18,I1,A5)") 'set output "jj.000',ndomain,'.eps"'
elseif ( (ndomain<100) .and. (ndomain>=10) ) then
write(10,"(A17,I2,A5)") 'set output "jj.00',ndomain,'.eps"'
elseif ( (ndomain<1000) .and. (ndomain>=100) ) then
write(10,"(A16,I3,A5)") 'set output "jj.0',ndomain,'.eps"'
else 
write(10,"(A15,I4,A5)") 'set output "jj.',ndomain,'.eps"'
endif





!***************************************************
! plot domain
write(10,*) 'unset key'
write(10,*) 'unset xtics'
write(10,*) 'unset ytics'
write(10,*) 'set tmargin at screen 1'
write(10,*) 'set bmargin at screen 0'
write(10,*) 'set lmargin at screen 0'
write(10,*) 'set rmargin at screen 1'
write(10,*) 'set border 15 lw 5'

write(10,*) 'set style arrow 1 head filled size screen 0.015,20,40 linecolor rgb "black" linewidth 0.2'
write(10,*) 'set style arrow 2 head filled size screen 0.015,20,40 linecolor rgb "white" linewidth 0.8'
write(10,*) 'set style arrow 3 head filled size screen 0.015,20,40 linecolor rgb "black" linewidth 0.2'
write(10,*) 'set style arrow 4 head filled size screen 0.015,20,40 linecolor rgb "black" linewidth 0.2'
write(10,*) 'set style arrow 5 head filled size screen 0.015,20,40 linecolor rgb "white" linewidth 0.8'
write(10,*) 'set object 1 rectangle from 9,9 to 72,72  fs empty border rgb "gold" lw 3'  

write(10,*) 'plot ''qui1.txt'' using 1:2 with points pt 7 ps 1.6 linecolor rgb "red",&
	    ''qui2.txt'' using 1:2 with points pt 7 ps 1.6 linecolor rgb "blue",& 
	    ''qui3.txt'' using 1:2 with points pt 7 ps 1.6 linecolor rgb "green",&
            ''qui4.txt'' using 1:2 with points pt 7 ps 1.6 linecolor rgb "yellow",&
            ''qui5.txt'' using 1:2 with points pt 7 ps 1.6 linecolor rgb "black",& 
            ''qui1.txt'' using 1:2:3:4 with vectors arrowstyle 1, &   
            ''qui2.txt'' using 1:2:3:4 with vectors  arrowstyle 2, &
            ''qui3.txt'' using 1:2:3:4 with vectors  arrowstyle 3, &
            ''qui4.txt'' using 1:2:3:4 with vectors  arrowstyle 4, &
            ''qui5.txt'' using 1:2:3:4 with vectors  arrowstyle 2 '
            
write(10,*) 'unset multiplot'

CLOSE(10,STATUS='KEEP')

! ! END OF DOMAIN CONFIGURATION FIGURE
!***************************************************
!***************************************************
! ret=SYSTEM('gnuplot gp3.gnu')



  deallocate(qu1)
  deallocate(qu2)
  deallocate(qu3)
  deallocate(qu4)
enddo
end