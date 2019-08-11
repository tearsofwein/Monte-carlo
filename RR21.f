  






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
parameter( L = 63 , R = 8 , num_loc = 99 )
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
parameter( SJ = 2.0 , SE = 1.0 , SAnr = 6.8012 , SA2 = 6.8012 , SAr = 6.8012 )


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
parameter(kT = 2.8 )

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


!***************************************************
!***************************************************
module mgeneration

use global
contains
subroutine mtotal4(steps,m_choose)

integer :: site_choose,mtotal1(L*L/ ((R+1)* (R+1)),2),mt1,nt1


real(kind=8) :: rnd1,rnd2,rnd3,rnd4
! site_x is the how many sites are independent.
integer :: site_x,site_y
integer :: q,k,i,j,li

! number index of independent sites
integer :: II,site_number(L*L/(R+1)/(R+1)),steps,temp_site_number,site_sum

! coordinates of independent sites
integer :: m_choose(Judgestep * t / threads,threads,2)
steps = Judgestep * t / threads
call init_random_seed()

! Calculate the steps and how many site in one dimension
site_x = L/ (R+1)
site_y = L/ (R+1)

site_sum = site_x * site_y


! generate all the possible site number
do II = 1 , site_sum
site_number(II) = II
end do



! Dependening on the first 'm' to generate the rest 'm' in the whole lattice

! Pass 'm' to each thread




do q = 1, steps

     ! generate the initial site coordinates
      call random_number(rnd1)
      call random_number(rnd2)
      mt1 = rnd1*(R +1 ) + 1 + R
      nt1 = rnd2*(R +1 ) + 1 + R


      ! at the same electric field, in different theads  judge whether the site
      ! number is repeated. If it happens, regenerate the site number
      do II = 1, site_sum
         call random_number(rnd3)
         site_choose = rnd3 * site_x * site_y + 1
         temp_site_number = site_number(site_choose)
         site_number(site_choose) = site_number(II)
         site_number(II) = temp_site_number
      end do

do k = 1,threads

      ! obtain all the non-interactive parts
      do i = 1,site_x
        do j = 1,site_y
            mtotal1(j+(i-1)*site_y,1) = mt1 + ( R + 1 ) * ( i -1 )
            mtotal1(j+(i-1)*site_y,2) = nt1 + ( R + 1 ) * ( j -1 )
        end do
      end do



      ! We choose the site,which is not repeated in site_number serial, so that
      ! we can pass it to the site_number serial


      ! pass unrepeated site number to each thead
      m_choose(q,k,1) = mtotal1(site_number(k),1)
      m_choose(q,k,2) = mtotal1(site_number(k),2)


end do
end do
return
end subroutine

end module
!***************************************************
!********************************************f*******





! ***************************************************
! ***************************************************
module gaussianrndGenerator
use global
contains
function gaussianrnd( mean1 ,deviation1 )
! ***************************************************
! definition of pi should be done in a module
real(kind=8) :: random_number1 , random_number2 ! two random numbers with uniform distribution

real(kind=4) :: mean1        ! mean of the gaussian distribution we want a random number from

real(kind=4) :: deviation1   ! standard_deviation of the gaussian distribution

real(kind=8) :: gaussianrnd

! first step, get two random numbers with uniform distribution. This is done easily with the fortran 90 intrinsic function random_number

call random_number ( random_number1 )

call random_number ( random_number2 )

! convert to normal distribution, which is the purpose of this code :)
gaussianrnd = deviation1 * sqrt ( -2.0 * log ( random_number1 ) ) * cos ( 2.0 * pi * random_number2 ) + mean1
return
end function gaussianrnd

function Periodic(Pinput)
real(kind=8) :: Pinput(L,L),Periodic( L+2*R , L+2*R )
Periodic( (R+1):(R+L), (R+1):(R+L) ) = Pinput
Periodic( (R+1):(R+L) , 1:R ) = Periodic( (R+1):(R+L) ,(L+1):(R+L) )
Periodic( (R+1):(R+L) , (R+L+1):(L+2*R) ) = Periodic( (R+1):(R+L) ,(R+1):(2*R) )
Periodic( 1:R, : ) = Periodic( (L+1):(R+L),: )
Periodic( (R+L+1):(L+2*R) , : ) = Periodic( (R+1):(2*R) , : )
return
end function Periodic

end module
! ***************************************************
! ***************************************************






! ***************************************************
! ***************************************************
module EnergyOfNodes
use global
contains

!***************************************************
function EnergyHe()

real(kind=8) :: EnergyHe
! Electrostatic energy including one node
!if (Ex == 0.0) then
!EnergyHe = 0.0
!else
EnergyHe =  - sum( SE * PxTotal * Ex  ) 
!end if
return
end function EnergyHe

!***************************************************
function EnergyHfld()
real(kind=8) :: EnergyHfld
! Landau double-well potential

EnergyHfld = sum( &
  ( &
    a  * ( PxTotal**2.0 + PyTotal**2.0 ) + b * ( PxTotal**4.0 + PyTotal**4.0 ) &
  ) &
  * ( 1.0 - rndmap) + &
  (  &
    a1 * ( PxTotal**2.0 + PyTotal**2.0 ) + b1* ( PxTotal**4.0 + PyTotal**4.0 ) &
  ) &
  * rndmap &
) 

return
end function EnergyHfld

!***************************************************
function EnergyHdip()

real(kind=8) :: EnergyHdip
real(kind=8) :: Pxjk,Pyjk,fdipjk( L , L )
real(kind=8) :: fdip( 2*R+1 , 2*R+1 )
integer :: range11,range12,range21,range22
integer :: j2,k2
real(kind=8) :: Px2( 2*R+1 , 2*R+1 ),Py2( 2*R+1 , 2*R+1 )
!***************************************************
! Dipole-dipole interaction


do j2=RR,RL
do k2=RR,RL

range11 = j2-R
range12 = j2+R
range21 = k2-R
range22 = k2+R



Px2 = Px( range11:range12 ,range21:range22 )
Py2 = Py( range11:range12 ,range21:range22 )

Pxjk =  Px2(R+1,R+1)
Pyjk =  Py2(R+1,R+1)

! interaction energy
fdip =  ( Pxjk * Px2  +  Pyjk * Py2 ) / dists**3.0 &
-3* ( Pxjk * dxcoords + Pyjk * dycoords ) * ( Px2 * dxcoords +  Py2 * dycoords ) / dists**5.0
fdip(RR,RR) = 0.0

if (Ex_rnd(j2,k2)==1.0) then
fdip = fdip * Ex_rnd(range11:range12 ,range21:range22 ) * SAr &
 + fdip * ( 1.0- Ex_rnd(range11:range12 ,range21:range22 ) ) * SA2
else
fdip = fdip * ( 1.0 - Ex_rnd(range11:range12 ,range21:range22 ) ) * SAnr &
 + fdip * Ex_rnd(range11:range12 ,range21:range22 ) * SA2
endif


fdipjk(j2-R,k2-R) = sum(fdip)

end do
end do

EnergyHdip = sum(fdipjk) / 2.0
return
end function EnergyHdip

!***************************************************
function EnergyHgr()



real(kind=8) :: Hgrx(L,L),Hgry(L,L),pxxT(L,L),pyyT(L,L),pxyT(L,L),pyxT(L,L)

! Gradient of polarization field
if (SJ == 0.0) then
Hgr = 0.0
else

pxxT = ( Px( RR:RL , RRR:RLL ) - Px(RR:RL,RR:RL) )/ det_x

pyyT = ( Py( RRR:RLL , RR:RL ) - Py(RR:RL,RR:RL) )/ det_y

pxyT = ( Px( RRR:RLL , RR:RL ) - Px(RR:RL,RR:RL) )/ det_y

pyxT = ( Py( RR:RL , RRR:RLL ) - Py(RR:RL,RR:RL) )/ det_x


Hgrx = pxxT **2 + g2 *  pyxT **2 + g22 * pyxT **2
Hgry = pyyT **2 + g2 *  pxyT **2 + g22 * pxyT **2

EnergyHgr = ( sum(Hgrx) + sum(Hgry) ) * SJ / 2.0

end if
return
end function EnergyHgr







!***************************************************
!***************************************************
function EnergyOfSelectNodes(m,n,Pxmn,Pymn,Px2,Py2)

! energy(Electrostatic energy,Landau double-well potential,Dipole-dipole interaction,Gradient of polarization field,total energy)
real(kind=8) :: fEmn,fldmn,fdipmn,fgrmn
real(kind=8) :: fgrx(2),fgry(2),pxx(2),pyy(2),pxy(2),pyx(2)

real(kind=8) :: fdip( 2*R+1 , 2*R+1 ),Px2( 2*R+1 , 2*R+1 ),Py2( 2*R+1 , 2*R+1 )

real(kind=8) :: Pxmn,Pymn

real(kind=8) :: EnergyOfSelectNodes

integer  :: m,n,range11,range12,range21,range22
!***************************************************

! Electrostatic energy including one node
!if (Ex == 0.0) then
!fEmn = 0.0
!else
fEmn = - SE *  Pxmn * Ex 

!! Dipole-dipole interaction
! dipole interaction energy
range11 = m-R
range12 = m+R
range21 = n-R
range22 = n+R
fdip =   ( Pxmn * Px2  +  Pymn * Py2 ) / dists**3.0 &
-3.0*( Pxmn * dxcoords + Pymn * dycoords ) * ( Px2 * dxcoords +  Py2 * dycoords )  / dists**5.0
! Landau double-well potential
if (Ex_rnd(m,n)==1.0) then
fldmn = a1  * ( Pxmn**2 + Pymn**2 )  + b1 * ( Pxmn**4 + Pymn**4 ) 
!***************************************************
! Dipole-dipole interaction
! dipole interaction energy
fdip = fdip * Ex_rnd(range11:range12 ,range21:range22 ) * SAr &
 + fdip * ( 1.0 - Ex_rnd(range11:range12 ,range21:range22 ) ) * SA2

fdip(RR,RR) = 0.0
fdipmn = sum(fdip)




else
fldmn = a  * ( Pxmn**2 + Pymn**2 )  + b * ( Pxmn**4 + Pymn**4 ) 
!***************************************************
! Dipole-dipole interaction
! dipole interaction energy
fdip = fdip * ( 1.0 - Ex_rnd(range11:range12 ,range21:range22 ) ) * SAnr &
 + fdip * Ex_rnd(range11:range12 ,range21:range22 ) * SA2
fdip(RR,RR) = 0.0
fdipmn = sum(fdip)
endif




!***************************************************
! Gradient of polarization field
if (SJ == 0.0) then
fgrmn = 0.0
else

pxx = ( Px2( R+1 , (R+1):(R+2) ) - Px2( R+1 , R:(R+1) ) )/ det_x

pyy = ( Py2( (R+1):(R+2) , R+1 ) - Py2( R:(R+1) , R+1 ) )/ det_y

pxy = ( Px2( (R+1):(R+2) , R+1 ) - Px2( R:(R+1) , R+1 ) )/ det_y

pyx = ( Py2( R+1 , (R+1):(R+2) ) - Py2( R+1 , R:(R+1) ) )/ det_x



fgrx = Pxx **2 + g2 *  Pyx **2 + g22 * Pyx **2
fgry = Pyy **2 + g2 *  Pxy **2 + g22 * Pxy **2

fgrmn = ( sum(fgrx) + sum(fgry) ) *SJ / 2.0


end if

! Total Energy
EnergyOfSelectNodes = fldmn  + fdipmn + fgrmn + fEmn

return
end function EnergyOfSelectNodes

! The end of function for energy calculation of selected nodes
!***************************************************
!***************************************************

end module
! The end of  for energy calculation of whole system
!****************************************************
!****************************************************
!****************************************************
 


! ***************************************************
! ***************************************************
! call the data and time in computer to generate random number
subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms

            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
!            open(newunit=un, file="/dev/urandom", access="stream", &
!                 form="unformatted", action="read", status="old", iostat=istat)
!            if (istat == 0) then
!               read(un) seed
!               close(un)
!            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
!            end if
            call random_seed(put=seed)
          end subroutine init_random_seed
! The end of  random seed
!****************************************************
!****************************************************
!****************************************************




!***************************************************
!***************************************************
Program Re6

use global
use gaussianrndGenerator
use mgeneration
use EnergyOfNodes


include 'omp_lib.h'



real(kind=8) :: Px2( 2*R+1 , 2*R+1 ),Py2( 2*R+1 , 2*R+1 )
integer :: m,n,step,range11s,range12s,range22s,range21s

integer :: direction(4,2)
data direction /1,-1,0,0,0,0,1,-1/

real(kind=8) :: M1(L,L),N1(L,L),EM1(L,L),EN1(L,L)
integer :: xdirection(L,L),ydirection(L,L),Exdirection(L,L),Eydirection(L,L)
integer :: media(L,L),variant,Emedia(L,L)

! random number for selecting site
real(kind=8) :: rnd1,rnd2,rnd3


!  electric field
integer(kind=8) :: Et(Estep+1)= (/(I,I=0,Estep)/)
real(kind=8) :: Emax, E(Estep+1)

! Energy of choosed nodes : H0 - previous one and Hnew - current one after switching. Random switching dipole value. detaH - differences between the previous and current energy
real(kind=8) :: H0,Hnew,Pmnrnd,detaH,Hnew_0

! energy
real(kind=8) :: HeProcess( Estep+1 ),HfldProcess( Estep+1 ),HdipProcess( Estep+1 )
real(kind=8) :: HgrProcess( Estep+1 ),HtotalProcess( Estep+1 )


! probability
real(kind=8) :: prob

! new polarization fields
real(kind=8) :: Pxmnold, Pymnold,Pxmnnew,Pymnnew,Pxmnnew_0=0.0,Pymnnew_0=0.0
real(kind=8) :: PxProcess(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: PyProcess(2*R+L , 2*R+L , Estep+1 )
real(kind=8) :: Energy_profile_Process(2*R+L , 2*R+L , Estep+1 )

integer :: aaa,bbb

integer :: steps
integer :: m_choose(Judgestep * t / threads,threads,2)

! The default value for cpu_time_type is used
REAL(kind=8) T1, T2 , ave_He

! RandomNumber groups
integer :: mtotal5(Judgestep * t / threads,threads,2) , ntotal5(Judgestep * t / threads,threads,2)
!integer :: m_choose(Judgestep * t / threads,threads,2) , m_choose(Judgestep * t / threads,threads,2)
integer :: mtotal6(Judgestep * t , 1  ) , ntotal6(Judgestep * t ,1 )

integer(kind=8) :: mm(L*L,2),aa,me(1,2)

! end of declaration
!***************************************************
! For some functions
! set the number of used threads
CALL OMP_SET_NUM_THREADS(threads)
! get the time
CALL CPU_TIME(T1)
call init_random_seed()

write(*,*)OMP_GET_MAX_THREADS()
!***************************************************
! some invariables
RR = R+1
RRR = R+2
RL = R+L
RLL = R+L+1
rangetotal = 2*R+1


steps = Judgestep * t / threads
!***************************************************
!  electric field
Emax = 10.0
E = Emax * sin( 2.0*pi*f*Et )

!***************************************************



! x,y coordinates
do I=1,(L+2*R),1
xcoord(:,I)=I*det_x
ycoord(I,:)=I*det_y
end do


!***************************************************
! save parameter
open( unit = 14, file ='L1.bin',status = 'replace' , access = 'direct'  , recl = 4 )
write(unit=14, rec =1)L
close(unit=14)
open( unit = 16, file ='R1.bin',status = 'replace' ,access = 'direct'  , recl = 4 )
write(unit=16,rec=1)R
close(unit=16)


open( unit = 14, file ='KT1.bin',status = 'replace' , access = 'direct'  , recl = 8 )
write(unit=14, rec =1)kt
close(unit=14)

open( unit = 14, file ='SJ1.bin',status = 'replace' , access = 'direct'  , recl = 8 )
write(unit=14, rec =1)SJ
close(unit=14)


! save the xcoord and ycoord data in file 'axis.bin'
open( unit = 10, file ='xaxis1.bin' ,status = 'replace'  ,access = 'direct'  , recl =  (L+2*R) * (L+2*R)*4 )
write(unit=10,rec=1)xcoord
close(unit=10)

open( unit = 12, file ='yaxis1.bin' ,status = 'replace' , access = 'direct'  , recl =  (L+2*R) * (L+2*R) *4 )
write(unit=12,rec=1)ycoord
close(unit=12)



!***************************************************
!call random_seed()
call random_number(M1)
call random_number(N1)
!write(*,*)M1,N1
! generate random integer matrix
media = 4*M1+1

! polarization direction variants in x axis and y axis
do I=1,L
do J=1,L
xdirection(I,J) = direction (  media(I,J) ,1 )
ydirection(I,J) = direction (  media(I,J) ,2 )
end do
end do


!***************************************************
! polarization in x axis and y axis
Px((R+1):(R+L) ,(R+1):(R+L) ) = xdirection * N1
Py((R+1):(R+L) ,(R+1):(R+L) ) = ydirection * N1





do i=1,L+2*R,1
do j=1,L+2*R,1
if (Px(i,j)/=0.0 .and. Py(i,j)/=0.0) then
write(*,*)'very old wrong'
end if
end do
end do

!WRITE(*,*)Px,Py

!Px = 0.0
!Py = 0.0


!***************************************************
! save the PxD and PyD data in file
open( unit = 10, file ='PxD1.bin' ,status = 'old'  ,access = 'direct'  , recl =  (L+2*R) * (L+2*R)*8 )
read(unit=10,rec=1)Px
close(unit=10)

open( unit = 12, file ='PyD1.bin' ,status = 'old' , access = 'direct'  , recl =  (L+2*R) * (L+2*R) *8 )
read(unit=12,rec=1)Py
close(unit=12)


!***************************************************

! Generate all the necessary random numbers and reshape them


call mtotal4(steps,m_choose)


!mtotal6 = reshape(mtotal5, (/ Judgestep * t , 1 /))
!ntotal6 = reshape(ntotal5, (/ Judgestep * t , 1 /))

!***************************************************
! Distance between the select node and other nodes for calculation of the dipole-dipole interaction
dxcoords = xcoord(R+1,R+1) - xcoord( 1:(2*R+1) ,1:(2*R+1)  )
dycoords = ycoord(R+1,R+1) - ycoord( 1:(2*R+1) ,1:(2*R+1)  )
dist2s = dxcoords**2 + dycoords**2
dists = reshape( sqrt( reshape(dist2s,(/rangetotal**2,1/))), (/rangetotal,rangetotal/) )
dists(RR,RR)=1

!***************************************************

! Boundary condition for random electric field


call random_number(EM1)


!***************************************************
! polarization in x axis and y axis
! random angle
EM1 = 2.0 * pi * EM1

!write(*,*)EN1
! For ALL dipole influenced by random field
call random_number(EN1)





!***************************************************
!***************************************************
!decide how many percent without random field



! initial sequent array
do i=1,L
do j=1,L
mm(j+(i-1)*L,1) = i
mm(j+(i-1)*L,2) = j
enddo
enddo
! write(*,*)mm(63+54,1)
! write(*,*)mm(63+54,2)



! random array
do i=1,L
do j=1,L
call random_number(rnd1)
aa = (L*L) * rnd1 + 1

me(1,1) = mm(aa,1)
me(1,2) = mm(aa,2)

mm(aa,1) =  mm(j+(i-1)*L,1)
mm(aa,2) = mm(j+(i-1)*L,2)

mm(j+(i-1)*L,1) = me(1,1)
mm(j+(i-1)*L,2) = me(1,2)

enddo
enddo



do i = 1,num_loc




end do




! Periodic boundary condition for Ex_Rnd and Ey_Rnd



rndmap = Ex_rnd( RR:RL , RR:RL )


! save the PxD and PyD data in file
open( unit = 10, file ='Exrnd.bin' ,status = 'old'  ,access = 'direct'  , recl =  (L+2*R) * (L+2*R)*8 )
read(unit=10,rec=1)Ex_Rnd
close(unit=10)




! Periodic boundary condition for Px and Py



!***************************************************
!***************************************************


E(1:(Estep/9*2)) = 3.2
E((Estep/9*2+1):(Estep/9*4)) = 0.0
E((Estep/9*4):(Estep/3*2)) = 3.2
E( (Estep/3*2+1):(Estep+1) ) = 0.0


open( unit = 10, file ='E1.bin' ,status = 'replace'  ,access = 'direct'  , recl =  (Estep+1)*8 )
write(unit=10,rec=1)E
close(unit=10)






do k=1,Estep/3*2

Ex = E(k)


!***************************************************
!***************************************************
! begin parallel process


do j=1,(Judgestep * t / threads)
!$omp parallel
!$omp do schedule(static, 1) private(i,m,n,H0,variant,Pmnrnd,Pxmnold,Pymnold,Hnew,detaH,rnd1 &
!$omp , rnd2,rnd3,prob,range11s,range12s,range21s,range22s,Pxmnnew,Pymnnew,Px2,Py2,fdip,fdipmn,fgrmn &
!$omp , fEmn,fldmn,pxx,pyy,pxy,pyx,fgrx,fgry ,Pxs ,Pys &
!$omp , random_number1 , random_number2 , mean1 , deviation1 ,Hnew_0 )

do i=1,threads


!***************************************************
! shared value of dipoles before and after switching

! Select the node randomly
m = m_choose( j ,  i , 1 )
n = m_choose( j ,  i , 2 )
!write(*,*)m,n

if (Ex_Rnd(m,n)==0.0)then

! rangeset for selected nodes4
range11s = m-R
range12s = m+R
range21s = n-R
range22s = n+R

!***************************************************
! Judge whether the data transfer process is correct


!***************************************************
! Copy Px and Py to Pxs and Pys
Pxs = Px
Pys = Py


!***************************************************
! Boundary conditions
Pxs(  RR:RL , 1:R ) = Pxs( RR:RL ,(L+1):RL )
Pxs(  RR:RL , RLL:(L+2*R) ) = Pxs( RR:RL ,RR:(2*R) )
Pxs( 1:R, : ) = Pxs( (L+1):RL,: )
Pxs( RLL:(L+2*R) , : ) = Pxs( RR:(2*R) , : )

!***************************************************
Pys(  RR:RL , 1:R ) = Pys( RR:RL ,(L+1):RL )
Pys(  RR:RL , RLL:(L+2*R) ) = Pys( RR:RL ,RR:(2*R) )
Pys( 1:R, : ) = Pys( (L+1):RL,: )
Pys( RLL:(L+2*R) , : ) = Pys( RR:(2*R) , : )



! the necessary Px and Py for calculation the selected nodes
Px2 = Pxs( range11s:range12s ,range21s:range22s )
Py2 = Pys( range11s:range12s ,range21s:range22s )


!***************************************************
! Energy of choosed nodes : previous one

!  polarization before switching
Pxmnold = Px2(R+1,R+1)
Pymnold = Py2(R+1,R+1)


H0 = EnergyOfSelectNodes(m,n,Pxmnold,Pymnold,Px2,Py2)
!***************************************************




!***************************************************
! Energy after switching

! random switching dipole value
call random_number(Pmnrnd)


! random switching variants
call random_number(rnd1)
variant = rnd1*4+1

! random switching dople
Pxmnnew = direction(variant,1) * Pmnrnd
Pymnnew = direction(variant,2) * Pmnrnd

! change Px2 and Py2 for
Px2(R+1,R+1) = Pxmnnew
Py2(R+1,R+1) = Pymnnew

! Energy of choosen nodes : after switching
Hnew = EnergyOfSelectNodes(m,n,Pxmnnew,Pymnnew,Px2,Py2)
!***************************************************

!***************************************************
! do monte carlo (refusal or approval)
! generate judge probability
call random_number(rnd2)
!WRITE(*,*)rnd2,rnd3
! detaH - differences between the previous and current energy
detaH = Hnew - H0
! switching probability
prob = exp(-detaH/kT)


!***************************************************
! judge whether switch. If the condition is not fulfilled, return to old value of Px and Py.

if ( (detaH<0.0 .or. prob > rnd2) )then
! change Px and Py
Px(m,n) = Pxmnnew
Py(m,n) = Pymnnew
!Energy_profile(m,n) = Hnew
endif


endif


end do
!$omp end do
!$omp end parallel
! End parallel process
!***************************************************
!***************************************************

end do





!***************************************************
!***************************************************
! Boundary condition

Px(  (R+1):(R+L) , 1:R ) = Px( (R+1):(R+L) ,(L+1):(R+L) )
Px(  (R+1):(R+L) , (R+L+1):(L+2*R) ) = Px( (R+1):(R+L) ,(R+1):(2*R) )
Px( 1:R, : ) = Px( (L+1):(R+L),: )
Px( (R+L+1):(L+2*R) , : ) = Px( (R+1):(2*R) , : )

Py(  (R+1):(R+L) , 1:R ) = Py( (R+1):(R+L) ,(L+1):(R+L) )
Py(  (R+1):(R+L) , (R+L+1):(L+2*R) ) = Py( (R+1):(R+L) ,(R+1):(2*R) )
Py( 1:R, : ) = Py( (L+1):(R+L),: )
Py( (R+L+1):(L+2*R) , : ) = Py( (R+1):(2*R) , : )

!Energy_profile(  (R+1):(R+L) , 1:R ) = Energy_profile( (R+1):(R+L) ,(L+1):(R+L) )
!Energy_profile(  (R+1):(R+L) , (R+L+1):(L+2*R) ) = Energy_profile( (R+1):(R+L) ,(R+1):(2*R) )
!Energy_profile( 1:R, : ) = Energy_profile( (L+1):(R+L),: )
!Energy_profile( (R+L+1):(L+2*R) , : ) = Energy_profile( (R+1):(2*R) , : )

!***************************************************
!***************************************************
! Px and Py in every different electric field
PxProcess(: , : , k) = Px
PyProcess(: , : , k) = Py
!Energy_profile_process(: , : , k) = Energy_profile


! Calculation of total energy



!***************************************************
PxTotal =  Px( RR:RL , RR:RL )
PyTotal =  Py( RR:RL , RR:RL )




HeProcess(k) = EnergyHe()
HfldProcess(k) = EnergyHfld()
HdipProcess(k) = EnergyHdip()
HgrProcess(k) = EnergyHgr()
HtotalProcess(k) = HeProcess(k)+HfldProcess(k)+HdipProcess(k)+HgrProcess(k)
!HtotalProcess(k) = sum( Energy_profile( (R+1):(R+L) ,(R+1):(R+L) ) )
!write(*,*)'locmn'
CALL CPU_TIME(T2)
write(*,*) 'Time for second chunk of code was ', (T2-T1)/threads, 'seconds.'
end do


ave_He = HtotalProcess(Estep/3)/L**2









! to get better statistics, calibrate the result
Demon = kT + ( sum( HtotalProcess(1001:1200) ) / 200 - HtotalProcess(1200) ) / (L*L*df)
DemonProcess(: , : , 1:(Estep/3*2) ) = df*kT




do k= (Estep/3*2+1),(Estep+1)

Ex = E(k)


!***************************************************
!***************************************************
! begin parallel process


do j=1,(Judgestep * t / threads)
!$omp parallel
!$omp do schedule(static, 1) private(i,m,n,H0,variant,Pmnrnd,Pxmnold,Pymnold,Hnew,detaH,rnd1 &
!$omp , rnd2,rnd3,prob,range11s,range12s,range21s,range22s,Pxmnnew,Pymnnew,Px2,Py2,fdip,fdipmn,fgrmn &
!$omp , fEmn,fldmn,pxx,pyy,pxy,pyx,fgrx,fgry ,Pxs ,Pys  &
!$omp , random_number1 , random_number2 , mean1 , deviation1 ,Hnew_0 )

do i=1,threads


!***************************************************
! shared value of dipoles before and after switching

! Select the node randomly
m = m_choose( j ,  i , 1 )
n = m_choose( j ,  i , 2 )
!write(*,*)m,n

if (Ex_Rnd(m,n)==0.0)then

! rangeset for selected nodes4
range11s = m-R
range12s = m+R
range21s = n-R
range22s = n+R

!***************************************************
! Judge whether the data transfer process is correct


!***************************************************
! Copy Px and Py to Pxs and Pys

Pxs = Px
Pys = Py

!if (Px(m,n)/=0.0 .and. Py(m,n)/=0.0) then
!write(*,*)'input,wrong',i,j,m,n,Px(m,n),Py(m,n)
!endif
!***************************************************
! Boundary conditions
Pxs(  RR:RL , 1:R ) = Pxs( RR:RL ,(L+1):RL )
Pxs(  RR:RL , RLL:(L+2*R) ) = Pxs( RR:RL ,RR:(2*R) )
Pxs( 1:R, : ) = Pxs( (L+1):RL,: )
Pxs( RLL:(L+2*R) , : ) = Pxs( RR:(2*R) , : )

!***************************************************
Pys(  RR:RL , 1:R ) = Pys( RR:RL ,(L+1):RL )
Pys(  RR:RL , RLL:(L+2*R) ) = Pys( RR:RL ,RR:(2*R) )
Pys( 1:R, : ) = Pys( (L+1):RL,: )
Pys( RLL:(L+2*R) , : ) = Pys( RR:(2*R) , : )



! the necessary Px and Py for calculation the selected nodes
Px2 = Pxs( range11s:range12s ,range21s:range22s )
Py2 = Pys( range11s:range12s ,range21s:range22s )


!***************************************************
! Energy of choosed nodes : previous one

!  polarization before switching
Pxmnold = Px2(R+1,R+1)
Pymnold = Py2(R+1,R+1)

! Judge whether the data transfer process is correct
!if (Pxmnold/=0.0 .and. Pymnold/=0.0) then
!write(*,*)'Pass from original to calculate the old energy, wrong',i,j,m,n,Pxmnold,Pymnold
!endif
!if (k== Estep/3+1) then
!H0 = EnergyOfSelectNodes(m,n,Pxmnold,Pymnold,Px2,Py2)
!!***************************************************
!else
H0 = EnergyOfSelectNodes(m,n,Pxmnold,Pymnold,Px2,Py2)
!endif

!H0 = Energy_profile(m,n)

!***************************************************
! Energy after switching

! random switching dipole value
call random_number(Pmnrnd)

! random switching variants
call random_number(rnd1)
variant = rnd1*4+1

! random switching dople
Pxmnnew = direction(variant,1) * Pmnrnd
Pymnnew = direction(variant,2) * Pmnrnd

! change Px2 and Py2 for
Px2(R+1,R+1) = Pxmnnew
Py2(R+1,R+1) = Pymnnew

! Energy of choosen nodes : after switching
Hnew = EnergyOfSelectNodes(m,n,Pxmnnew,Pymnnew,Px2,Py2)
!***************************************************

!***************************************************
! do monte carlo (refusal or approval)
! generate judge probability
call random_number(rnd2)
!WRITE(*,*)rnd2,rnd3
! detaH - differences between the previous and current energy
detaH = Hnew - H0

!***************************************************
! judge whether switch. If the condition is not fulfilled, return to old value of Px and Py.

variant = df * rnd2 + 1
if ( Demon(m,n,variant)>detaH )then
! change Px and Py
Px(m,n) = Pxmnnew
Py(m,n) = Pymnnew
Demon(m,n,variant) = Demon(m,n,variant) - detaH
endif




endif

end do ! i=1,threads
!$omp end do
!$omp end parallel
! End parallel process
!***************************************************
!***************************************************
! Calculate the demon sum for every step
DemonSum( j , k ) = sum ( Demon( RR:RL , RR:RL , :) ) / (L * L) / df

end do  ! j=1,(Judgestep * t / threads)

!***************************************************
!***************************************************
! Boundary condition

Px(  (R+1):(R+L) , 1:R ) = Px( (R+1):(R+L) ,(L+1):(R+L) )
Px(  (R+1):(R+L) , (R+L+1):(L+2*R) ) = Px( (R+1):(R+L) ,(R+1):(2*R) )
Px( 1:R, : ) = Px( (L+1):(R+L),: )
Px( (R+L+1):(L+2*R) , : ) = Px( (R+1):(2*R) , : )

Py(  (R+1):(R+L) , 1:R ) = Py( (R+1):(R+L) ,(L+1):(R+L) )
Py(  (R+1):(R+L) , (R+L+1):(L+2*R) ) = Py( (R+1):(R+L) ,(R+1):(2*R) )
Py( 1:R, : ) = Py( (L+1):(R+L),: )
Py( (R+L+1):(L+2*R) , : ) = Py( (R+1):(2*R) , : )


!***************************************************
!***************************************************
! Px and Py in every different electric field
PxProcess(: , : , k) = Px
PyProcess(: , : , k) = Py
DemonProcess(: , : , k) = sum(Demon,dim=3)



! Calculation of total energy



!***************************************************
PxTotal =  Px( RR:RL , RR:RL )
PyTotal =  Py( RR:RL , RR:RL )





HeProcess(k) = EnergyHe()
HfldProcess(k) = EnergyHfld()
HdipProcess(k) = EnergyHdip()
HgrProcess(k) = EnergyHgr()
HtotalProcess(k) = HeProcess(k)+HfldProcess(k)+HdipProcess(k)+HgrProcess(k)
!HtotalProcess(k) = sum( Energy_profile( (R+1):(R+L) ,(R+1):(R+L) ) )
!write(*,*)'locmn'
CALL CPU_TIME(T2)
write(*,*) 'Time for second chunk of code was ', (T2-T1)/threads, 'seconds.'
end do




write(*,*)HtotalProcess(Estep)

!***************************************************
! save the PxD, PyD, He, Hfld, Hdip, Hgr, Htotal data in files
open( unit = 10, file ='PxProcess1.bin' ,status = 'replace'  ,access = 'direct'  , recl =  (L+2*R) * (L+2*R) * (Estep+1) *8 )
write(unit=10,rec=1)PxProcess
close(unit=10)

open( unit = 12, file ='PyProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (L+2*R) * (L+2*R) * (Estep+1) *8 )
write(unit=12,rec=1)PyProcess
close(unit=12)

open( unit = 12, file ='DemonProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (L+2*R) * (L+2*R) * (Estep+1) *8 )
write(unit=12,rec=1)DemonProcess
close(unit=12)


open( unit = 12, file ='HeProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (Estep+1) *8 )
write(unit=12,rec=1)HeProcess
close(unit=12)

open( unit = 12, file ='DemonSum1.bin' ,status = 'replace' , access = 'direct'  , recl =  (Estep+1) * Judgestep * t / threads  *8 )
write(unit=12,rec=1)DemonSum
close(unit=12)


open( unit = 12, file ='HfldProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (Estep+1) *8 )
write(unit=12,rec=1)HfldProcess
close(unit=12)

open( unit = 12, file ='HdipProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (Estep+1) *8 )
write(unit=12,rec=1)HdipProcess
close(unit=12)

open( unit = 12, file ='HgrProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (Estep+1) *8 )
write(unit=12,rec=1)HgrProcess
close(unit=12)


open( unit = 12, file ='HtotalProcess1.bin' ,status = 'replace' , access = 'direct'  , recl =  (Estep+1) *8 )
write(unit=12,rec=1)HtotalProcess
close(unit=12)


stop
end program  Re6

! The end of the whole program
!***************************************************
!***************************************************
