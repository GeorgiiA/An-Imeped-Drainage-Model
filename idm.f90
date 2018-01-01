program main
  character (len=3):: hi
  character (len=5) :: str
  integer :: n !number of x steps
  integer :: m !number of time steps
  integer :: k ! counter for time steps
  integer :: nyrs ! number of years, years
  integer :: i ! counter for x steps
  integer :: narg ! reading the number of arguments in command 
  real :: tau  ! time step
  real :: xstep ! x step
  real :: U ! the equivalent effective rainfall, in km of water layer per year
  real :: u0 ! effective rainfall, in mm of water layer
  real :: cf0! hydraulic coefficient, in km/yr
  real :: dlim ! distance between draining streams, in km
  real, dimension(1000) :: y1,y2 !water table heights, in km
  real, dimension(1000) :: z  ! mineral ground heights, in km
  real, dimension(1000) :: p  ! peat layer thickness, in km
  real, dimension(1000) :: s  ! peat surface heights, in km
  
  nyrs=100 ! the numbers of years (time horizon)
  str='empty'
  narg=command_argument_count()
  
  do k=1,narg
   call get_command_argument(k,str)
   !Process string str to determine what this argument is requesting
   if (str=='nn500') nyrs=500
   if (str=='n1000') nyrs=1000
   write(*,*) str
  end do
  
  open (20, file='idmout.dat')
  n=100    ! the number of x steps (spatial resolution)
  tau=0.001 ! the length of time step
  dlim=1.0 ! distance between draining streams, in km
  cf0=0.35 ! hydraulic coefficient, in km/yr
  u0=35.0  ! effective rainfall, in mm of water layer per year
  xstep=1.0/real(n) ! the length of x step in km
  m=nyrs*nint(1.0/tau) ! the number of time steps
  U=(dlim**2)*u0/(10**6) !the equivalent effective rainfall, in km
  
! set up mineral ground surface, initial
! water table surface, peat thickness, and peat surface  
  call initial (n,y1,y2,z,p,s) 
  write(*,*) s(n/2) ! peat surface height in the center, in km

!  determine the shape of water table surface and peat surface
!  after the period given by 'nyrs'
  
  do k=1,m
    call nextstep (tau,xstep,n,y1,y2,z,p,s,U,cf0)
  end do
  
! write down the table of water table and peat surface heights
! at each x step, in meters
   
  do i=1,n
    write(20,'(i4,f15.2, f8.2)') i,y2(i)*10**3, s(i)*10**3
  end do
  
! give a message that calculations was done 
  call first(hi)
  write (*,*) hi

! definitions for subroutines used in calculations  
  contains
  
   subroutine first(hi)
   ! the message appeared on screen in the end of calculations
     character (len=3):: hi
     hi='Hi!'
   end subroutine first
  
   subroutine initial(n,y1,y2,z,p,s)
   ! this is to set up mineral ground surface, initial
   ! water table surface, peat thickness, and peat surface 
     integer :: n, i, j
     real, dimension(1000) :: y1,y2 !water table heights
     real, dimension(1000) :: z  ! mineral ground heights
     real, dimension(1000) :: p  ! peat layer thickness
     real, dimension(1000) :: s  ! peat heights
	 real, parameter :: pi = 3.1415927
	 do i=1,n+1
	  y1(i)=0
	  y2(i)=0
	  p(i)=0
	  z(i)=(4.5+1.0*sin(3.0*pi*(i-1)/n))/10**3
	  s(i)=z(i)+p(i)
	 end do
   end subroutine initial
   
   subroutine nextstep(tau,xstep,n,y1,y2,z,p,s,U,cf0)
     ! this is to determine the shape of water table surface and peat surface
     !  after the period given by 'nyrs'
     integer :: n !number of x nodes
	 real :: tau ! time step
	 real :: xstep ! x step
	 real :: cf0! hydraulic coefficient
     real :: U ! effective rainfall
	 real, dimension(1000) :: y1,y2 !water table heights
     real, dimension(1000) :: z  ! mineral ground heights
     real, dimension(1000) :: p  ! peat layer thickness
     real, dimension(1000) :: s  ! peat heights
	 ! local variables
	 integer :: i, j
	 real, dimension(1000) :: A,B,C,F,alpha, beta, q
	 q(1)=0.5*Cf(cf0,y1(1),s(1))
	 do i=2,n
	   q(i)=0.5*(Cf(cf0,y1(i-1),s(i-1))+Cf(cf0,y1(i),s(i)))
	 end do
	
	 do i=1,n-1
	   A(i)=-(tau/(xstep*xstep))*q(i)
	   B(i)=-(tau/(xstep*xstep))*q(i+1)
	   C(i)=-(1+(tau/(xstep*xstep))*(q(i+1)+q(i)))
	   F(i)=-(y1(i)+tau*U)  
	 end do
	 
	 alpha(1)=0
	 beta(1)=0
	 
	 do i=1,n-1
	   alpha(i+1)=B(i)/(C(i)-alpha(i)*A(i))
	   beta(i+1)=(A(i)*beta(i)+F(i))/(C(i)-alpha(i)*A(i))
	 end do
	 
	 y2(n)=0
	 
	 do i=1,n-1
	  j=n-i  ! j from n-1 to 1; y2(0)=0
	  y2(j)=alpha(j+1)*y2(j+1)+beta(j+1)
	 end do
	 
	 do i=1,n
	  y1(i)=y2(i)
	  s(i)=s(i)+tau*ds(y2(i),s(i))
	 end do
	
   end subroutine nextstep
   
   function Cf(cf0,y,s)
    ! this is to define the hydraulic conductivity as a function
	! of water table
    real :: y,s,cf0, b, d
	b=100000.0 ! km -> cm
	d=15.0 !cm
	Cf=min(500000.0,cf0*(y+exp(b*(y-s)+d)/b-exp(d-b*s)/b))
	!Cf=cf0*y !assuming only subsurface runoff 
   end function Cf
   
   function ds(y,s)
    ! this is to define the rate of peat growth as a function 
	! of water table
     real :: y,s,peatinc,ubv
	 ubv=(s-y)*10**5
	 peatinc=(10.0/10**6)*(exp(-(1.0-ubv/20.0)**2)-exp(-1.0))
	 if ((ubv < 0.0).or.(ubv > 40.0)) peatinc=0.0
	 ds=peatinc
   end function ds
   
end program main
