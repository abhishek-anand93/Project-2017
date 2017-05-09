!===================================================================================================================================
!definition of some other functions
real*8 function e(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  e = (0.5)*cos((sqrt(3.0)*y)) + 0.75 + cos(3*x/2)*cos(sqrt(3.0)*y/2) + (0.75)*(t**2)                  !verified
end function e

real*8 function ox(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  ox = (cos(sqrt(3.0)*y/2)**2) + (0.5)*cos((3*x - sqrt(3.0)*y)/2) - (t/2)*sin((3*x + sqrt(3.0)*y)/2)   !verified
end function ox

real*8 function oy(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  oy = (cos(sqrt(3.0)*y/2)**2) + (0.5)*cos((3*x + sqrt(3.0)*y)/2) + (t/2)*sin((3*x - sqrt(3.0)*y)/2)   !verified
end function oy

real*8 function oz(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  oz = cos(3*x/2)*cos(sqrt(3.0)*y/2) + (0.5) + (t/2)*sin(sqrt(3.0)*y)                                  !verified
end function oz

!====================================================================================================================================
!polar form to general form conversion

complex*16 function com(arg)
  implicit none
  real*8, intent(in)::arg
  com = dcmplx(dcos(arg),dsin(arg))
end function com

real*8 function k1(x,y)
  implicit none
  real*8, intent(in)::x,y
  k1=(x -sqrt(3.0)*y)/2.0
end function k1

real*8 function k2(x,y)
  implicit none
  real*8, intent(in)::x,y
  k2=(x +sqrt(3.0)*y)/2.0
end function k2

real*8 function k3(x)
  implicit none
  real*8, intent(in)::x
  k3 = -x
end function k3

!======================================================================================================================================
real*8 function o_mod(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  real*8::ox,oy,oz
  o_mod = sqrt(ox(x,y,t)**2 + oy(x,y,t)**2 + oz(x,y,t)**2)           !verified
end function o_mod

real*8 function epp(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  real*8::e,o_mod
  epp = sqrt(e(x,y,t) - t*o_mod(x,y,t)) 
end function epp

real*8 function theta(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  real*8::o_mod,oz
  theta = acos(oz(x,y,t)/o_mod(x,y,t))   
end function theta

real*8 function phi(x,y,t)
  implicit none
  real*8, intent(in)::x,y,t
  real*8::ox,oy,o_mod,theta
  phi = atan(oy(x,y,t)/ox(x,y,t)) 
end function phi
!===========================================================================================================================================
!main inner product function
subroutine eigen(x,y,t,state)
   implicit none
   real*8, intent(in)::x,y,t
   real*8::theta,phi,k1,k2,k3,epp
   complex*16::state(4,1),com
   complex*16, parameter::iota = (0.0,1.0d+0)
   state(1,1) = cos(theta(x,y,t)/2.0)
   state(2,1) = exp(iota*phi(x,y,t))*sin(theta(x,y,t)/2.0d0)
   state(3,1) = 0.5d0*(com(-k1(x,y)) + com(-k2(x,y)) + com(-k3(x)))*cos(theta(x,y,t)/2) &
              + 0.5d0*t*(com(phi(x,y,t)-k1(x,y))*sin(theta(x,y,t)/2) &
              - iota*com(phi(x,y,t)-k2(x,y))*sin(theta(x,y,t)/2) &
              + com(-k3(x))*cos(theta(x,y,t)/2))
   state(4,1) = 0.5d0*(com(-k1(x,y)) + com(-k2(x,y)) + com(-k3(x)))*sin(theta(x,y,t)/2)*com(phi(x,y,t)) &
              + 0.5d0*t*(com(-k1(x,y))*cos(theta(x,y,t)/2) &
              + iota*com(-k2(x,y))*cos(theta(x,y,t)/2) &
              - com(phi(x,y,t)-k3(x))*sin(theta(x,y,t)/2))
   state = (1/sqrt(2.0d0))*state
end subroutine eigen

 !================================================================================================================================================
!derivative
subroutine hamil(x,y,t,h)
  implicit none
  real*8, intent(in)::x,y,t
  real*8::k1,k2,k3
  complex*16, parameter::iota=(0.0d+0,1.0d+0)
  integer::i,j

  complex*16::com,h(4,4)

  
  h = (0.0d+0,0.0d+0)

  do i=1,2
  h(i,i) = h(i,i) + 0.01d+0
  h(i+2,i+2) = h(i+2,i+2) - 0.01d+0
  enddo

  h(1,3)= 0.5*(com(-k1(x,y))+com(-k2(x,y))+com(-k3(x))) + (t/2)*com(-k3(x))
  h(1,4)= (t/2)*(com(-k1(x,y)) - iota*com(-k2(x,y)))
  h(2,3)= (t/2)*(com(-k1(x,y)) + iota*com(-k2(x,y)))
  h(2,4)= 0.5*(com(-k1(x,y))+com(-k2(x,y))+com(-k3(x))) - (t/2)*com(-k3(x))
  h(3,1)= 0.5*(com(k1(x,y))+com(k2(x,y))+com(k3(x))) + (t/2)*com(k3(x))
  h(3,2)= (t/2)*(com(k1(x,y)) - iota*com(k2(x,y)))
  h(4,1)= (t/2)*(com(k1(x,y)) + iota*com(k2(x,y)))
  h(4,2)= 0.5*(com(k1(x,y))+com(-k2(x,y))+com(-k3(x))) - (t/2)*com(k3(x))
end subroutine hamil
!===================================================================================================================================================
function innerprod(u,v)
implicit none
complex*16::u(4),v(4),p,innerprod
integer::i
p=(0.d0,0.d0)
do i=1,4
p=p+conjg(u(i))*v(i)
enddo
innerprod=p
end function

!====================================================================================================================================================

program chern
  implicit none
  real*8,parameter::pi=acos(-1.d0)
  real*8, parameter::x_range=4.d0*pi/3.d0,y_range=2.d0*pi/sqrt(3.d0)
  integer*8::nx,ny,i,j,ei
  real*8::x,y,t,hx,hy,v1,v2,sum2 =0.0d+0,epp,k1,k2,k3,phi,hv1,hv2
  complex*16, parameter::iota=(0.0,1.0) 
  complex*16::a_x,a_y,sum1=(0.0d+0,0.0d+0),inner,d_x,test,inner1,com
  integer,parameter::n=20
  real*8::kx(n*n),ky(n*n),en(4),en2(2)
  complex*16::state_k(n*n,4),st_k_derx(n*n,4),st_k_dery(n*n,4),innerprod,fukui(n*n)
  complex*16::state_c(1,4),state(4,1),h(4,4),eigv(4,4),c1,c2,c3,c4

  !t=0.0001d0
  Print*,'Enter the t_dash value :'
  read*,t
  print*,'Enter the eigenstate (1 to 4) :'
500  read*,ei
  if (ei < 1 .OR. ei>4) then
     print*,'Wrong input. Number should be between 0-4! Enter again :'
     goto 500
  endif
  
  open(10,file='curvature.dat')
  open(20,file='guage.dat')
  hv1 = (x_range/n)
  hv2 = (y_range/n)
  hx = hv1 + (1.0d0/sqrt(3.0d0))*hv2
  hy = hv2
  do i = 1,n
     do j=1,n
        v1 = (i-0.5)*hv1
        v2 = (j-0.5)*hv2
        x = (v1 + (1/sqrt(3.0))*v2)
        y = v2
        kx((i-1)*n+j)=x
	ky((i-1)*n+j)=y
        call hamil(x,y,t,h(:,:))
	call Diagon_Alley(h,4,en,eigv)
        state_k((i-1)*n+j,:) = eigv(:,ei)
     enddo
  enddo


fukui = (0.0d0,0.0d0)
do i=1,n
do j=1,n
 if (i==n .and. j/=n) then
 c1 = innerprod(state_k((i-1)*n+j,:),state_k((0)*n+j,:))/abs(innerprod(state_k((i-1)*n+j,:),state_k((0)*n+j,:)))
 c2 = innerprod(state_k((0)*n+j,:),state_k((0)*n+j+1,:))/abs(innerprod(state_k((0)*n+j,:),state_k((0)*n+j+1,:)))
 c3 = innerprod(state_k((0)*n+j+1,:),state_k((i-1)*n+j+1,:))/abs(innerprod(state_k((0)*n+j+1,:),state_k((i-1)*n+j+1,:)))
 c4 = innerprod(state_k((i-1)*n+j+1,:),state_k((i-1)*n+j,:))/abs(innerprod(state_k((i-1)*n+j+1,:),state_k((i-1)*n+j,:)))
 fukui((i-1)*n+j) = log(c1*c2*c3*c4)
 elseif (j==n .and. i/=n) then
 c1 = innerprod(state_k((i-1)*n+j,:),state_k((i)*n+j,:))/abs(innerprod(state_k((i-1)*n+j,:),state_k((i)*n+j,:)))
 c2 = innerprod(state_k((i)*n+j,:),state_k((i)*n+1,:))/abs(innerprod(state_k((i)*n+j,:),state_k((i)*n+1,:)))
 c3 = innerprod(state_k((i)*n+1,:),state_k((i-1)*n+1,:))/abs(innerprod(state_k((i)*n+1,:),state_k((i-1)*n+1,:)))
 c4 = innerprod(state_k((i-1)*n+1,:),state_k((i-1)*n+j,:))/abs(innerprod(state_k((i-1)*n+1,:),state_k((i-1)*n+j,:)))
 fukui((i-1)*n+j) = log(c1*c2*c3*c4)
 elseif (i==n .and. j==n) then
 c1 = innerprod(state_k((i-1)*n+j,:),state_k((0)*n+j,:))/abs(innerprod(state_k((i-1)*n+j,:),state_k((0)*n+j,:)))
 c2 = innerprod(state_k((0)*n+j,:),state_k((0)*n+1,:))/abs(innerprod(state_k((0)*n+j,:),state_k((0)*n+1,:)))
 c3 = innerprod(state_k((0)*n+1,:),state_k((i-1)*n+1,:))/abs(innerprod(state_k((0)*n+1,:),state_k((i-1)*n+1,:)))
 c4 = innerprod(state_k((i-1)*n+1,:),state_k((i-1)*n+j,:))/abs(innerprod(state_k((i-1)*n+1,:),state_k((i-1)*n+j,:)))
 fukui((i-1)*n+j) = log(c1*c2*c3*c4)
 else
 c1 = innerprod(state_k((i-1)*n+j,:),state_k((i)*n+j,:))/abs(innerprod(state_k((i-1)*n+j,:),state_k((i)*n+j,:)))
 c2 = innerprod(state_k((i)*n+j,:),state_k((i)*n+j+1,:))/abs(innerprod(state_k((i)*n+j,:),state_k((i)*n+j+1,:)))
 c3 = innerprod(state_k((i)*n+j+1,:),state_k((i-1)*n+j+1,:))/abs(innerprod(state_k((i)*n+j+1,:),state_k((i-1)*n+j+1,:)))
 c4 = innerprod(state_k((i-1)*n+j+1,:),state_k((i-1)*n+j,:))/abs(innerprod(state_k((i-1)*n+j+1,:),state_k((i-1)*n+j,:)))
 fukui((i-1)*n+j) = log(c1*c2*c3*c4)
endif
enddo
enddo


st_k_derx=(0.d0,0.d0)
st_k_dery=(0.d0,0.d0)


do i=2,n-1
do j=2,n-1
st_k_derx((i-1)*n+j,:)=(state_k((i)*n+j,:)-state_k((i-2)*n+j,:))/(2.0d0*hx)
enddo
enddo


do i=2,n-1
do j=2,n-1
st_k_dery((i-1)*n+j,:)=(state_k((i-1)*n+j+1,:)-state_k((i-1)*n+j-1,:))/(2.0d0*hy)
enddo
enddo



sum1=(0.d0,0.d0)
do i=1,n*n

sum2= sum2 + aimag(fukui(i))
write(10,*)(kx(i) - (1.0/sqrt(3.0))*ky(i)),ky(i),aimag(fukui(i))
enddo
print*,'Chern number is',sum2/(2.0d0*pi)
end program chern


	subroutine Diagon_Alley(R,N,W,Mat)
	integer::N,info,Iwork(3+5*N)
	complex*16:: R(N,N)
	Complex*16::Mat(N,N),Work(2*N+N**2)
	real*8::W(N),RWORK(1+5*N+2*N**2)
	Mat=(0.0d0,0.0d0)
	Mat=R

	W=0.0d0
 

	call zheevd ('V','U',N,Mat,N,W,WORK,2*N+N**2,RWORK,1+5*N+2*N**2,IWORK,3+5*N, INFO)


	end subroutine
