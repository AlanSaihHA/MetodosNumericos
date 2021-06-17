program source
implicit none
real pi,a,u0,u,rho,Area, Ct
real r, dx, dv, epsi, xcc, xc
real k, f, m
integer i
pi = acos(-1.0)
!a = (u0-u)/u0 
Ct = 4.*a*(1.-a)
!f = ((0.5*rho*Area*u0*0*Ct*r)/(rho*u0*u0)*Area*dx))*dv
!epsi=1.*dx
!d = xcc-xc
!k = ((1./(epsi*(pi**(0.5))))*exp(-(d/epsi)**2.))/((1./(epsi*(pi**(0.5))))*exp(-(0./epsi)**2.))
epsi = 7./150.
f = 540.

do i =0, 100
m = float(i)*0.01 
k = ((1./(epsi*(pi**(0.5))))*exp(-(m/epsi)**2.))/((1./(epsi*(pi**(0.5))))*exp(-(0./epsi)**2.))
write(*,*) i,f , k , f*k
enddo

endprogram

