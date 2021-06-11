module grid
 real*8 hxi,hyi,xl,yl
 integer nxp1,nyp1,nx,ny,nx1,ny1
end module

module bounds
 integer ibdry,jbdry,nxlast,nylast
end module

module front
 integer  np,ffp,lfp,fep,ne,ffe,lfe,fee
 real*8 amax,amin
end module

module press
 integer itmax
 real*8 xerr, beta
end module

module flprop
 real*8 r1,r2,ifi1,ifi2,gx,gy,rro
end module

module fparam
 real*8 xmv,mv,fxl,fyl
end module 

module source
 integer num_srce,isource_pt(20)
 real*8 fsource(20)
end module   
