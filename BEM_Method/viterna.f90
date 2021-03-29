!function [AlphaTot,CLTot,CDTot,last_err]=ViternaMethod(filename,wtperfilename) 
![alpha,CL,CD,last_err] = Read_size_polares(filename);
!!$if last_err==1;
!!$    AlphaTot(:,1)=zeros(180,1);
!!$    CLTot(:,1)=AlphaTot;
!!$    CDTot(:,1)=AlphaTot;
!!$    CnTot(:,1)=AlphaTot;%%%%%NOTA: AGREGAR    
!!$else
program Viterna
  implicit none
  integer::ind(1),ind1(1)
  Real:: alpha(9)=(/-8.1,-4.2,-2.6,-0.5,1.5,3.6,5.1,7.1,9.2/)
  Real:: CL(9)=(/-0.574,-0.416,-0.198,0.056,0.286,0.535,0.737,0.914,0.914/)
  Real:: CD(9)=(/0.0742,0.0117,0.0114,0.0109,0.0102,0.0091,0.0069,0.0088,0.0525/)
  Real:: Cn(size(CL))=0.0
  Real, allocatable::CD1(:),CL1(:),CD2(:),CL2(:),CD3(:),CL3(:),CD4(:),CL4(:),Cn1(:),Cn2(:),Cn3(:),Cn4(:),alpha1(:),alpha2(:),alpha3(:),alpha4(:)
  Real, allocatable::CD5(:),CL5(:),Cn5(:),alpha5(:),CD6(:),CL6(:),Cn6(:),alpha6(:),AlphaTot(:)
  Real::CLs(1),alphas(1),CDs(1),high,low,A1(1),A2(2),B1(1),B2(1)
  Real::AR,CDmax(1),alphag,alphah,alphai,alphaf,alphaj,alphal
  integer::i,j,k
  allocate(CD1(7),CL1(7),Cn1(7),CD2(6),CL2(6),CD3(2),CL3(2),CD4(80),CL4(80),Cn2(6),Cn3(2),Cn4(80),alpha1(7),alpha2(6),alpha3(2),alpha4(80))
  allocate(CD5(60),CL5(60),Cn5(60),alpha5(60),CD6(50),CL6(10),Cn6(10),alpha6(10),AlphaTot(500))
AR=10
ind=MAXLOC(alpha) !find the position of the maximum lift coefficient
CLs=CL(ind) !stall lift coefficient
alphas=alpha(ind) !stall angle of attack
CDs=CD(ind); !stall drag coefficiente
high=maxval(alpha)
low=minval(alpha)
! write(*,*)low
ind1=minloc(alpha)

!cladj = 0.7;
open(1,file='perfil.dat',status='replace')
write(1,'(5X,a5,3X,a2,6X,a2,7X,a2)')"Alpha","CL","CD","Cn"
k=1
do alphal=-180.0,-170.0,10.0
   CDmax(1)=1.11+0.018*AR   
   B1(1)=CDmax(1)
   B2(1)=(CDs(1)-B1(1)*(sind(alphas(1)))**2.0)/cosd(alphas(1))  
   CD6(k)=B1(1)*(sind(alphal+180.0))**2+B2(1)*cosd(alphal+180.0)   
   CL6(k)= (alphal+180.0)/high*CLs(1)*0.7   
   alpha6(k)=alphal   
   Cn6(k)=0.0
     write(1,'(4F10.3)')alpha6(k),CL6(k),CD6(k),Cn6(k)
   k=k+1     
end do

i=1
do alphaj=-160.0,-100.0,10.0   
   CDmax(1)=1.11+0.018*AR
   B1=CDmax
   B2(1)=(CDs(1)-B1(1)*(sind(alphas(1)))**2.0)/cosd(alphas(1))      
   CD5(i)=B1(1)*(sind(alphaj+180.0))**2+B2(1)*cosd(alphaj+180.0)   
   !Coeficiente de sustentación
   A1(1)=B1(1)/2.0
   A2(1)=(CLs(1)-B1(1)*sind(alphas(1))*cosd(alphas(1)))*(sind(alphas(1))/(cosd(alphas(1)))**2.0)
   CL5(i)=0.7*(A1(1)*sind(2.0*(alphaj+180.0))+A2(1)*((cosd(alphaj+180.0))**2/sind(alphaj+180.0)))
   Cn5(i)=0.0
   alpha5(i)=alphaj
     write(1,'(4F10.3)')alpha5(i),CL5(k),CD5(k),Cn5(k)
!   write(*,*)alpha5(i),CL5(i),CD5(i)   
  i=i+1        
end do
  i=1
!  if (low .LE. 0.0) .OR. (low .GE. 0.0) then     
     do alphaf=-90.0,low-1.0,10.0        
        CDmax(1)=1.11+0.018*AR
        if((alphaf.EQ.0.0) .OR.( alphaf.GE.-10.0))then
           CL4(i) = -CLs(1)*0.7 + (alphaf + high)/(low + high)*(CL(ind1(1))+CLs(1)*0.7)
           CD4(i) = CD(ind1(1)) + (alphaf-low)/(-high-low)*(CDs(1)-CD(ind1(1)))
           Cn4(i)=0.0
           alpha4(i)=alphaf
           write(1,'(4F10.3)')alpha4(i),CL4(i),CD4(i),Cn4(i)
        else
           !Coeficiente de arrastre
           B1(1)=CDmax(1)
           B2(1)=(CDs(1)-B1(1)*(sind(alphas(1)))**2.0)/cosd(alphas(1))
           CD4(i)=B1(1)*(sind(alphaf))**2.0+B2(1)*cosd(alphaf)
           !Coeficiente de sustentación
           A1(1)=B1(1)/2.0
           A2(1)=(CLs(1)-B1(1)*sind(alphas(1))*cosd(alphas(1)))*(sind(alphas(1))/(cosd(alphas(1)))**2.0)
           CL4(i)=0.7*(A1(1)*sind(2.0*alphaf) + A2(1)*((cosd(alphaf))**2.0/sind(alphaf)))
           Cn4(i)=0.0
           alpha4(i)=alphaf
           write(1,'(4F10.3)')alpha4(i),CL4(i),CD4(i),Cn4(i)
        end if       
     end do    
     i=i+1
 ! end if

  do i=1,size(CL)
     write(1,'(4F10.3)')alpha(i),CL(i),CD(i),Cn(i)
  end do
  
i=1;
if (high .LT. 20.0) then
   do alphag=20.0,90.0,10.0
      CDmax(1)=1.11+0.018*AR                       !  write(*,*)CDmax
      !Coeficiente de arrastre      
      B1(1)=CDmax(1)      
      B2(1)=(CDs(1)-B1(1)*(sind(alphas(1)))**2.0)/cosd(alphas(1))      
      CD1(i)=B1(1)*(sind(alphag))**2.0+B2(1)*cosd(alphag)
      !Coeficiente de sustentación      
      A1(1)=B1(1)/2.0
      A2(1)=(CLs(1)-B1(1)*sind(alphas(1))*cosd(alphas(1)))*(sind(alphas(1))/(cosd(alphas(1)))**2.0)
      CL1(i)=A1(1)*sind(2.0*alphag)+A2(1)*((cosd(alphag))**2.0/sind(alphag))
      Cn1(i)=0;
      alpha1(i)=alphag
      write(1,'(4F10.3)')alpha1(i),CL1(i),CD1(i),Cn1(i)
      i=i+1;  
   end do
end if
i=1;
do alphah=100.0,160.0,10.0
   CDmax(1)=1.11+0.018*AR;  
     B1(1)=CDmax(1)
     B2(1)=(CDs(1)-B1(1)*(sind(alphas(1)))**.02)/cosd(alphas(1))
     CD2(i)=B1(1)*(sind(180.0-alphah))**2+B2(1)*cosd(180.0-alphah)
     !Coeficiente de sustentación
     A1(1)=B1(1)/2.0
     A2(1)=(CLs(1)-B1(1)*sind(alphas(1))*cosd(alphas(1)))*(sind(alphas(1))/(cosd(alphas(1)))**2)
     CL2(i)=-0.7*(A1(1)*sind(2.0*(180.0-alphah))+A2(1)*((cosd(180.0-alphah))**2/sind(180.0-alphah)))     
     Cn2(i)=0
     alpha2(i)=alphah
          write(1,'(4F10.3)')alpha2(i),CL2(i),CD2(i),Cn2(i)
     i=i+1;
  end do

  j=1;  
  do alphai=170.0,180.0,10.0     
     CDmax(1)=1.11+0.018*AR
     B1(1)=CDmax(1)
     B2(1)=(CDs(1)-B1(1)*(sind(alphas(1)))**2.0)/cosd(alphas(1))      
     CL3(j)=(alphai-180.0)/high*CLs(1)*0.7
     CD3(j)=B1(1)*(sind(180.0-alphai))**2+B2(1)*cosd(180.0-alphai)   
     Cn3(j)=0
     alpha3(j)=alphai
       write(1,'(4F10.3)')alpha3(j),CL3(j),CD3(j),Cn3(j)
     j=j+1
  end do
!AlphaTot(:)=[alpha6,alpha5,alpha4,alpha,alpha1,alpha2,alpha3]
!!$CLTot(:,1)=[CL6;CL5;CL4;CL;CL1;CL2;CL3];
!!$CDTot(:,1)=[CD6;CD5;CD4;CD;CD1;CD2;CD3];
!!$CnTot(:,1)=[Cn6;Cn5;Cn4;Cn;Cn1;Cn2;Cn3];%%%%%NOTA: AGREGAR
!write(*,*)AlphaTot
close(1)
end program Viterna

