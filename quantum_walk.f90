program qrw_finite_lattice

use parametros 
use funcoes

implicit none
	
!----------------------------------------------------------------------------------!
!Variáveis                                                                         !
!----------------------------------------------------------------------------------!	

 complex*16,dimension(1:3,-1100:1100,-1100:1100)::ket0,ket
 real*8 :: ptotal,prob,probx,proby,pbulk,pborda,pxborda,pyborda,pr
 real*8 :: rmed,r2med,xmed,x2med,ymed,y2med
 integer:: i,j,k,s,w, jlim, klim,m,jj,kk
 character (len=6):: nome

ptotal=0.0d0
prob=0.0d0
probx=-0.0d0
proby=0.0d0
rmed=0.0d0
r2med=0.0d0
xmed=0.0d0
x2med=0.0d0
ymed=0.0d0
y2med=0.0d0
ket0=(0.0d0,0.0d0)
ket=(0.0d0,0.0d0)
pbulk=0.0d0
pxborda=0.0d0
pyborda=0.0d0
pborda=0.0d0
pr=0.0d0
!---------------------------------------------!
!Estado inicial do sistema.
! ket0(3,0,0)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(2,0,1)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(1,1,2)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(3,2,1)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(2,2,0)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(1,1,-1)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(2,0,0)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(3,1,-1)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(1,2,0)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(2,2,1)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(3,1,2)=(1.0d0,0.0d0)/sqrt(12.0d0)
! ket0(1,0,1)=(1.0d0,0.0d0)/sqrt(12.0d0)


 ket0(2,0,0)=(1.0d0,0.0d0)/sqrt(2.0d0)
 ket0(2,0,-3)=(1.0d0,0.0d0)/sqrt(2.0d0)
!---------------------------------------------!




do w=0,step 

nome=char(48+mode(int(w/100000),10))//char(48+mode(int(w/10000),10))//char(48+mode(int(w/1000),10))//&
    &char(48+mode(int(w/100),10))//char(48+mode(int(w/10),10))//char(48+mode(w,10))


 if(w<jmax .and. w<kmax)then
  jlim=w+5
  klim=w+5
 elseif(w<jmax .and. w>kmax)then
  jlim=w+5
  klim=kmax+5
 elseif(w>jmax .and. w<kmax)then
  jlim=jmax+5
  klim=w+5
 else
  jlim=jmax+5
  klim=kmax+5
 endif
 !------------------------------
 !Evolução temporal de um passo.
 if(w>0)then
  do j=-jlim,jlim
   do k=-klim,klim
    do s=1,3
     ket(1,f(1,j,k),g(1,j,k))=ket(1,f(1,j,k),g(1,j,k))+coef(1,s,j,k)*ket0(s,j,k)
     ket(2,f(2,j,k),g(2,j,k))=ket(2,f(2,j,k),g(2,j,k))+coef(2,s,j,k)*ket0(s,j,k)
     ket(3,f(3,j,k),g(3,j,k))=ket(3,f(3,j,k),g(3,j,k))+coef(3,s,j,k)*ket0(s,j,k)
    enddo
   enddo
  enddo
 else
  ket=ket0
 endif 
 

 do j=jmin,jmax
  do k=kmin,kmax  
   if(((mode(k,4)==0.or.mode(k,4)==1).and.mode(j,2)==0).or.&
    &((mode(k,4)==2.or.mode(k,4)==3).and.mode(j,2)==1))then
    !-------------------------------------------------------------------------------------------------!
    !Condição de normalização
    ptotal=ptotal+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2 
    !-------------------------------------------------------------------------------------------------!
    !Distribuição espacial de probabilidades de encontrar a partícula convergindo para o vértice (j,k)
!     if(w==Nz.or.w==2*Nz.or.w==3*Nz.or.w==4*Nz.or.w==5*Nz)then
      open(20,file="Prob"//nome//".dat")
      write(20,*)j,k,abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2
!     endif
    endif
    
    if(j==jmin.or.j==jmin+1.or.j==jmin+2.or.j==jmax.or.j==jmax-1.or.j==jmax-2)then
     pyborda=pyborda+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2
    endif
    if(k==kmin.or.k==kmin+1.or.k==kmin+2.or.k==kmin+3.or.k==kmax.or.k==kmax-1.or.k==kmax-2.or.k==kmax-3)then
     pxborda=pxborda+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2
    endif
    if(j==jmin.or.j==jmin+1.or.j==jmin+2.or.j==jmax.or.j==jmax-1.or.j==jmax-2 .or. & 
    & k==kmin.or.k==kmin+1.or.k==kmin+2.or.k==kmin+3.or.k==kmax.or.k==kmax-1.or.k==kmax-2.or.k==kmax-3)then
     pborda=pborda+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2
    else
     pbulk=pbulk+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2
    endif
   
   !-------------------------------------------------------------------------------------------------!
   !<r>
   rmed = rmed + (abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2)*sqrt(Xv(j,k)**2+Yv(j,k)**2)
   !<r**2>
   r2med = r2med + (abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2)*(Xv(j,k)**2+Yv(j,k)**2)
   !<x>
   xmed = xmed + (abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2)*Xv(j,k)
   !<x**2>
   x2med = x2med + (abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2)*Xv(j,k)**2
   !<y>
   ymed = ymed + (abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2)*Yv(j,k)
   !<y**2>
   y2med = y2med + (abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2)*Yv(j,k)**2
  
   probx=probx+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2

  enddo
!  if(w==Nz.or.w==2*Nz.or.w==3*Nz.or.w==4*Nz.or.w==5*Nz)then
   open(33,file="Probx"//nome//".dat")
   write(33,*)j,probx
!  endif
  probx=0.0d0
 
 enddo 

 do k=kmin,kmax
  do j=jmin,jmax  
   proby=proby+abs(ket(1,j,k))**2+abs(ket(2,j,k))**2+abs(ket(3,j,k))**2
  enddo
!  if(w==Nz.or.w==2*Nz.or.w==3*Nz.or.w==4*Nz.or.w==5*Nz)then
   open(37,file="Proby"//nome//".dat")
   write(37,*)k,proby
!  endif
  proby=0.0d0
 enddo
  
 if(w==Nz.or.w==2*Nz.or.w==3*Nz.or.w==4*Nz.or.w==5*Nz)then
  do m=0,(jmax-1)/2
   do jj=-(2*m+1),(2*m+1)
    do kk=-(2*m+3),2*m 
     if(jj==-2*m-1.or.jj==-2*m.or.jj==2*m.or.jj==2*m+1.or.kk==-2*m-3.or.kk==-2*m-2.or.kk==2*m-1.or.kk==2*m)then
      pr=pr+abs(ket(1,jj,kk))**2+abs(ket(2,jj,kk))**2+abs(ket(3,jj,kk))**2
     endif
    enddo
   enddo
   open(39,file="Probr"//nome//".dat")
   write(39,*)m,pr
   pr=0.0d0
  enddo
 endif

 open(30,file="msdR.dat")
 write(30,*)w,abs(r2med-rmed**2.0)

 open(31,file="msdX.dat")
 write(31,*)w,abs(x2med-xmed**2.0)

 open(32,file="msdY.dat")
 write(32,*)w,abs(y2med-ymed**2.0)

 open(34,file="pxborda.dat")
 write(34,*)w,pxborda

 open(35,file="pyborda.dat")
 write(35,*)w,pyborda

 open(36,file="pborda.dat")
 write(36,*)w,pborda
 
 open(38,file="pbulk.dat")
 write(38,*)w,pbulk


 print*, w,ptotal
 ptotal=0.0d0

 proby=0.0d0 
 ket0=ket
 ket=(0.0d0,0.0d0)
 rmed=0.0d0
 r2med=0.0d0
 xmed=0.0d0
 x2med=0.0d0
 ymed=0.0d0
 y2med=0.0d0
 pbulk=0.0d0
 pxborda=0.0d0
 pyborda=0.0d0
 pborda=0.0d0
enddo

end program qrw_finite_lattice
