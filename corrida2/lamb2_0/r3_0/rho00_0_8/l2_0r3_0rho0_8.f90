Module MathConstants
 Real(8), Parameter :: pi    = 3.14159265358979324
 Real(8), Parameter :: dospi = 2.0_8*pi
 Real(8), Parameter :: pi2   = pi/2.0_8
End Module MathConstants

Module Parameters
 Use MathConstants

 !*Namelist*!
 Namelist /input/ kmax,kt,kf,ks,check
 !*input namelist*!
 Integer            :: kmax,kt,kf,ks,check 
 !*Namelist2*!
! Namelist /input2/ Rcyl,rho
! !*input namelist2*!
! Real(8), Parameter ::Rcyl,rho

 Real(8),parameter  :: Rcyl   =3.0
 Real(8), Parameter :: rho    =0.8
 Real(8), Parameter :: xA     = 0.0
 Real(8), Parameter :: dT     = 0.01
 Real(8), Parameter :: To     = 1.0
 Real(8), Parameter :: T1     = 0.1
 Real(8), Parameter :: sigma  = 1.0
 Real(8), Parameter :: epsilon= 1.0
 Real(8), Parameter :: lambdaA= 2.0
 Real(8), Parameter :: lambdaB= 2.0
 Real(8), Parameter :: lambda_max = max(lambdaA,lambdaB)
 Real(8), Parameter :: Rcut  = 1.5*sigma*lambda_max
 Real(8), Parameter :: Lx     = 2.0*Rcyl + 2.0*Rcut 
 Real(8), Parameter :: Ly     = 2.0*Rcyl + 2.0*Rcut 
 Real(8), Parameter :: Lz     = 40.0_8      
 Character(len=20), Parameter :: efile   = 'energy18.dat'    !file name
 Character(len=20), Parameter :: pfile   = 'config'    !file name
 Character(len=20), Parameter :: ifile   = 'data18.dat'      !file name
 Character(len=20), Parameter :: sfile   = 'l2_0r3_0rho0_8.xyz'    !file name
 Character(len=20), Parameter :: finput  = 'input.in'      !file name 
 Integer, Parameter :: nx    = int(Lx/Rcut)            !cell-list
 Integer, Parameter :: ny    = int(Ly/Rcut)            !cell-list
 Integer, Parameter :: nz    = int(Lz/Rcut)            !cell-list 
 Real(8), Parameter :: dLx   = Lx/float(nx)
 Real(8), Parameter :: dLy   = Ly/float(ny)
 Real(8), Parameter :: dLz   = Lz/float(nz) 
 Integer, Parameter :: ntcell  = nx*ny*nz
 Integer, Parameter :: mapsize = 13*ntcell
 Real(8), Parameter :: Vol     = pi*Rcyl*Rcyl*Lz  !Lx*Ly*Lz
 Integer, Parameter :: Np      = int(rho*Vol)
 Integer, Parameter :: NB      = int(Np/(1.0 + xA))
 Integer, Parameter :: NA      = Np - NB
 Real(8), Parameter :: dmax    = 0.2
End Module Parameters


Module Particulas

 Type Particles
  Real(8) :: x
  Real(8) :: y
  Real(8) :: z
  Real(8) :: lambda 
  Integer :: tipo
 End Type Particles 


End Module Particulas


Module Globals
 Use Parameters
 Use Particulas

 Type (Particles) :: part(1:Np)
 Real(8)          :: E,Temp
 Integer          :: list(1:Np),head(1:ntcell)
 Integer          :: cfile,ki
 Real(8)          :: inicio,fin,cont
 Integer          :: map(1:mapsize)
End Module Globals

!#########################################    Main Program    #########################################!
Program square
 Use Parameters
 Use Globals
 Implicit None

 Integer :: i,k

 Call SetUp 
 Call Initialization
 Call Neighbor_Cells
 Call Linked_Lists
 Call Hnb

 Call cpu_time(inicio)
 Do k = ki+1,kmax
   Do i = 1,Np
    Call MonteCarlo
   End Do

  ! if(mod(k,kf) == 0) then
  !   Call Save_xyz
  ! end if

   if(mod(k,ks) == 0) then
     Write(11,FMT='(3(F20.10))')float(k),E/float(Np),Temp
   end if

   if(mod(k,check) == 0) then
     Call Save_System(k)
   end if

   if(mod(k,kt) == 0) then
      Temp = Temp - dT
      if(Temp .lt. T1) then
        Temp = T1
      end if
   end if

 End Do
 Call cpu_time(fin)
 Call SimulationData
 !k=0
 Call Save_System(k)
 Call Save_xyz

End Program square
!####################   End Main Program  ########################################!


Subroutine SetUp
 Use Parameters
 Use Globals
 Implicit None

!***Read input file***!
 Open(8,FILE=finput,STATUS='Old')
 Read(8,nml = input)
 Close(8)
!*********************!

End Subroutine SetUp


Subroutine Initialization
  Use Parameters
  Use Globals 
  Implicit None
  
  Logical :: lexist
  Integer :: i,j,k

  Call RandomSeed

  Inquire(FILE=pfile,Exist=lexist)

  if(lexist) then
   Open(10,FILE=pfile,STATUS='Old',Action='Read')
   Read(10,*)ki
   Read(10,*)cfile 
   Read(10,*)Temp
   Do i=1,Np
     Read(10,*)part(i)%x,part(i)%y,part(i)%z,part(i)%lambda,part(i)%tipo
   End Do
   Close(10)
   if(ki .ne. 0) then
     Open(11,FILE=efile,STATUS='Old',Position='append',Action='Write')
   else  
     Open(11,FILE=efile,STATUS='replace')
     Call Save_xyz 
     Temp = To
   end if
 else
   Open(11,FILE=efile,STATUS='replace')
   Call Initial
   Call Types
   ki     = 0
   cfile  = 0
   Temp   = To
   Call Save_xyz
 end if

End Subroutine Initialization


Subroutine RandomSeed
 Implicit None

 Integer, Dimension(:), Allocatable :: seed_random
 Integer                            :: i,seed

 Open(8,File='random_seed.in',status='unknown')
  
 Call Random_Seed
 Call Random_Seed(SIZE = i)
 Allocate(seed_random(i))
 Call Random_Seed(GET = seed_random) 
 Read(8,*)seed
 seed_random = ( seed_random + seed )/2
 Call Random_Seed(PUT = seed_random)   
 Deallocate(seed_random)
 Close(8)

End Subroutine RandomSeed


Subroutine Save_xyz
  Use Parameters
  Use Globals
  Implicit None

  Integer          :: i
  Character(len=6) :: namefile

  write(namefile,'(I6)')cfile
  Open(15+cfile,FILE=sfile,STATUS='replace')
  Write(15+cfile,*) Np
  Write(15+cfile,*)
  Do i = 1,Np
   if(part(i)%tipo == -1) then
     Write(15+cfile,*)"F   ", part(i)%x,part(i)%y,part(i)%z
   else
     Write(15+cfile,*)"H   ", part(i)%x,part(i)%y,part(i)%z
   end if
  End Do       
  Close(15+cfile)  
  cfile=cfile+1

End Subroutine Save_xyz


Subroutine Initial
 Use Parameters
 Use Globals
 Implicit None

 Real(8) :: p,rr,rxp,ryp
 Integer :: i,j,flag  
                                                    
 
 Do i = 1,Np
   flag = 1
   Do While(flag == 1)  
     Call Random_Number(p)
     part(i)%x = p*Lx
     Call Random_Number(p)
     part(i)%y = p*Ly
     rxp = part(i)%x - Lx/2
     ryp = part(i)%y - Ly/2
     rr = rxp*rxp + ryp*ryp
     if(rr .lt. (Rcyl*Rcyl)) then
       flag = 0
     end if
   End Do
 
   Call Random_Number(p)
   part(i)%z = p*Lz
 End Do

End Subroutine Initial


Subroutine Types
 Use Parameters
 Use Globals
 Implicit None
 
 Integer :: i,j,cont1

 cont1 = 0
 Do i = 1,Np
  cont1 = cont1 + 1
  if(cont1 <= NA) then
    part(i)%tipo   = -1
    part(i)%lambda = lambdaA
  else
    part(i)%tipo = +1
    part(i)%lambda = lambdaB
  end if
 End Do

End Subroutine Types


Subroutine Neighbor_Cells
 Use Parameters
 Use Globals
 Implicit None

 Integer :: ix,iy,iz,imap,icell

 Do iz = 1,nz
  Do iy = 1,ny
   Do ix = 1,nx
     Call Ixcell(ix,iy,iz,icell)
     imap = ( icell - 1 ) * 13
     Call Ixcell(ix+1,iy,iz,icell)
     map(imap + 1) = icell
     Call Ixcell(ix+1,iy+1,iz,icell)
     map(imap + 2) = icell
     Call Ixcell(ix,iy+1,iz,icell) 
     map(imap + 3) = icell
     Call Ixcell(ix-1,iy+1,iz,icell) 
     map(imap + 4) = icell
     Call Ixcell(ix+1,iy,iz-1,icell)
     map(imap + 5) = icell
     Call Ixcell(ix+1,iy+1,iz-1,icell) 
     map(imap + 6) = icell
     Call Ixcell(ix,iy+1,iz-1,icell)
     map(imap + 7) = icell
     Call Ixcell(ix-1,iy+1,iz-1,icell)
     map(imap + 8) = icell
     Call Ixcell(ix+1,iy,iz+1,icell)
     map(imap + 9) = icell
     Call Ixcell(ix+1,iy+1,iz+1,icell)
     map(imap + 10) = icell
     Call Ixcell(ix,iy+1,iz+1,icell)
     map(imap + 11) = icell
     Call Ixcell(ix-1,iy+1,iz+1,icell) 
     map(imap + 12) = icell
     Call Ixcell(ix,iy,iz+1,icell) 
     map(imap + 13) = icell
   End Do
  End Do
 End Do

Contains

Subroutine Ixcell(ix,iy,iz,icell)
 Use MathConstants   
 Implicit None
  
 Integer, Intent(In)  :: ix
 Integer, Intent(In)  :: iy 
 Integer, Intent(In)  :: iz 
 Integer, Intent(Out) :: icell
 
 icell = 1 + mod(ix-1+nx,nx) + mod(iy-1+ny,ny)*nx + mod(iz-1+nz,nz)*nx*ny 

End Subroutine Ixcell

End Subroutine Neighbor_Cells


Subroutine Linked_Lists
 Use Parameters
 Use Globals
 Implicit None

 Integer :: i,k,l,q,icell

 head = 0

 Do i = 1,Np
   k = int(part(i)%x/dLx)
   l = int(part(i)%y/dLy)
   q = int(part(i)%z/dLz)

   icell = 1 + k + l*nx + q*nx*ny

   list(i)     = head(icell)
   head(icell) = i
 End Do

End Subroutine Linked_Lists


Subroutine Save_System(iter)
 Use Parameters
 Use Globals
 Implicit None

 Integer, Intent(In) :: iter

 Integer :: i,j,k

 if(iter == 0) then
   Open(10,FILE=pfile,STATUS='replace')
   Write(10,*) '0'
   Write(10,*) '0'
   Write(10,*) '0'
   Do i=1,Np
     Write(10,*) part(i)%x,part(i)%y,part(i)%z,part(i)%lambda,part(i)%tipo
   End Do
   Close(10)
 else 
  if(iter == (kmax+1)) then 

   Open(10,FILE=pfile,STATUS='replace')
   Write(10,*) iter-1
   Write(10,*) cfile
   Write(10,*) Temp
   Do i=1,Np
     Write(10,*) part(i)%x,part(i)%y,part(i)%z,part(i)%lambda,part(i)%tipo
   End Do
   Close(10)
 
  else

   Open(10,FILE=pfile,STATUS='replace')
   Write(10,*) iter
   Write(10,*) cfile
   Write(10,*) Temp
   Do i=1,Np
     Write(10,*) part(i)%x,part(i)%y,part(i)%z,part(i)%lambda,part(i)%tipo
   End Do
   Close(10) 
  end if
 end if

End Subroutine Save_System


Subroutine SimulationData
 Use Parameters
 Use Globals
 Implicit None

 Real(8) :: Ebb 

 Ebb = E
 Call Hnb
 cont = cont*100.0

 Open(13,FILE=ifile,STATUS='replace') 
 
 Write(13,*)'==============================================================================='
 Write(13,*)'=         Monte Carlo Simulation of Soft-Nanoparticles In Bulk v1.0           ='
 Write(13,*)'=                                                                             =' 
 Write(13,*)'=                                                                             ='
 Write(13,*)'=                          by A. Ramirez-Hernandez                            ='
 Write(13,*)'=                                                                             =' 
 Write(13,*)'==============================================================================='
 Write(13,*)'========================   Block Copolymer Parameters  ========================'
 Write(13,100)Np,NA,NB,rho,nx,ny,nz,lambdaA,lambdaB,Rcut
 100 Format(4X,'Number of Particles           = ',4X,I6/,&
          & 4X,'Number of A Particles         = ',4X,I6/,&
          & 4X,'Number of B Particles         = ',4X,I6/,&
          & 4X,'Number Density                = ',F10.4/,& 
          & 4X,'nx                            = ',4X,I6/,& 
          & 4X,'ny                            = ',4X,I6/,&
          & 4X,'nz                            = ',4X,I6/,&
          & 4X,'LambdaA                       = ',F10.4/,&   
          & 4X,'LambdaB                       = ',F10.4/,& 
          & 4X,'Rcut                          = ',F10.4)
  Write(13,*)'==============================================================================='
  Write(13,*)'=============================   Simulation Data   ============================='
  Write(13,300)Lx,Ly,Lz,Ebb-E,Temp,kmax,kf,check,cont/(float(kmax)*float(Np))*100.0,(fin-inicio)/60.0
 300 Format(4X,'Lx (sigma)                    = ',F10.4/,&
          & 4X,'Ly (sigma)                    = ',F10.4/,&
          & 4X,'Lz (sigma)                    = ',F10.4/,&
          & 4X,'Energy error (absolute)       = ',ES10.2/,&
          & 4X,'Temp                          = ',F10.4/,&
          & 4X,'Monte Carlo steps             = ',4X,I6/,& 
          & 4X,'MC steps for snapshots        = ',4X,I6/,&
          & 4X,'MC steps for checkpoint       = ',4X,I6/,&  
          & 4X,'acceptance ratio              = ',F10.4/,&
          & 4X,'run time (min)                = ',F10.4)
  Write(13,*)'==============================================================================='
  Close(13)
 
End Subroutine SimulationData


Subroutine distance(xi,xj,L,xr)
 Implicit None

 Real(8), Intent(In) :: xi
 Real(8), Intent(In) :: xj
 Real(8), Intent(In) :: L
 Real(8), Intent(Out):: xr
  
 xr = xi - xj
 xr = xr - L*nint(xr/L)

End Subroutine distance


Subroutine Inbox(x,y,z,flag)
 Use Parameters
 Use Globals
 Implicit None

 Real(8), Intent(InOut) :: x
 Real(8), Intent(InOut) :: y
 Real(8), Intent(InOut) :: z
 Integer, Intent(InOut) :: flag

 Real(8)  :: rr,rxp,ryp

 rxp = x - Lx/2
 ryp = y - Ly/2
 rr = rxp*rxp + ryp*ryp
 if(rr .lt. (Rcyl*Rcyl)) then
   x = x
   y = y
   z = Modulo(z,Lz)
   flag = 1
 else
   flag = 0
 end if

End Subroutine Inbox


Subroutine Hnb
 Use Parameters
 Use Globals
 Implicit None

 Integer :: icell,i,j,jcell0,jcell,jj 
 Real(8) :: rxi,ryi,rzi,fxi,fyi,fzi,rxij,ryij,rzij,rij2,rij,Unb,dd,fij,dmin

 E = 0.0_8

!** Loop over all cells **!
 Do icell = 1,ntcell

   i = head(icell)
   !** Loop over all particles in the cell **!
   Do While(i .gt. 0)
     rxi = part(i)%x
     ryi = part(i)%y
     rzi = part(i)%z
     !** Loop over all the particles in the current cell **!
     j = list(i)
     Do While(j .gt. 0)
       rxij = rxi - part(j)%x
       ryij = ryi - part(j)%y 
       rzij = rzi - part(j)%z

       rxij = rxij - anint(rxij/Lx)*Lx
       ryij = ryij - anint(ryij/Ly)*Ly
       rzij = rzij - anint(rzij/Lz)*Lz 

       rij2 = rxij*rxij + ryij*ryij + rzij*rzij
       if(rij2 .lt. (Rcut*Rcut)) then
         rij = sqrt(rij2)
         Call delta(part(i)%tipo,part(j)%tipo,dd)
         dmin = (part(i)%lambda + part(j)%lambda)*sigma/2.0
         if(rij .lt. sigma) then
           Unb = 100.0
         else
           if(rij .lt. dmin) then
             Unb = dd*epsilon
           else
             Unb = 0.0
           end if
         end if
         E = E + Unb
       end if
       j = list(j)
     End Do
     !** Loop over all neighbouring cells **!
     jcell0 = 13*(icell - 1)
     Do jj = 1, 13
       jcell = map(jcell0 + jj)
       !** Loop over all the particles in neighbouring cells **!
       j = head(jcell)
       Do While(j .gt. 0)
         rxij = rxi - part(j)%x
         ryij = ryi - part(j)%y 
         rzij = rzi - part(j)%z

         rxij = rxij - anint(rxij/Lx)*Lx
         ryij = ryij - anint(ryij/Ly)*Ly
         rzij = rzij - anint(rzij/Lz)*Lz 

         rij2 = rxij*rxij + ryij*ryij + rzij*rzij
         if(rij2 .lt. (Rcut*Rcut)) then
           rij = sqrt(rij2)
           Call delta(part(i)%tipo,part(j)%tipo,dd)
           dmin = (part(i)%lambda + part(j)%lambda)*sigma/2.0
           if(rij .lt. sigma) then
             Unb = 100.0
           else
             if(rij .lt. dmin) then
               Unb = dd*epsilon
             else
               Unb = 0.0
             end if
           end if
           E = E + Unb
         end if
         j = list(j)
       End Do
     End Do

     i = list(i)

   End Do

 End Do

Contains

 Subroutine delta(i,j,dd)
  Use MathConstants   
  Implicit None
  
  Integer, Intent(In) :: i
  Integer, Intent(In) :: j 
  Real(8), Intent(Out):: dd
 
  if(i == j) then
    dd = +1.0
  else
    dd = -1.0
  end if

 End Subroutine delta

End Subroutine Hnb


Subroutine MonteCarlo
 Use Parameters
 Use Globals
 Implicit None

 Real(8)  :: ix,iy,iz
 Integer  :: i,flag
 Real(8)  :: d,p,rx,ry,rz,dx,dy,dz,Eo,En,dE
 Real(8)  :: Hf,sum1,sum2,pacc

 Call Random_Number(p)
 i = nint(p*float(Np-1))+1

 ix = part(i)%x
 iy = part(i)%y
 iz = part(i)%z

 Call Energy(i,Eo)

 Call Random_Number(p)
 dx = dmax*(2.0_8*p-1.0_8)
 Call Random_Number(p)
 dy = dmax*(2.0_8*p-1.0_8)
 Call Random_Number(p)
 dz = dmax*(2.0_8*p-1.0_8)
 rx = part(i)%x + dx
 ry = part(i)%y + dy
 rz = part(i)%z + dz
 Call Inbox(rx,ry,rz,flag)
 if(flag == 0) return

 part(i)%x = rx
 part(i)%y = ry
 part(i)%z = rz

 Call Energy(i,En)

 dE  = En - Eo
 pacc= exp(-dE/Temp)
 Call Random_Number(p)
 if(p < pacc) then
   E = E + dE
   cont = cont + 0.01
   Call Linked_Lists
 else
   part(i)%x = ix
   part(i)%y = iy
   part(i)%z = iz
 end if

End Subroutine MonteCarlo


Subroutine Energy(i,U)
 Use Parameters
 Use Globals   
 Implicit None
    
 Integer,  Intent(In) :: i
 Real(8),  Intent(Out):: U
   
 Real(8)  :: Unb,dd,rij,rij2,rxij,ryij,rzij,rxi,ryi,rzi,dmin
 Integer  :: j,k,l,q,ncell,jcell0,icell,jj,jcell,ix,iy,iz,xcell,ycell,zcell,xxcell,yycell,zzcell

 U = 0.0

 k = int(part(i)%x/dLx) + 1
 l = int(part(i)%y/dLy) + 1
 q = int(part(i)%z/dLz) + 1

 rxi = part(i)%x
 ryi = part(i)%y
 rzi = part(i)%z

 Do xcell = k-1,k+1
  Do ycell = l-1,l+1
   Do zcell = q-1,q+1
      Call Ixcell(xcell,ycell,zcell,icell)
      !xxcell = Modulo(xcell,nx)
      !yycell = Modulo(ycell,ny)
      !zzcell = Modulo(zcell,nz)
      !icell  = 1 + xxcell + yycell*nx + zzcell*nx*ny
      j = head(icell)
      Do While(j .gt. 0)
        if(i .ne. j) then
           rxij = rxi - part(j)%x
           ryij = ryi - part(j)%y 
           rzij = rzi - part(j)%z

           rxij = rxij - anint(rxij/Lx)*Lx
           ryij = ryij - anint(ryij/Ly)*Ly
           rzij = rzij - anint(rzij/Lz)*Lz 

           rij2 = rxij*rxij + ryij*ryij + rzij*rzij
           if(rij2 .lt. (Rcut*Rcut)) then
             rij = sqrt(rij2)
             Call delta(part(i)%tipo,part(j)%tipo,dd)
             dmin = (part(i)%lambda + part(j)%lambda)*sigma/2.0
             if(rij .lt. sigma) then
                Unb = 100.0
             else
               if(rij .lt. dmin) then
                 Unb = dd*epsilon
               else
                 Unb = 0.0
               end if
             end if
             U = U + Unb
           end if
        end if
        j = list(j)
      End Do
   End Do
  End Do
 End Do

 
Contains

 Subroutine delta(i,j,dd) 
  Implicit None
  
  Integer, Intent(In) :: i
  Integer, Intent(In) :: j 
  Real(8), Intent(Out):: dd
 
  if(i == j) then
    dd = +1.0
  else
    dd = -1.0
  end if

 End Subroutine delta

 Subroutine Ixcell(ix,iy,iz,icell)
  Use MathConstants   
  Implicit None
  
  Integer, Intent(In)  :: ix
  Integer, Intent(In)  :: iy 
  Integer, Intent(In)  :: iz 
  Integer, Intent(Out) :: icell
 
  icell = 1 + mod(ix-1+nx,nx) + mod(iy-1+ny,ny)*nx + mod(iz-1+nz,nz)*nx*ny 

 End Subroutine Ixcell

End Subroutine Energy
