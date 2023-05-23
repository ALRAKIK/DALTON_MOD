!
!  Dalton, a molecular electronic structure program
!  Copyright (C) The Dalton Authors (see AUTHORS file for details).
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
      module turin 
      
      ! character ! 
      
      CHARACTER*8 ANGO , space 
      
      ! integer ! 
            
      
      ! logical ! 
      
      LOGICAL    toura , debug ,Omatrix
      LOGICAL    First_c , Second_c 
      
      ! real ! 
      
       
      
      REAL*8     LX, LY, LZ
      REAL*8     ax,ay,az
      REAL*8     rx,ry,rz 
      REAL*8     nh1 , nh2 
      
      ! 1 D !
      
      REAL*8     x_par1 , x_par2
      REAL*8     x11, x12
      REAL*8     x21, x22
      INTEGER    ix4c , jx4c 
      REAL*8     ix_center1(5000) , ix_center2(5000)
      REAL*8     jx_center1(5000) , jx_center2(5000) ,x_stupid,y_stupid
      
      
      ! 1Dmm ! 
      
      REAL*8     xx_h , xy_h 
      
      ! 2 D m ! 
      
      REAL*8     y_par1 , y_par2
      REAL*8     y11, y12
      REAL*8     y21, y22    
      INTEGER    iy4c , jy4c 
      REAL*8     iy_center1(5000) , iy_center2(5000)
      REAL*8     jy_center1(5000) , jy_center2(5000) 
       
      
      ! 3 D m ! 
      
      REAL*8     z_par1 , z_par2
      REAL*8     z11, z12
      REAL*8     z21, z22    
      INTEGER    iz4c , jz4c 
      REAL*8     iz_center1(5000) , iz_center2(5000)
      REAL*8     jz_center1(5000) , jz_center2(5000) 
      
      
      
      ! 2 D ! 
      
      REAL*8     R_o, r_i 
      
      REAL*8     Xt1 , Xt2 , Yt1 , Yt2 , Zt1 , Zt2 
      
      REAL*8     x12D , y12D , z12D
      REAL*8     x22D , y22D , z22D
      
      REAL*8, parameter   ::  an = 0.5291772108D0
      
      end module turin 






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! _____                                _                 _                                     _                                _      
!(  _  )                             (_ )               ( )_                         _        ( )_                            (_ )    
!| ( ) |  ___     __             __   | |    __     ___ | ,_) _ __    _     ___     (_)  ___  | ,_)   __     __   _ __    _ _  | |    
!| | | |/' _ `\ /'__`\         /'__`\ | |  /'__`\ /'___)| |  ( '__) /'_`\ /' _ `\   | |/' _ `\| |   /'__`\ /'_ `\( '__) /'_` ) | |    
!| (_) || ( ) |(  ___/        (  ___/ | | (  ___/( (___ | |_ | |   ( (_) )| ( ) |   | || ( ) || |_ (  ___/( (_) || |   ( (_| | | |    
!(_____)(_) (_)`\____)        `\____)(___)`\____)`\____)`\__)(_)   `\___/'(_) (_)   (_)(_) (_)`\__)`\____)`\__  |(_)   `\__,_)(___)   
!                                                                                                         ( )_) |                     
!                                                                                                          \___/'                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DISTANA(x1,x2,X,y1,y2,Y,z1,z2,Z,DIST)
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! it is realted to the overlap matrix, comparsion  between the 
!linear chain with the ring one.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      use turin 
      
      implicit none 

#include "pi.h"
#include "priunit.h"
      
      ! input ! 
      
      REAL*8, intent(in)  :: x1 , x2
      REAL*8, intent(in)  :: y1 , y2
      REAL*8, intent(in)  :: z1 , z2
      
      ! local ! 
      
      REAL*8 :: X
      REAL*8 :: Y
      REAL*8 :: Z
      
      ! output ! 
      
      REAL*8, intent(out) :: DIST
      
      
!      from her1car.F ! 

!     DIFABX = CORAX - CORBX
!     DIFABY = CORAY - CORBY
!     DIFABZ = CORAZ - CORBZ
!     DISTAB = DIFABX*DIFABX + DIFABY*DIFABY + DIFABZ*DIFABZ

      	  
      if (toura) then 
      
        if (space == "1D") then 
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "1Dm") then 
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if
        
         if (space == "1Dmm") then 
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if
        
        if (space == "1Dm2") then 
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if
        
        if (space == "2Dm") then 
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if
        
        if (space == "3Dm") then 
        
          X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
          Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
          Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if
                
        if (space == "2D") then 
        
           Xt1 = (R_o+r_i*cos(ax*x1))*cos(ay*y1)
           Yt1 = (R_o+r_i*cos(ax*x1))*sin(ay*y1)
           Zt1 = r_i*sin(ax*x1)            
                                                   
           Xt2 = (R_o+r_i*cos(ax*x2))*cos(ay*y2)
           Yt2 = (R_o+r_i*cos(ax*x2))*sin(ay*y2)
           Zt2 = r_i*sin(ax*x2)   
                       
           X = Xt1 - Xt2
           Y = Yt1 - Yt2
           Z = Zt1 - Zt2        
           
           if (ANGO == "Angstrom") then 
           
              X =  X/an 
              Y =  Y/an 
              Z =  Z/an
              
           end if 
        
        end if 
        
      else 
            
      X = x1 - x2
      Y = y1 - y2
      Z = z1 - z2 
              
      end if 
      
      DIST = X*X+Y*Y+Z*Z

!      write(lupri,'(a,3f12.6)') 'D  :',DIST
       
      END 
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE BarCENT(e1,e2,
     &                   x1,x2,X,
     &                   y1,y2,Y,
     &                   z1,z2,Z)
      
      use turin 
      
      implicit none 
           
#include "priunit.h"
      
      ! define the variables ! 
      
      real*8,intent(in)       :: e1 , e2 
      real*8,intent(in)       :: x1 , x2
      real*8,intent(in)       :: y1 , y2
      real*8,intent(in)       :: z1 , z2
      
      
      real*8,intent(out)      :: X , Y , Z

      ! local ! 
      
      real*8                  :: theta 
      
      
!       from her1car.F ! 
      
!      CORPX  = EXPAPI*CORAX + EXPBPI*CORBX
!      CORPY  = EXPAPI*CORAY + EXPBPI*CORBY
!      CORPZ  = EXPAPI*CORAZ + EXPBPI*CORBZ
           
      
      if (toura) then 
         
         if (space == "1D") then
         
         theta = ax*x1
         x11 = rx*cos(theta)
         x12 = rx*sin(theta)
         
         theta = ax*x2
         x21 = rx*cos(theta)
         x22 = rx*sin(theta)
         
         x_par1 = e1*x11 + e2*x21
         x_par2 = e1*x12 + e2*x22 
                  
         X  = x_par1
         
         Y  = x_par2 
         
         Z  = e1*z1 + e2*z2
                  
         end if
         
         if (space == "1Dm") then
         
         theta = ax*x1
         x11 = rx*cos(theta)
         x12 = rx*sin(theta)
         
         theta = ax*x2
         x21 = rx*cos(theta)
         x22 = rx*sin(theta)
         
         x_par1 = e1*x11 + e2*x21
         x_par2 = e1*x12 + e2*x22 
                  
         X  = e1*x1 + e2*x2 
         
         Y  = e1*y1 + e2*y2
         
         Z  = e1*z1 + e2*z2
         
               
         end if
         
         if (space == "1Dmm") then
         
         theta = ax*x1
         x11 = rx*cos(theta)
         x12 = rx*sin(theta)
         
         theta = ax*x2
         x21 = rx*cos(theta)
         x22 = rx*sin(theta)
         
         x_par1 = e1*x11 + e2*x21
         x_par2 = e1*x12 + e2*x22 
                  
         X  = e1*x1 + e2*x2 
         
         Y  = e1*y1 + e2*y2
         
         Z  = e1*z1 + e2*z2
               
         end if
         
         if (space == "1Dm2") then
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !           - y -              !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         theta = ay*y1
         y11 = ry*cos(theta)
         y12 = ry*sin(theta)
         
         theta = ay*y2
         y21 = ry*cos(theta)
         y22 = ry*sin(theta)
         
         y_par1 = e1*y11 + e2*y21
         y_par2 = e1*y12 + e2*y22 
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         X  = e1*x1 + e2*x2 
         
         Y  = e1*y1 + e2*y2
         
         Z  = e1*z1 + e2*z2
              
         end if
         
         if (space == "2Dm") then
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !           - x -              !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         theta = ax*x1
         x11 = rx*cos(theta)
         x12 = rx*sin(theta)
         
         theta = ax*x2
         x21 = rx*cos(theta)
         x22 = rx*sin(theta)
         
         x_par1 = e1*x11 + e2*x21
         x_par2 = e1*x12 + e2*x22
         
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !           - y -              !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         theta = ay*y1
         y11 = ry*cos(theta)
         y12 = ry*sin(theta)
         
         theta = ay*y2
         y21 = ry*cos(theta)
         y22 = ry*sin(theta)
         
         y_par1 = e1*y11 + e2*y21
         y_par2 = e1*y12 + e2*y22 
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         X  = e1*x1 + e2*x2 
         
         Y  = e1*y1 + e2*y2
         
         Z  = e1*z1 + e2*z2
   
         end if
         
         if (space == "3Dm") then
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !           - x -              !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         theta = ax*x1
         x11 = rx*cos(theta)
         x12 = rx*sin(theta)
         
         theta = ax*x2
         x21 = rx*cos(theta)
         x22 = rx*sin(theta)
         
         x_par1 = e1*x11 + e2*x21
         x_par2 = e1*x12 + e2*x22
         
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !           - y -              !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         theta = ay*y1
         y11 = ry*cos(theta)
         y12 = ry*sin(theta)
         
         theta = ay*y2
         y21 = ry*cos(theta)
         y22 = ry*sin(theta)
         
         y_par1 = e1*y11 + e2*y21
         y_par2 = e1*y12 + e2*y22 
         
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !           - z -              !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         theta = az*z1
         z11 = rz*cos(theta)
         z12 = rz*sin(theta)
         
         theta = az*z2
         z21 = rz*cos(theta)
         z22 = rz*sin(theta)
         
         z_par1 = e1*z11 + e2*z21
         z_par2 = e1*z12 + e2*z22 
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         X  = e1*x1 + e2*x2 
         
         Y  = e1*y1 + e2*y2
         
         Z  = e1*z1 + e2*z2
   
         end if
         
         if (space == "2D") then 
         
         x12D = (1/an)*(R_o+r_i*cos(ax*x1))*cos(ay*y1)
         y12D = (1/an)*(R_o+r_i*cos(ax*x1))*sin(ay*y1)
         z12D = (1/an)*r_i*sin(ax*x1)            
         
         x22D = (1/an)*(R_o+r_i*cos(ax*x2))*cos(ay*y2)
         y22D = (1/an)*(R_o+r_i*cos(ax*x2))*sin(ay*y2)
         z22D = (1/an)*r_i*sin(ax*x2) 
         
         X  = e1*x12D + e2*x22D
         Y  = e1*y12D + e2*y22D
         Z  = e1*z12D + e2*z22D 
         
         if (ANGO /= "Angstrom") then 
         
           X = X*an 
           Y = Y*an 
           Z = Z*an 
            
         end if 
         
         end if 
                    
      else 
      
      X  = e1*x1 + e2*x2
      Y  = e1*y1 + e2*y2
      Z  = e1*z1 + e2*z2
       
      end if 
            
      END 
 
! done ! 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE DISTARCP(x1,x2,X,y1,y2,Y,z1,z2,Z)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the distance between the nucleus C and the gaussian Px
! Py,Pz as in the book RPC

! this play a role on one_electron hamiltonian integral not the overlap 
! nor the kinetic energy 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use turin 
      
      implicit none 
      
#include "pi.h"
#include "priunit.h"      
      
      REAL*8, intent(in)  :: x1 , x2
      REAL*8, intent(in)  :: y1 , y2
      REAL*8, intent(in)  :: z1 , z2
      
      REAL*8, intent(out) :: X
      REAL*8, intent(out) :: Y
      REAL*8, intent(out) :: Z

!    from her1car.F ! 

!    DIFCPX = COORC(1,IATOMC) - CORPX
!    DIFCPY = COORC(2,IATOMC) - CORPY
!    DIFCPZ = COORC(3,IATOMC) - CORPZ


      ! local ! 
      
      REAL*8              :: theta 

      

      if (toura) then 
        
        if (space == "1D") then 
             
          theta = x1*ax
          
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          X = x11-x_par1 
          Y = x12-x_par2
          
          Z = z1 - z2
                  
        end if
        
        if (space == "1Dm") then 
             
          theta = x1*ax
          
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          X = sqrt((x11-x_par1)**2+(x12-x_par2)**2)
          
          Y = y1 - y2
          Z = z1 - z2
          
                  
        end if
        
          if (space == "1Dmm") then 
             
          theta = x1*ax
          
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
                    
          xx_h = x11-x_par1 
          xy_h = x12-x_par2
          
          X = x1 - x2 
          Y = y1 - y2
          Z = z1 - z2
          
                  
        end if
        
        if (space == "1Dm2") then 
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - y -              !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
          theta = ay*y1
          
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
          
          Y = sqrt((y11-y_par1)**2+(y12-y_par2)**2)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
          
          X = x1 - x2 
          Z = z1 - z2
                  
        end if
        
        if (space == "2Dm") then 
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - x -              !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          theta = ax*x1
          
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          X = sqrt((x11-x_par1)**2+(x12-x_par2)**2)
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - y -              !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
          theta = ay*y1
          
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
          
          Y = sqrt((y11-y_par1)**2+(y12-y_par2)**2)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          Z = z1 - z2
                  
        end if
        
        if (space == "3Dm") then 
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - x -              !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          theta = ax*x1
          
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          X = sqrt((x11-x_par1)**2+(x12-x_par2)**2)
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - y -              !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
          theta = ay*y1
          
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
          
          Y = sqrt((y11-y_par1)**2+(y12-y_par2)**2)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - z -              !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
          theta = az*z1
          
          z11 = rz * cos(theta)
          z12 = rz * sin(theta)
          
          Z = sqrt((z11-z_par1)**2+(z12-z_par2)**2)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
        end if
        
        if (space == "2D") then 
        
          x12D = (1/an)*(R_o+r_i*cos(ax*x1))*cos(ay*y1)
          y12D = (1/an)*(R_o+r_i*cos(ax*x1))*sin(ay*y1)
          z12D = (1/an)*r_i*sin(ax*x1)
        
        if (ANGO /= "Angstrom") then 
         
           x12D = x12D*an 
           y12D = y12D*an 
           z12D = z12D*an 
            
         end if 
        
          X = x12D - x2 
          Y = y12D - y2 
          Z = z12D - z2 
                 
        end if 
                    
!          write(lupri,'(a,3f12.8)') "X distance ", X ,Y, Z     
                    
        else 
        
        X = x1 - x2
        Y = y1 - y2
        Z = z1 - z2
      
!          write(lupri,'(a,3f12.8)') "X distance ", X , Y, Z 
          
      end if 
               
      END
 
! done ! 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NCCOUNT()
      
      use turin 
      
      implicit none 
          
      if (toura) then
            stop 'Non-classical cont. not implemented for TURIN_TORUS'
      end if
      
      END 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DISTANAD(xp,xa,Xpa,
     &                    yp,ya,Ypa,
     &                    zp,za,Zpa,
     &                       xb,Xpb,
     &                       yb,Ypb,
     &                       zb,Zpb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! this subroutine used to calculate the distance between the gaussian 
! and the bary-center which is important into the kin integral 
! and the one - elec hamiltonian integral 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      use turin 
      
      implicit none 
      
#include "pi.h"
#include "priunit.h"      
      
      REAL*8, intent(in)  :: xp , xa ,xb 
      REAL*8, intent(in)  :: yp , ya ,yb
      REAL*8, intent(in)  :: zp , za ,zb
      
      REAL*8, intent(out) :: Xpa , Xpb
      REAL*8, intent(out) :: Ypa , Ypb
      REAL*8, intent(out) :: Zpa , Zpb


      ! local ! 
      
      REAL*8              :: theta 


!    from her1odc.F ! 

!      PAX = CORPX - CORAX
!      PAY = CORPY - CORAY
!      PAZ = CORPZ - CORAZ
!      PBX = CORPX - CORBX
!      PBY = CORPY - CORBY
!      PBZ = CORPZ - CORBZ

!      the coor = center - gaussian 


      if (toura) then 
      
        if (space == "1D") then 
                                
        theta = xa*ax
      
        x11 = rx * cos(theta)
        x12 = rx * sin(theta)
             
        Xpa = x_par1-x11
        Ypa = x_par2-x12
                
        theta = xb*ax
      
        x21 = rx * cos(theta)
        x22 = rx * sin(theta)
              
        Xpb = x_par1-x21
        Ypb = x_par2-x22
                        
        Zpa = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(zp-za))))**(0.5d0)
       
        Zpb = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(zp-zb))))**(0.5d0)
        
        end if

        if (space == "1Dm") then 
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - x -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = ax*xa
      
        x11 = rx * cos(theta)
        x12 = rx * sin(theta)
             
        Xpa = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
                
        theta = ax*xb
      
        x21 = rx * cos(theta)
        x22 = rx * sin(theta)
              
        Xpb = sqrt((x_par1-x21)**2+(x_par2-x22)**2)
        
        Ypa = yp - ya
        Ypb = yp - yb                
        Zpa = zp - za
        Zpb = zp - zb
          
        end if
                
        if (space == "1Dm2") then 
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - y -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = ay*ya
      
        y11 = ry * cos(theta)
        y12 = ry * sin(theta)
             
        Ypa = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
        
        theta = ay*yb
      
        y21 = ry * cos(theta)
        y22 = ry * sin(theta)
              
        Ypb = sqrt((y_par1-y21)**2+(y_par2-y22)**2)
        

        Xpa = xp - xa
        Xpb = xp - xb 
        
        Zpa = zp - za
        Zpb = zp - zb
          
        end if
        
        if (space == "2Dm") then 
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - x -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = ax*xa
      
        x11 = rx * cos(theta)
        x12 = rx * sin(theta)
             
        Xpa = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
        
        theta = ax*xb
      
        x21 = rx * cos(theta)
        x22 = rx * sin(theta)
              
        Xpb = sqrt((x_par1-x21)**2+(x_par2-x22)**2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - y -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = ay*ya
      
        y11 = ry * cos(theta)
        y12 = ry * sin(theta)
             
        Ypa = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
        
        theta = ay*yb
      
        y21 = ry * cos(theta)
        y22 = ry * sin(theta)
              
        Ypb = sqrt((y_par1-y21)**2+(y_par2-y22)**2)
              
        Zpa = zp - za
        Zpb = zp - zb
          
        end if
         
        if (space == "3Dm") then 
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - x -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = ax*xa
      
        x11 = rx * cos(theta)
        x12 = rx * sin(theta)
             
        Xpa = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
        
        theta = ax*xb
      
        x21 = rx * cos(theta)
        x22 = rx * sin(theta)
              
        Xpb = sqrt((x_par1-x21)**2+(x_par2-x22)**2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - y -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = ay*ya
      
        y11 = ry * cos(theta)
        y12 = ry * sin(theta)
             
        Ypa = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
        
        theta = ay*yb
      
        y21 = ry * cos(theta)
        y22 = ry * sin(theta)
              
        Ypb = sqrt((y_par1-y21)**2+(y_par2-y22)**2)
              
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !               - z -               !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        theta = az*za
      
        z11 = rz * cos(theta)
        z12 = rz * sin(theta)
             
        Zpa = sqrt((z_par1-z11)**2+(z_par2-z12)**2)
        
        theta = az*zb
      
        z21 = rz * cos(theta)
        z22 = rz * sin(theta)
              
        Zpb = sqrt((z_par1-z21)**2+(z_par2-z22)**2)
        
        end if 
          
        if (space == "2D") then 
        
          x12D = (1/an)*(R_o+r_i*cos(ax*xa))*cos(ay*ya)
          y12D = (1/an)*(R_o+r_i*cos(ax*xa))*sin(ay*ya)
          z12D = (1/an)*r_i*sin(ax*xa)
        
          if (ANGO /= "Angstrom") then 
         
           x12D = x12D*an 
           y12D = y12D*an 
           z12D = z12D*an 
            
         end if 
          
          Xpa = xp - x12D 
          Ypa = yp - y12D 
          Zpa = zp - z12D 
          
          x22D = (1/an)*(R_o+r_i*cos(ax*xb))*cos(ay*yb)
          y22D = (1/an)*(R_o+r_i*cos(ax*xb))*sin(ay*yb)
          z22D = (1/an)*r_i*sin(ax*xb)
        
          if (ANGO /= "Angstrom") then 
         
           x22D = x22D*an 
           y22D = y22D*an 
           z22D = z22D*an 
            
         end if
                 
          Xpb = xp - x22D 
          Ypb = yp - y22D 
          Zpb = zp - z22D 
        
        
        end if 

      else 
      
        Xpa = xp - xa
        Ypa = yp - ya
        Zpa = zp - za
        
        Xpb = xp - xb
        Ypb = yp - yb
        Zpb = zp - zb
        
      end if 
                        
      END

! done !       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! _____                               _                 _                                     _                                _      
!(_   _)                            (_ )               ( )_                         _        ( )_                            (_ )    
!  | |   _   _   _    _         __   | |    __     ___ | ,_) _ __    _     ___     (_)  ___  | ,_)   __     __   _ __    _ _  | |    
!  | |  ( ) ( ) ( ) /'_`\     /'__`\ | |  /'__`\ /'___)| |  ( '__) /'_`\ /' _ `\   | |/' _ `\| |   /'__`\ /'_ `\( '__) /'_` ) | |    
!  | |  | \_/ \_/ |( (_) )   (  ___/ | | (  ___/( (___ | |_ | |   ( (_) )| ( ) |   | || ( ) || |_ (  ___/( (_) || |   ( (_| | | |    
!  (_)  `\___x___/'`\___/'   `\____)(___)`\____)`\____)`\__)(_)   `\___/'(_) (_)   (_)(_) (_)`\__)`\____)`\__  |(_)   `\__,_)(___)   
!                                                                                                        ( )_) |                     
!                                                                                                         \___/'                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The Overlap distributions for both electron, check ODCVEC in her2drv.F
! which is in the her2odc.F  and her2her.F 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE DISTANAX1X2(x1,x2,X,
     &                       y1,y2,Y,
     &                       z1,z2,Z,
     &                       D)
      
      use turin 
      
      implicit none 
      
#include "pi.h"
#include "priunit.h"
     
      REAL*8, intent(in)  :: x1 , x2
      REAL*8, intent(in)  :: y1 , y2
      REAL*8, intent(in)  :: z1 , z2
      
      ! local ! 
      
      REAL*8, intent(out) :: X
      REAL*8, intent(out) :: Y
      REAL*8, intent(out) :: Z
      
      
      ! output ! 
      
      REAL*8, intent(out) :: D
      
!      from her2odc.F 

!      DIFX = CRX1 - CRX2
!      DIFY = CRY1 - CRY2
!      DIFZ = CRZ1 - CRZ2

               
      if (toura) then 
        
        if (space == "1D") then 
      
        X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
        Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
        Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "1Dm") then 
      
        X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
        Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
        Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "1Dmm") then 
      
        X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
        Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
        Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "1Dm2") then 
      
        X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
        Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
        Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "2Dm") then 
      
        X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
        Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
        Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "3Dm") then 
      
        X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
        Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
        Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
        
        end if 
        
        if (space == "2D") then 
        
           Xt1 = (R_o+r_i*cos(ax*x1))*cos(ay*y1)
           Yt1 = (R_o+r_i*cos(ax*x1))*sin(ay*y1)
           Zt1 = r_i*sin(ax*x1)            
                                                   
           Xt2 = (R_o+r_i*cos(ax*x2))*cos(ay*y2)
           Yt2 = (R_o+r_i*cos(ax*x2))*sin(ay*y2)
           Zt2 = r_i*sin(ax*x2)   
                       
           X = Xt1 - Xt2
           Y = Yt1 - Yt2
           Z = Zt1 - Zt2        
           
           if (ANGO == "Angstrom") then 
           
              X =  X/an 
              Y =  Y/an 
              Z =  Z/an
              
           end if 
        
        end if 
        
      else 
      
        X =  x1 - x2
        Y =  y1 - y2
        Z =  z1 - z2
        
      end if 
          

        D = X*X + Y*Y + Z*Z
          
!      write(lupri,'(a,3f12.6))') 'D :',D
!      write(lupri,'()')
       
      END 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE COORP12(e1,e2,
     &                   x1,x2,X,
     &                   y1,y2,Y,
     &                   z1,z2,Z)
      
      
      use turin 
      
      implicit none
      
#include "pi.h"
#include "priunit.h" 


      REAL*8, intent(in)  :: e1 , e2
      REAL*8, intent(in)  :: x1 , x2
      REAL*8, intent(in)  :: y1 , y2
      REAL*8, intent(in)  :: z1 , z2
      
      REAL*8, intent(out) :: X
      REAL*8, intent(out) :: Y
      REAL*8, intent(out) :: Z
       
      ! local !  
      
      REAL*8              :: theta 
       
       
!       from her2odc.F
      
!       CORPX  = EXP1PI*CRX1 + EXP2PI*CRX2
!       CORPY  = EXP1PI*CRY1 + EXP2PI*CRY2
!       CORPZ  = EXP1PI*CRZ1 + EXP2PI*CRZ2
      
      
      if (toura) then 
        
        if (space == "1D") then 
        
          theta = ax*x1
          x11 = rx*cos(theta)
          x12 = rx*sin(theta)
        
          theta = ax*x2
          x21 = rx*cos(theta)
          x22 = rx*sin(theta)
        
          x_par1 = e1*x11 + e2*x21
          x_par2 = e1*x12 + e2*x22
        
          X  = x_par1
          Y  = x_par2        
          Z  = e1*z1 + e2*z2

        end if 
        
        if (space == "1Dm") then 
        
          theta = ax*x1
          x11 = rx*cos(theta)
          x12 = rx*sin(theta)
        
          theta = ax*x2
          x21 = rx*cos(theta)
          x22 = rx*sin(theta)
          
          if (Second_c) then 
            jx4c = jx4c + 1 
            jx_center1(jx4c) = e1*x11 + e2*x21
            jx_center2(jx4c) = e1*x12 + e2*x22           
          end if 
          
          if (First_c) then 
            ix4c = ix4c + 1 
            ix_center1(ix4c) = e1*x11 + e2*x21
            ix_center2(ix4c) = e1*x12 + e2*x22
          end if
          
          x_stupid = e1*x11 + e2*x21
          y_stupid = e1*x12 + e2*x22
          
          X  = e1*x1 + e2*x2
          Y  = e1*y1 + e2*y2        
          Z  = e1*z1 + e2*z2
          
        end if 
        
        if (space == "1Dmm") then 
        
          theta = ax*x1
          x11 = rx*cos(theta)
          x12 = rx*sin(theta)
        
          theta = ax*x2
          x21 = rx*cos(theta)
          x22 = rx*sin(theta)
          
          if (Second_c) then 
            jx4c = jx4c + 1 
            jx_center1(jx4c) = e1*x11 + e2*x21
            jx_center2(jx4c) = e1*x12 + e2*x22           
          end if 
          
          if (First_c) then 
            ix4c = ix4c + 1 
            ix_center1(ix4c) = e1*x11 + e2*x21
            ix_center2(ix4c) = e1*x12 + e2*x22
          end if
          
          x_stupid = e1*x11 + e2*x21
          y_stupid = e1*x12 + e2*x22
          
          X  = e1*x1 + e2*x2
          Y  = e1*y1 + e2*y2        
          Z  = e1*z1 + e2*z2
          
        end if 
        
        if (space == "1Dm2") then 
        
          theta = ay*y1
          y11 = ry*cos(theta)
          y12 = ry*sin(theta)
        
          theta = ay*y2
          y21 = ry*cos(theta)
          y22 = ry*sin(theta)
          
          if (Second_c) then 
            jy4c = jy4c + 1 
            jy_center1(jy4c) = e1*y11 + e2*y21
            jy_center2(jy4c) = e1*y12 + e2*y22
          end if 
          
          if (First_c) then 
            iy4c = iy4c + 1 
            iy_center1(iy4c) = e1*y11 + e2*y21
            iy_center2(iy4c) = e1*y12 + e2*y22            
          end if
          
          X  = e1*x1 + e2*x2
          Y  = e1*y1 + e2*y2        
          Z  = e1*z1 + e2*z2
          
        end if 
        
        
         if (space == "2Dm") then 

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - x -             !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          theta = ax*x1
          x11 = rx*cos(theta)
          x12 = rx*sin(theta)
        
          theta = ax*x2
          x21 = rx*cos(theta)
          x22 = rx*sin(theta)
          
          if (Second_c) then 
            jx4c = jx4c + 1 
            jx_center1(jx4c) = e1*x11 + e2*x21
            jx_center2(jx4c) = e1*x12 + e2*x22           
          end if 
          
          if (First_c) then 
            ix4c = ix4c + 1 
            ix_center1(ix4c) = e1*x11 + e2*x21
            ix_center2(ix4c) = e1*x12 + e2*x22
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - y -             !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
          theta = ay*y1
          y11 = ry*cos(theta)
          y12 = ry*sin(theta)
        
          theta = ay*y2
          y21 = ry*cos(theta)
          y22 = ry*sin(theta)
          
          if (Second_c) then 
            jy4c = jy4c + 1 
            jy_center1(jy4c) = e1*y11 + e2*y21
            jy_center2(jy4c) = e1*y12 + e2*y22
          end if 
          
          if (First_c) then 
            iy4c = iy4c + 1 
            iy_center1(iy4c) = e1*y11 + e2*y21
            iy_center2(iy4c) = e1*y12 + e2*y22            
          end if
          
          X  = e1*x1 + e2*x2
          Y  = e1*y1 + e2*y2        
          Z  = e1*z1 + e2*z2
          
        end if 
       
        if (space == "3Dm") then 

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - x -             !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          theta = ax*x1
          x11 = rx*cos(theta)
          x12 = rx*sin(theta)
        
          theta = ax*x2
          x21 = rx*cos(theta)
          x22 = rx*sin(theta)
          
          if (Second_c) then 
            jx4c = jx4c + 1 
            jx_center1(jx4c) = e1*x11 + e2*x21
            jx_center2(jx4c) = e1*x12 + e2*x22           
          end if 
          
          if (First_c) then 
            ix4c = ix4c + 1 
            ix_center1(ix4c) = e1*x11 + e2*x21
            ix_center2(ix4c) = e1*x12 + e2*x22
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - y -             !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
          theta = ay*y1
          y11 = ry*cos(theta)
          y12 = ry*sin(theta)
        
          theta = ay*y2
          y21 = ry*cos(theta)
          y22 = ry*sin(theta)
          
          if (Second_c) then 
            jy4c = jy4c + 1 
            jy_center1(jy4c) = e1*y11 + e2*y21
            jy_center2(jy4c) = e1*y12 + e2*y22
          end if 
          
          if (First_c) then 
            iy4c = iy4c + 1 
            iy_center1(iy4c) = e1*y11 + e2*y21
            iy_center2(iy4c) = e1*y12 + e2*y22            
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           - z -             !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    
          theta = az*z1
          z11 = rz*cos(theta)
          z12 = rz*sin(theta)
        
          theta = az*z2
          z21 = rz*cos(theta)
          z22 = rz*sin(theta)
          
          if (Second_c) then 
            jz4c = jz4c + 1 
            jz_center1(jz4c) = e1*z11 + e2*z21
            jz_center2(jz4c) = e1*z12 + e2*z22
          end if 
          
          if (First_c) then 
            iz4c = iz4c + 1 
            iz_center1(iz4c) = e1*z11 + e2*z21
            iz_center2(iz4c) = e1*z12 + e2*z22            
          end if
          
          X  = e1*x1 + e2*x2
          Y  = e1*y1 + e2*y2        
          Z  = e1*z1 + e2*z2
          
        end if
        
        
        if (space == "2D") then 
         
         x12D = (1/an)*(R_o+r_i*cos(ax*x1))*cos(ay*y1)
         y12D = (1/an)*(R_o+r_i*cos(ax*x1))*sin(ay*y1)
         z12D = (1/an)*r_i*sin(ax*x1)            
         
         x22D = (1/an)*(R_o+r_i*cos(ax*x2))*cos(ay*y2)
         y22D = (1/an)*(R_o+r_i*cos(ax*x2))*sin(ay*y2)
         z22D = (1/an)*r_i*sin(ax*x2) 
         
         X  = e1*x12D + e2*x22D
         Y  = e1*y12D + e2*y22D
         Z  = e1*z12D + e2*z22D 
         
          if (ANGO /= "Angstrom") then 
         
           X = X*an 
           Y = Y*an 
           Z = Z*an 
            
          end if 
         
         end if
            
      else 
      
        X  = e1*x1 + e2*x2
        Y  = e1*y1 + e2*y2
        Z  = e1*z1 + e2*z2
            
      end if 
            
      END 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE COORD12(xp,xa,Xpa,
     &                   yp,ya,Ypa,
     &                   zp,za,Zpa,
     &                      xb,Xpb,
     &                      yb,Ypb,
     &                      zb,Zpb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
! this subroutine used to calculate the distance between the gaussian 
! and the bery-center for the screening effect 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      use turin 
      
      implicit none 
      

#include "pi.h"
#include "priunit.h"      
      
      REAL*8, intent(in)  :: xp , xa ,xb 
      REAL*8, intent(in)  :: yp , ya ,yb
      REAL*8, intent(in)  :: zp , za ,zb
      
      REAL*8, intent(out) :: Xpa , Xpb
      REAL*8, intent(out) :: Ypa , Ypb
      REAL*8, intent(out) :: Zpa , Zpb


      ! local ! 
      
      REAL*8              :: theta 


!    from her2odc.F ! 

!      COOR12(IPRIM,1,2,NSYMOP) = CORPX - CRX1
!      COOR12(IPRIM,2,2,NSYMOP) = CORPY - CRY1
!      COOR12(IPRIM,3,2,NSYMOP) = CORPZ - CRZ1
!      COOR12(IPRIM,1,3,NSYMOP) = CORPX - CRX2
!      COOR12(IPRIM,2,3,NSYMOP) = CORPY - CRY2
!      COOR12(IPRIM,3,3,NSYMOP) = CORPZ - CRZ2

!      the coor = center - gaussian 


      if (toura) then 
      
        if (space == "1D") then 
        
          theta = xa*ax
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
        
          Xpa = x_par1 - x11
          Ypa = x_par2 - x12 
                    
          theta = xb*ax
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
        
          Xpb = x_par1 - x11
          Ypb = x_par2 - x12 
                  
          Zpa = zp - za
          Zpb = zp - zb 
          
        end if
        
        if (space == "1Dm") then 
        
          theta = ax*xa
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          Xpa = sqrt((x_stupid-x11)**2+(y_stupid-x12)**2)
          
          theta = ax*xb
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
         
          Xpb = sqrt((x_stupid-x11)**2+(y_stupid-x12)**2)
          
          Ypa = yp - ya
          Ypb = yp - yb         
          Zpa = zp - za
          Zpb = zp - zb 
          
          
        end if
        
        if (space == "1Dmm") then 
        
          theta = ax*xa
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          Xpa = sqrt((x_stupid-x11)**2+(y_stupid-x12)**2)
          
          theta = ax*xb
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
         
          Xpb = sqrt((x_stupid-x11)**2+(y_stupid-x12)**2)
          
          Ypa = yp - ya
          Ypb = yp - yb         
          Zpa = zp - za
          Zpb = zp - zb 
          
          
        end if
        
        if (space == "1Dm2") then 
        
          theta = ay*ya
        
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
          
          Ypa = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
          
          theta = ay*yb
        
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
         
          Ypb = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
                    
          Xpa = xp - xa
          Xpb = xp - xb         
          Zpa = zp - za
          Zpb = zp - zb 
          
        end if
        
        if (space == "2Dm") then 
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          !          - x -         !
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          
          theta = ax*xa
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          Xpa = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
          
          theta = ax*xb
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
         
          Xpb = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
            
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          !          - y -         !
          !!!!!!!!!!!!!!!!!!!!!!!!!!
        
          theta = ay*ya
        
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
          
          Ypa = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
          
          theta = ay*yb
        
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
         
          Ypb = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
                           
          Zpa = zp - za
          Zpb = zp - zb 
          
        end if
        
        if (space == "3Dm") then 
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          !          - x -         !
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          
          theta = ax*xa
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
          
          Xpa = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
          
          theta = ax*xb
        
          x11 = rx * cos(theta)
          x12 = rx * sin(theta)
         
          Xpb = sqrt((x_par1-x11)**2+(x_par2-x12)**2)
            
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          !          - y -         !
          !!!!!!!!!!!!!!!!!!!!!!!!!!
        
          theta = ay*ya
        
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
          
          Ypa = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
          
          theta = ay*yb
        
          y11 = ry * cos(theta)
          y12 = ry * sin(theta)
         
          Ypb = sqrt((y_par1-y11)**2+(y_par2-y12)**2)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          !          - z -         !
          !!!!!!!!!!!!!!!!!!!!!!!!!!
        
          theta = az*za
        
          z11 = rz * cos(theta)
          z12 = rz * sin(theta)
          
          Zpa = sqrt((z_par1-z11)**2+(z_par2-z12)**2)
          
          theta = az*zb
        
          z11 = rz * cos(theta)
          z12 = rz * sin(theta)
         
          Zpb = sqrt((z_par1-z11)**2+(z_par2-z12)**2)
          
        end if
        
        if (space == "2D") then 
        
          x12D = (1/an)*(R_o+r_i*cos(ax*xa))*cos(ay*ya)
          y12D = (1/an)*(R_o+r_i*cos(ax*xa))*sin(ay*ya)
          z12D = (1/an)*r_i*sin(ax*xa)
        
           if (ANGO /= "Angstrom") then 
         
            x12D = x12D*an 
            y12D = y12D*an 
            z12D = z12D*an 
            
           end if 
        
          Xpa = xp - x12D 
          Ypa = yp - y12D 
          Zpa = zp - z12D 
        
          x22D = (1/an)*(R_o+r_i*cos(ax*xb))*cos(ay*yb)
          y22D = (1/an)*(R_o+r_i*cos(ax*xb))*sin(ay*yb)
          z22D = (1/an)*r_i*sin(ax*xb)
        
           if (ANGO /= "Angstrom") then 
         
            x22D = x22D*an 
            y22D = y22D*an 
            z22D = z22D*an 
            
           end if         
        
            Xpb = xp - x22D 
            Ypb = yp - y22D 
            Zpb = zp - z22D 
        
        
        end if 
        
!        write(lupri,'(a,3f12.8)') 'Distatnce  ' , Xpa , Ypa , Zpa  
         
      else 
      
        Xpa = xp - xa
        Ypa = yp - ya
        Zpa = zp - za
        
        Xpb = xp - xb
        Ypb = yp - yb
        Zpb = zp - zb
        
!        write(lupri,'(a,3f12.8)') 'Distatnce  ' , Xpa , Ypa , Zpa  
        
      end if 
                        
      END
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      _                                                      _                          
!                     (_ )                                                   (_ )        _               
!  ___   _   _    ___  | |    __     _ _  _ __     _ __    __   _ _    _   _  | |   ___ (_)   _     ___  
!/' _ `\( ) ( ) /'___) | |  /'__`\ /'_` )( '__)   ( '__) /'__`\( '_`\ ( ) ( ) | | /',__)| | /'_`\ /' _ `\
!| ( ) || (_) |( (___  | | (  ___/( (_| || |      | |   (  ___/| (_) )| (_) | | | \__, \| |( (_) )| ( ) |
!(_) (_)`\___/'`\____)(___)`\____)`\__,_)(_)      (_)   `\____)| ,__/'`\___/'(___)(____/(_)`\___/'(_) (_)
!                                                              | |                                       
!                                                              (_)                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          SUBROUTINE DISTNUC(x1,x2,y1,y2,z1,z2,D)
          
          use turin 
          
          implicit none 

#include "pi.h"
#include "priunit.h"
      
          REAL*8, intent(in)  :: x1 , x2
          REAL*8, intent(in)  :: y1 , y2
          REAL*8, intent(in)  :: z1 , z2
          
            
          REAL*8, intent(out) :: D
          
          ! local ! 
          
          REAL*8              :: X   , Y   , Z 
          
              
          if (toura) then 
            
             if (space == "1D") then 
                       
              X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
              Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
              Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
            
             end if 
             
              if (space == "1Dm") then 
                       
              X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
              Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
              Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
                          
              end if 
              
              if (space == "1Dmm") then 
                       
              X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
              Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
              Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
                          
              end if 
             
              if (space == "1Dm2") then 
                       
              X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
              Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
              Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
                          
              end if
              
              if (space == "2Dm") then 
                       
              X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
              Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
              Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
                          
              end if
              
              if (space == "3Dm") then 
                       
              X = (ax**(-2.d0)*(2.d0 - 2.d0*cos(ax*(x1-x2))))**(0.5d0)
              Y = (ay**(-2.d0)*(2.d0 - 2.d0*cos(ay*(y1-y2))))**(0.5d0)
              Z = (az**(-2.d0)*(2.d0 - 2.d0*cos(az*(z1-z2))))**(0.5d0)
                          
              end if
             
             
             if (space == "2D") then 
                                
                  Xt1 = (R_o+r_i*cos(ax*x1))*cos(ay*y1)
                  Yt1 = (R_o+r_i*cos(ax*x1))*sin(ay*y1)
                  Zt1 = r_i*sin(ax*x1)            
                                                          
                  Xt2 = (R_o+r_i*cos(ax*x2))*cos(ay*y2)
                  Yt2 = (R_o+r_i*cos(ax*x2))*sin(ay*y2)
                  Zt2 = r_i*sin(ax*x2)   
                              
                  X = Xt1 - Xt2
                  Y = Yt1 - Yt2
                  Z = Zt1 - Zt2
                         
                if (ANGO == "Angstrom") then
                              
                  X =  X/an 
                  Y =  Y/an 
                  Z =  Z/an
                  
                end if 
                
             end if 
            
          else 
                  
          X = x1 - x2
          Y = y1 - y2
          Z = z1 - z2 
                                                      
          end if 
          
          D = X*X+Y*Y+Z*Z
          
          !write(lupri,'(a,3f12.8)') "Distance  :", X/10 , Y/10 , Z/10 
          !write(lupri,'(a,f12.8)') "Distance  :", D/10
          
          END 
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! change the coordinates to circle one 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          SUBROUTINE CORDCHN(X,Y,Z)
          

          use turin 
          
          implicit none 

#include "pi.h"
#include "priunit.h"
      
          REAL*8, intent(inout)  :: X, Y , Z  
                    
          if (toura) then

            if (space == "1D") then 
              Y = rx*sin(ax*X)
              X = rx*cos(ax*X)
            end if 
            
            if (space == "1Dm") then 
              Y = rx*sin(ax*X)
              X = rx*cos(ax*X)
            end if 
            
            if (space == "1Dmm") then 
              Y = rx*sin(ax*X)
              X = rx*cos(ax*X)
            end if 
            
            if (space == "1Dm2") then 
              Y = ry*sin(ay*X)
              X = ry*cos(ay*X)
              X = sqrt(X**2+Y**2)
              Y = 0.d0
            end if 
            
            if (space == "2Dm") then 
              Y = ry*sin(ax*X)
              X = ry*cos(ax*X)
              Y = sqrt(X**2+Y**2)
              Y = 0.d0
            end if
            
            if (space == "3Dm") then 
              Y = ry*sin(ax*X)
              X = ry*cos(ax*X)
              Y = sqrt(X**2+Y**2)
              Y = 0.d0
            end if
                        
            if (space == "2D") then 
            
              Xt1 = (R_o+r_i*cos(ax*X))*cos(ay*Y)
              Yt1 = (R_o+r_i*cos(ax*X))*sin(ay*Y)
              Zt1 = r_i*sin(ax*X)
              
              X = Xt1 
              Y = Yt1 
              Z = Zt1 
              
            end if 
                      
          else 
                  
          
          end if 
          
!          write(lupri,'(a,3f12.8)') "X,Y,Z  ", X , Y , Z
          
          END 
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




