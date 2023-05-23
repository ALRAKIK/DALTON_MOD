!
      REAL*8     LX, LY, LZ
      REAL*8     ax,ay,az, r 
      REAL*8     parx(2)
      REAL*8     x11, x12
      REAL*8     x21, x22
      LOGICAL    toura , debug ,Omatrix
      CHARACTER*8 ANGO
!
      COMMON /TORUS/  Lx,Ly,Lz,toura,debug,ax,ay,az,r
      COMMON /TORUS/  parx,x11,x12,x21,x22,Omatrix
      COMMON /TORUS/  ANGO 
! 