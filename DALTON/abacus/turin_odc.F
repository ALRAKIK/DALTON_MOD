
C  /* Deck getodc */
      SUBROUTINE GETODC_Tx(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,DIFODC
     &               ,KINODC,ONECEN,EXPA,EXPB,IPRINT,SAAB13,EXPPI,WORK
     &               ,LWORK,CORPX,CORPY,CORPZ,DONUC1
     &               ,DOMOM1,ORIGIN,INTTYP)
     
#include "implicit.h"
#include "priunit.h"
      LOGICAL DIFODC, KINODC, ONECEN, DONUC1, DOMOM1
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &          WORK(LWORK), ORIGIN(3)
C
      IF (DONUC1) THEN
         JMAXAP = JMAXA + JMAXD
         JMAXBP = JMAXB
      ELSE
         JMAXAP = JMAXA
         JMAXBP = JMAXB + JMAXD
      END IF
      IF (DOMOM1) THEN
         JMAXAP = JMAXAP + JMAXM
      ELSE
         JMAXBP = JMAXBP + JMAXM
      END IF
      ITEX = MAX(1,JMAXD + JMAXM)
      JMAXTP = JMAXT + ITEX
C
      IF ((INTTYP .EQ. 49) .OR. (INTTYP .EQ. 51)) THEN
         IF (DOMOM1) THEN
            JMAXBP = JMAXBP + 1
         ELSE
            JMAXAP = JMAXAP + 1
         END IF
         JMAXTP = JMAXTP +1
      ENDIF
C
      KODCUN = 1
      KLAST  = 1 + 3*(JMAXAP + 4)*(JMAXBP + 3)*(JMAXTP + 3)*(JMAXD + 1)
     &              *(JMAXM + 1)
      IF (KLAST .GT. LWORK) CALL STOPIT('GETODC',' ',KLAST,LWORK)
      CALL GETOD1_Tx(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAXAP,JMAXBP,
     &            JMAXTP,WORK(KODCUN),DIFODC,KINODC,ONECEN,EXPA,EXPB,
     &            IPRINT,SAAB13,EXPPI,CORPX,CORPY,CORPZ,DONUC1,DOMOM1,
     &            ORIGIN,INTTYP)
      RETURN
      END
C  /* Deck getod1 */
      SUBROUTINE GETOD1_Tx(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,JMAXAP
     &            ,JMAXBP,JMAXTP,ODCUND,DIFODC,KINODC,ONECEN,EXPA,EXPB,
     &                  IPRINT,SAAB13,EXPPI,CORPX,CORPY,CORPZ,DONUC1,
     &                  DOMOM1,ORIGIN,INTTYP)
#include "implicit.h"
#include "priunit.h"
      LOGICAL DIFODC, KINODC, ONECEN, DONUC1, DOMOM1
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3),
     &      ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,0:JMAXD,0:JMAXM,3)
C
C     Clear arrays
C
      CALL DZERO(ODC,3*(JMAXA  + 1)*(JMAXB + 1)*(JMAXT + 1)*(JMAXD + 1)
     &                *(JMAXM + 1))
      CALL DZERO(ODCUND,3*(JMAXAP + 4)*(JMAXBP + 3)*(JMAXTP+3)
     &                   *(JMAXD+1)*(JMAXM + 1))
C
C     Undifferentiated expansion coefficients
C
      CALL ONEODC_Tx(ODCUND,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,SAAB13
     &       ,EXPPI,CORPX,CORPY,CORPZ,IPRINT)
      CALL COPODC_Tx(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT
     &            ,JMAXTP,JMAXD,JMAXM,0,0)
C
C     Expansion coefficients for derivatives
C
      IF (DIFODC) THEN
         IF (.NOT.ONECEN) THEN
            IF (DONUC1) THEN
               CALL DODCA(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                    JMAXTP,JMAXD,JMAXM,EXPA,IPRINT)
            ELSE
               CALL DODCB(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                    JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
            END IF
            IF (INTTYP .EQ. 51) CALL DEFODC(ODC,ODCUND,JMAXA,JMAXAP,
     $           JMAXB,JMAXBP,JMAXT,JMAXTP,JMAXD,JMAXM,EXPA,EXPB,
     $           IPRINT)
         END IF
      ELSE IF (KINODC) THEN
            CALL TODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,EXPA)
      END IF
C
C     Expansion coeffiecients for moments or electric derivatives
C
      IF (JMAXM .GT. 0) THEN
         IF ((INTTYP .EQ. 49) .OR. (INTTYP .EQ. 51)) THEN
            CALL EFODC(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,EXPA,EXPB,IPRINT)
         ELSE IF (INTTYP .EQ. 63 .OR. INTTYP .EQ. 83) THEN
            CALL DODCB2(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                  JMAXTP,JMAXD,JMAXM,EXPB,IPRINT)
         ELSE IF (DOMOM1) THEN
            CALL MODCA(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,ORIGIN,IPRINT,INTTYP)
         ELSE
            CALL MODCB(ODC,ODCUND,JMAXA,JMAXAP,JMAXBP,JMAXB,JMAXT,
     &                 JMAXTP,JMAXD,JMAXM,ORIGIN,IPRINT,INTTYP)
         END IF
      END IF
      RETURN
      END
C  /* Deck copodc */
      SUBROUTINE COPODC_TX(ODCUND,ODC,JMAXAP,JMAXBP,JMAXA,JMAXB,JMAXT,
     &                JMAXTP,JMAXD,JMAXM,ICOPYD,ICOPYM)
C
C     Copy expansion coefficients from ODCUND to ODC
C
#include "implicit.h"
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3),
     &          ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      DO 100 IC = 1, 3
         DO 100 IM = 0, ICOPYM
         DO 100 ID = 0, ICOPYD
            DO 100 IT = 0, JMAXT
               DO 100 IB = 0, JMAXB
               DO 100 IA = 0, JMAXA
                  ODC(IA,IB,IT,ID,IM,IC) = ODCUND(IA,IB,IT,ID,IM,IC)
  100 CONTINUE
      RETURN
      END
C  /* Deck oneodc */
      SUBROUTINE ONEODC_Tx(ODCUND,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,FAC,
     &                  EXPPI,CORPX,CORPY,CORPZ,IPRINT)
     
      use turin 
C
C     TUH 91
C      
#include "implicit.h"
#include "priunit.h"
      INTEGER T, A, B, AB
      DIMENSION ODCUND(-2:JMAXAP+1,-1:JMAXBP+1,-2:JMAXTP,
     &                 0:JMAXD,0:JMAXM,3)
#include "onecom.h"
C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! turin !  overlap matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! original !

!      PAX = CORPX - CORAX
!      PAY = CORPY - CORAY
!      PAZ = CORPZ - CORAZ
!      PBX = CORPX - CORBX
!      PBY = CORPY - CORBY
!      PBZ = CORPZ - CORBZ
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
       theta = ax*CORAX
      
       x11 = rx * cos(theta)
       x12 = rx * sin(theta)
            
       PAX = CORPX-x11
       PAY = CORPY-x12
               
       theta = ax*CORPY
      
       x21 = rx * cos(theta)
       x22 = rx * sin(theta)
             
       PBX = CORPX-x21
       PBY = CORPY-x22

       PAZ = CORPZ - CORAZ
       PBZ = CORPZ - CORBZ
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
      EXPPIH = EXPPI/2
      DO 100 A = 0, JMAXAP
         IF (A .EQ. 0) THEN
            ODCUND(0,0,0,0,0,1) = FAC
            ODCUND(0,0,0,0,0,2) = FAC
            ODCUND(0,0,0,0,0,3) = FAC                    
         ELSE
            DO 200 T = 0, A
               X_TP1 = T + 1
               
               ODCUND(A,0,T,0,0,1) = EXPPIH*ODCUND(A-1,0,T-1,0,0,1)
     *                                + PAX*ODCUND(A-1,0,T  ,0,0,1)
     *                              + X_TP1*ODCUND(A-1,0,T+1,0,0,1)
               ODCUND(A,0,T,0,0,2) = EXPPIH*ODCUND(A-1,0,T-1,0,0,2)
     *                                + PAY*ODCUND(A-1,0,T  ,0,0,2)
     *                              + X_TP1*ODCUND(A-1,0,T+1,0,0,2)
               ODCUND(A,0,T,0,0,3) = EXPPIH*ODCUND(A-1,0,T-1,0,0,3)
     *                                + PAZ*ODCUND(A-1,0,T  ,0,0,3)
     *                              + X_TP1*ODCUND(A-1,0,T+1,0,0,3)
                    
  200       CONTINUE
         END IF
         DO 300 B = 1, JMAXBP
            AB = A + B
C
            DO 400 T = 0, AB
               X_TP1 = T + 1
               
               ODCUND(A,B,T,0,0,1) = EXPPIH*ODCUND(A,B-1,T-1,0,0,1)
     &                                + PBX*ODCUND(A,B-1,T  ,0,0,1)
     &                              + X_TP1*ODCUND(A,B-1,T+1,0,0,1)
               ODCUND(A,B,T,0,0,2) = EXPPIH*ODCUND(A,B-1,T-1,0,0,2)
     &                                + PBY*ODCUND(A,B-1,T  ,0,0,2)
     &                              + X_TP1*ODCUND(A,B-1,T+1,0,0,2)
               ODCUND(A,B,T,0,0,3) = EXPPIH*ODCUND(A,B-1,T-1,0,0,3)
     &                                + PBZ*ODCUND(A,B-1,T  ,0,0,3)
     &                              + X_TP1*ODCUND(A,B-1,T+1,0,0,3)     
 
     
  400          CONTINUE
  300    CONTINUE
  100 CONTINUE
      RETURN
      END
      
      
      
      
      
      
C  /* Deck hernai */
      SUBROUTINE HERNAI_Tx(AHGTF,JMAX,PVAL,CPX,CPY,CPZ,FACINT,IADR,
     &                  ISTEPU,ISTEPV,NRTUV,IPRINT)
#include "implicit.h"
C
#include "priunit.h"
#include "pi.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER (D1 = 1.D0, D2 = 2.D0, D2PI = D2*PI)
      DIMENSION AHGTF(*)
#include "gamcom.h"

C
C     *************************************
C     ***** Incomplete Gamma Function *****
C     *************************************
C
      WVAL = PVAL*(CPX*CPX + CPY*CPY + CPZ*CPZ)
      
      JMAX0 = JMAX
      CALL GAMFUN
C
C     **********************************
C     ***** Special Case: JMAX = 0 *****
C     **********************************
C
      IF (JMAX .EQ. 0) THEN
         AHGTF(IADR) = FACINT*D2PI*FJW(0)/PVAL
         RETURN
      END IF
C
C     **********************************
C     ***** General Case: JMAX > 0 *****
C     **********************************
C
      IF (IAND(JMAX,1) .EQ. 0) THEN
         ISTRTJ = IADR
         ISTEPJ = NRTUV
      ELSE
         ISTRTJ = IADR + NRTUV
         ISTEPJ = - NRTUV
      END IF
      D2PVAL = PVAL + PVAL
      FACTOR = FACINT*D2PI/PVAL
      DO 100 JVAL = 0, JMAX
         FJW(JVAL) = FACTOR*FJW(JVAL)
         FACTOR    = - D2PVAL*FACTOR
  100 CONTINUE
C
C     ***** JVAL = 0 *****
C
      AHGTF(ISTRTJ) = FJW(JMAX)
C
C     ***** JVAL = 1 *****
C
      ISTRTJ                 =   ISTRTJ + ISTEPJ
      FJWMAX                 =   FJW(JMAX)
      AHGTF(ISTRTJ)          =   FJW(JMAX - 1)
      AHGTF(ISTRTJ +      1) = - CPX*FJWMAX
      AHGTF(ISTRTJ + ISTEPU) = - CPY*FJWMAX
      AHGTF(ISTRTJ + ISTEPV) = - CPZ*FJWMAX
C
      IF (JMAX .GT. 1) THEN
C
C        ***** JVAL > 1 *****
C
         ISTPTU = 1 - ISTEPU
         ISTEPJ =   - ISTEPJ
         DO 200 JVAL = 2, JMAX
            ISTRTJ = ISTRTJ + ISTEPJ
            ISTRTV = ISTRTJ
            ISTRTU = ISTRTV
            IPREV1 = ISTRTU - ISTEPJ - 1
            IPREV2 = IPREV1 - 1
C
C           RJ(0,0,0)
C
            AHGTF(ISTRTU) = FJW(JMAX - JVAL)
C
C           RJ(1,0,0)
C
            AHGTF(ISTRTU + 1) = - CPX*AHGTF(IPREV1 + 1)
C
C           RJ(t,0,0) for  t > 1
C
            TMIN1 = D1
            DO 300 IT = 2, JVAL
               AHGTF(ISTRTU + IT) = - CPX*AHGTF(IPREV1 + IT)
     *                            + TMIN1*AHGTF(IPREV2 + IT)
               TMIN1 = TMIN1 + D1
  300       CONTINUE
C
C           RJ(t,1,0)
C
            ISTRTU = ISTRTU + ISTEPU
            IPREV1 = IPREV1 + 1
            DO 400 IT = 0, JVAL - 1
               AHGTF(ISTRTU + IT) = - CPY*AHGTF(IPREV1 + IT)
  400       CONTINUE
C
C           RJ(t,u,0) for  u > 1
C
            UMIN1 = D1
            DO 500 IU = 2,JVAL
               ISTRTU = ISTRTU + ISTEPU
               IPREV1 = IPREV1 + ISTEPU
               IPREV2 = IPREV1 - ISTEPU
               DO 510 IT = 0, JVAL - IU
                  AHGTF(ISTRTU + IT) = - CPY*AHGTF(IPREV1 + IT)
     *                               + UMIN1*AHGTF(IPREV2 + IT)
  510          CONTINUE
               UMIN1 = UMIN1 + D1
  500       CONTINUE
C
C           RJ(t,u,1)
C
            ISTRTV = ISTRTV + ISTEPV
            ISTRTU = ISTRTV
            IPREV1 = ISTRTU - ISTEPJ - ISTEPV
            IUMAX = JVAL - 1
            DO 600 IU = 0, IUMAX
               DO 610 IT = 0, IUMAX - IU
                  AHGTF(ISTRTU + IT) = - CPZ*AHGTF(IPREV1 + IT)
  610          CONTINUE
               ISTRTU = ISTRTU + ISTEPU
               IPREV1 = IPREV1 + ISTEPU
  600       CONTINUE
C
C           RJ(t,u,v) for v > 1
C
            VMIN1 = D1
            DO 700 IV = 2,JVAL
               ISTRTV = ISTRTV + ISTEPV
               ISTRTU = ISTRTV
               IPREV1 = ISTRTU - ISTEPJ - ISTEPV
               IPREV2 = IPREV1 - ISTEPV
               IUMAX = JVAL - IV
               DO 710 IU = 0, IUMAX
                  DO 720 IT = 0, IUMAX - IU
                     AHGTF(ISTRTU + IT) = - CPZ*AHGTF(IPREV1 + IT)
     *                                  + VMIN1*AHGTF(IPREV2 + IT)
  720             CONTINUE
                  ISTRTU = ISTRTU + ISTEPU
                  IPREV1 = IPREV1 + ISTEPU
                  IPREV2 = IPREV2 + ISTEPU
  710          CONTINUE
               VMIN1 = VMIN1 + D1
  700       CONTINUE
            ISTEPJ = - ISTEPJ
  200    CONTINUE
      END IF    
            
      RETURN
      END
      
      
      SUBROUTINE CINT0_Tx(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,STDER0
     &                 ,SHGTF,ADER,AHGTF,NATOMC,IPRINT)
     
        use turin

C
C     TUH
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"

      PARAMETER (D0 = 0.0 D00, DP5 = 0.5 D00)
      DIMENSION STDER0(KCKTAB,2), AHGTF(*), ADER(*)
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)

#include "onecom.h"
#include "lmns.h"
      SFAC = SHGTF**3
      TFAC = -DP5*SFAC
      INT = 0
      DO 100 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 100 ICOMPB = 1,KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         
         
C     !!!!!!!!!!!!!!!!!! adding the integrals !!!!!!!!!!!!!!!!!
C     !          STDER0 1 is the overlap integral             !
C     !          STDER0 2 is the kinetic integral             !
C     !          ADER     is the nuclear attraction integral  !                
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
         
C
C    **********************************************************
C    ***** CALCULATE OVERLAP AND KINETIC ENERGY INTEGRALS *****
C    **********************************************************
C
         X0 = ODC(LVALA,LVALB,0,0,0,1)
         Y0 = ODC(MVALA,MVALB,0,0,0,2)
         Z0 = ODC(NVALA,NVALB,0,0,0,3)
         X2 = ODC(LVALA,LVALB,0,2,0,1)
         Y2 = ODC(MVALA,MVALB,0,2,0,2)
         Z2 = ODC(NVALA,NVALB,0,2,0,3)
         
         INT = INT + 1
         STDER0(INT,1) = STDER0(INT,1)+SFAC*X0*Y0*Z0
         STDER0(INT,2) = STDER0(INT,2)+TFAC*(X2*Y0*Z0+X0*Y2*Z0+X0*Y0*Z2)
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C
         IADRAV = 1
         AINT = D0
         DO 200 IV = 0, NVALA + NVALB
            EV = ODC(NVALA,NVALB,IV,0,0,3)
            IADRAU = IADRAV
            DO 300 IU = 0, MVALA + MVALB
               EE = ODC(MVALA,MVALB,IU,0,0,2)*EV
               DO 400 IT = 0, LVALA + LVALB
                  EEE = ODC(LVALA,LVALB,IT,0,0,1)*EE
                  IADR00 = IADRAU + IT
                  IADD = - NAHGTF
                  DO 500 IATOM = 1,NATOMC
                     IADD = IADD + NAHGTF
                     AINT = AINT + EEE*AHGTF(IADR00 + IADD)
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         ADER(INT) = ADER(INT) + AINT
    
  100 CONTINUE
      RETURN
      END