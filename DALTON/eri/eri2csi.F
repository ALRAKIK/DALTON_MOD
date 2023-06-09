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
!
C
C  /* Deck ericsm */
      SUBROUTINE ERICSM(SO,CSMAT,IPNTCR,INDXBT,WORK,LWORK,IPRINT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "eridst.h"
#include "eribuf.h"
      DIMENSION SO(NAOINT,*), CSMAT(*), IPNTCR(MAXBCH,4), INDXBT(*),
     &          WORK(LWORK) 
#include "cbieri.h"
#include "ericom.h"
#include "symmet.h"
#include "hertop.h"
#include "inftap.h"
C
C     Allocation for ERIOUT
C
      LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
      KNDORB = 1 
      KLAST  = KNDORB + (2*LBIN - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIOUT',' ',KLAST,LWORK)
      CALL ERICS1(SO,CSMAT,WORK(KNDORB),IPNTCR,INDXBT,IPRINT)
      RETURN
      END
C  /* Deck erics1 */
      SUBROUTINE ERICS1(SO,CSMAT,INDORB,IPNTCR,INDXBT,IPRINT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER (D0 = 0.0D0, TEN14 = 1.1D14)
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,*),
     &          INDORB(NCCS,2), IPNTCR(MAXBCH,4),
     &          IPNRST(0:7,3), INDXBT(MXSHEL*MXCONT,0:7),
     &          IADCMP(MXAQN,MXAQN,2), CSMAT(NBAST,NBAST)
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "aobtch.h"
#include "hertop.h"
#include "inforb.h"
#include "symmet.h"
C

      IBTEST(I,J,K,L) = IAND(I,IEOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERICSM',-1)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      NBITX = 16
      IBITA = 2**NBITX
      IBITB = 1
      IBITC = 2**NBITX
      IBITD = 1
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
         WRITE (LUPRI,'(2X,A,4I5)') 
     &                    ' NCTFX         ',NCTFA,NCTFB,NCTFC,NCTFD 
         WRITE (LUPRI,'(2X,A,2I5)') ' NCCS2  ',NCCS2
      END IF
C
C
      IF (NCCS2 .GT. 0) THEN
C
         DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            C = A
            D = B
            CD = IEOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 200
            IF (DIAGCD .AND. D.GT.C) GO TO 200
C
            IADR = 0
            DO N = 1, NPQBCS
               NA = KNDXBT(IPNTCR(N,1)) - 1
               NB = KNDXBT(IPNTCR(N,2)) - 1
               NC = KNDXBT(IPNTCR(N,3)) - 1
               ND = KNDXBT(IPNTCR(N,4)) - 1
               DO I = 1, NCTFA 
               DO J = 1, NCTFB 
               DO K = 1, NCTFC 
               DO L = 1, NCTFD 
                  IADR = IADR + 1
                  INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                           + (INDXBT(NB + J,B)-1)
                  INDORB(IADR,2) = (INDXBT(NC + K,C)-1)*IBITC
     &                           + (INDXBT(ND + L,D)-1)
               END DO 
               END DO 
               END DO 
               END DO
            END DO
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
C
            IA   = 0
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB = IAB + IBITB
                  IC = 0
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  DO ICMPC = 1, ICMPA
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC  = IC + IBITC
                     ICD = IC
                     IF (ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ICD  = ICD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IF (ICMPAB .EQ. ICMPCD) THEN
                           DO I = 1, NCCS
                              KAB = INDORB(I,1) + IAB
                              KCD = INDORB(I,2) + ICD
                              IF (KAB.EQ.KCD) THEN
                                 KA = IAND(ISHFT(KAB,-NBITX),IBITA-1)
                                 KB = IAND(       KAB,       IBITA-1)
                                 CSI = SQRT(SO(I,R,S,T,ICMPAB,ICMPCD,1))
                                 CSMAT(KA,KB) = CSI
                                 CSMAT(KB,KA) = CSI
                              END IF
                           END DO
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  200       CONTINUE
C
         END IF
         END DO
         END IF
         END DO
      END IF
C
      RETURN
      END
C  /* Deck ericso */
      SUBROUTINE ERICSO(SO,ITYPE,INDORB,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,INT,LBIN,NIBUF,INDXBT,NBITS,IPRINT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      PARAMETER (D0 = 0.0D0, TEN14 = 1.1D14)
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,*),
     &          INDORB(NCCS,NIBUF),
     &          IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          BIN(LBIN), IBIN(LBIN,NIBUF),
     &          IPNRST(0:7,3), INDXBT(MXSHEL*MXCONT,0:7),
     &          IADCMP(MXAQN,MXAQN,2)
#include "cbieri.h"
#include "ericom.h"
#include "erithr.h"
#include "aobtch.h"
#include "hertop.h"
#include "symmet.h"
C

      IBTEST(I,J,K,L) = IAND(I,IEOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERICSO',-1)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      NBITX = 16
      IBITA = 2**NBITX
      IBITB = 1
      IBITC = 2**NBITX
      IBITD = 1
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
         WRITE (LUPRI,'(2X,A,4I5)') 
     &                    ' NCTFX         ',NCTFA,NCTFB,NCTFC,NCTFD 
         WRITE (LUPRI,'(2X,A,2I5)') ' NCCS2  ',NCCS2
         WRITE (LUPRI,'(2X,A,2I5)') ' NIBUF, NBITX  ',NIBUF, NBITX 
      END IF
C
      INT = 0
C
      IF (NCCS2 .GT. 0) THEN
C
         DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            C = A
            D = B
            CD = IEOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 200
            IF (DIAGCD .AND. D.GT.C) GO TO 200
C
            IADR = 0
            DO N = 1, NPQBCS
               NA = KNDXBT(IPNTCR(N,1)) - 1
               NB = KNDXBT(IPNTCR(N,2)) - 1
               NC = KNDXBT(IPNTCR(N,3)) - 1
               ND = KNDXBT(IPNTCR(N,4)) - 1
               DO I = 1, NCTFA 
               DO J = 1, NCTFB 
               DO K = 1, NCTFC 
               DO L = 1, NCTFD 
                  IADR = IADR + 1
                  INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                           + (INDXBT(NB + J,B)-1)
                  INDORB(IADR,2) = (INDXBT(NC + K,C)-1)*IBITC
     &                           + (INDXBT(ND + L,D)-1)
               END DO 
               END DO 
               END DO 
               END DO
            END DO
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
C
            IA   = 0
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB = IAB + IBITB
                  IC = 0
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  DO ICMPC = 1, ICMPA
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC  = IC + IBITC
                     ICD = IC
                     IF (ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ICD  = ICD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IF (ICMPAB .EQ. ICMPCD) THEN
                           DO I = 1, NCCS
                              KAB = INDORB(I,1) + IAB
                              KCD = INDORB(I,2) + ICD
C                             KA = IAND(ISHFT(KAB,-NBITX),IBITA-1)
C                             KB = IAND(       KAB,       IBITA-1)
C                             KC = IAND(ISHFT(KCD,-NBITX),IBITA-1)
C                             KD = IAND(       KCD,       IBITA-1)
C                             IJ = MAX(KA,KB)*IBITA + KA + KB
C                             KL = MAX(KC,KD)*IBITA + KC + KD
C                             IF (IJ.EQ.KL) THEN
                              IF (KAB.EQ.KCD) THEN
                                 CSI = SQRT(SO(I,R,S,T,ICMPAB,ICMPCD,1))
                                 IF (ABS(CSI) .GT. THRSH) THEN
                                    INT = INT + 1
                                    BIN (INT  ) = CSI 
                                    IBIN(INT,1) = KAB 
                                 END IF
                              END IF
                           END DO
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  200       CONTINUE
C
         END IF
         END DO
         END IF
         END DO
      END IF
C
      RETURN
      END
C  /* Deck csmext */
      SUBROUTINE CSMEXT(CSMAT)
C
C     October 2004 tuh
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.D0)
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      DIMENSION CSMAT(NBAST,NBAST), CSMSHL(KMAX,KMAX)
#include "inforb.h"
#include "shells.h"
      IOFF = 0
      DO ISHELL = 1, KMAX
         JOFF = 0
         DO JSHELL = 1, KMAX
            CSMAX = D0
            DO ICMP = IOFF + 1, IOFF + KHKT(ISHELL)
               DO JCMP = JOFF + 1, JOFF + KHKT(JSHELL)
                  CSMAX = MAX(CSMAT(ICMP,JCMP),CSMAX)
               END DO
            END DO
            CSMSHL(ISHELL,JSHELL) = CSMAX
            JOFF = JOFF + KHKT(JSHELL)
         END DO
         IOFF = IOFF + KHKT(ISHELL)
      END DO
      IF (.FALSE.) THEN
         CALL HEADER('Subroutine CSMAT',-1) 
         CALL OUTPUT(CSMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Subroutine CSMSHL',-1)
         CALL OUTPUT(CSMSHL,1,KMAX,1,KMAX,KMAX,KMAX,1,LUPRI)
      END IF
      RETURN
      END
