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
C FILE: abacus/aba2r12.F
C
C  /* Deck eriher */
      SUBROUTINE ABEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,
     &                  EXP12,EXP34,NUC12,NUC34,NUABCD,
     &                  NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
C
C     Written by Wim Klopper (University of Karlsruhe, 14 November 2002).
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DM2 = -2.0D0)
      DIMENSION Q(*), R(*), W(*), A(*), EXP12(NUC12), EXP34(NUC34),
     &          IODDHR(*), INDHER(*), PQX(*), PQY(*), PQZ(*)
      IODS = 0
      DO  IOD12 = 1, NUC12
         EXPP = EXP12(IOD12)
         DO  IOD34 = 1, NUC34
            IODS = IODS + 1
            EXPQ = EXP34(IOD34)
            W(IODS) = DM2*EXPP*EXPQ/(EXPP + EXPQ)
         END DO
      END DO
      CALL WKEQ52(Q,R,W,A,PQX,PQY,PQZ,INDHER,JMAX,NUABCD,
     &            NTUV,IPQ0X,IPQ0Y,IPQ0Z,IODDHR,IPRINT)
      RETURN
      END
C  /* Deck r00G */
      SUBROUTINE R00G(RJ000,COOR12,COOR34,EXP12,EXP34,FAC12,FAC34,PQX,
     &               PQY,PQZ,JMAX,NOINT,NUABCD,NUC1,NUC2,NUC12,NUC3,
     &               NUC4,NUC34,THRESH,ONECEN,IPRINT,IPQ0X,IPQ0Y,IPQ0Z,
     &               SIGNT,FACINT,HEXPP,HEXPQ)
C     Copy of R0001 for use with R12EIN (WK/UniKA/20-11-2002).
#include "implicit.h"
#include "priunit.h"
#include "subdir.h"
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, DP25 = 0.25D0)
      LOGICAL ONECEN, NOINT
      DIMENSION RJ000(NUABCD,0:JMAX),
     &          PQX(NUABCD), PQY(NUABCD), PQZ(NUABCD),
     &          COOR12(NUC1*NUC2,3), COOR34(NUC3*NUC4,3),
     &          EXP12(*), EXP34(*), FAC12 (*), FAC34(*), SIGNT(3),
     &          FACINT(*), HEXPP(*) ,HEXPQ(*)
      NOINT = .FALSE.
      IF (ONECEN) THEN
         CALL DZERO(PQX,NUABCD)
         CALL DZERO(PQY,NUABCD)
         CALL DZERO(PQZ,NUABCD)
         IPQ0X = 1
         IPQ0Y = 1
         IPQ0Z = 1
         IODS = NUABCD
         NODS = NUABCD
      ELSE
         IF (.NOT.DPATH1) THEN
            SGN12X = - SIGNT(1)
            SGN12Y = - SIGNT(2)
            SGN12Z = - SIGNT(3)
            SGN34X = - D1
            SGN34Y = - D1
            SGN34Z = - D1
         ELSE
            SGN12X = D1
            SGN12Y = D1
            SGN12Z = D1
            SGN34X = SIGNT(1)
            SGN34Y = SIGNT(2)
            SGN34Z = SIGNT(3)
         END IF
C
         IODS  = 1
         NODS  = 1
         DO 300 IOD12 = 1, NUC12
            PX     = SGN12X*COOR12(IOD12,1)
            PY     = SGN12Y*COOR12(IOD12,2)
            PZ     = SGN12Z*COOR12(IOD12,3)
            DO 310 IOD34 = 1, NUC34
               PQXI = PX - SGN34X*COOR34(IOD34,1)
               PQYI = PY - SGN34Y*COOR34(IOD34,2)
               PQZI = PZ - SGN34Z*COOR34(IOD34,3)
               PQX(IODS) = PQXI
               PQY(IODS) = PQYI
               PQZ(IODS) = PQZI
               IODS = IODS + 1
               NODS = NODS + 1
  310       CONTINUE
  300    CONTINUE

         IPQ0X = 1
         IPQ0Y = 1
         IPQ0Z = 1
         IF (DASUM(NUABCD,PQX,1) .GT. THRESH) IPQ0X = 0
         IF (DASUM(NUABCD,PQY,1) .GT. THRESH) IPQ0Y = 0
         IF (DASUM(NUABCD,PQZ,1) .GT. THRESH) IPQ0Z = 0
      END IF
C
      IJ = 0
      DO IOD12 = 1, NUC12
         FAC = FAC12(IOD12)
         DO IOD34 = 1, NUC34
            IJ = IJ + 1
            FACINT(IJ) = FAC * FAC34(IOD34)
            HEXPP(IJ)  = EXP12(IOD12)
            HEXPQ(IJ)  = EXP34(IOD34)
         END DO
      END DO
      RETURN
      END
C  /* Deck r12wrt */
      SUBROUTINE R12WRT(BUF,LBUF,ICOUNT,ITYPE,INDA,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "ibtpar.h"
      DIMENSION BUF(LBUF,5)
#include "twosta.h"
#include "r12int.h"
#include "inftap.h"
#include "nuclei.h"
#include "symmet.h"

C
       IF (ITYPE .EQ. -1) THEN
         IDUM = 0
         I = 1
         CALL GPOPEN(LUR12(I),'AOXYZ','UNKNOWN',' ',' ',IDUM,.FALSE.)
         CALL NEWLAB('AOXYZINT',LUR12(I),LUPRI)
         IF (V12INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOV12','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOV12INT',LUR12(I),LUPRI)
         END IF
         IF (R12INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOR12','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOR12INT',LUR12(I),LUPRI)
         END IF
         IF (U12INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOU12','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOU12INT',LUR12(I),LUPRI)
         END IF
         IF (U21INT) THEN
            I = I + 1
            CALL GPOPEN(LUR12(I),'AOU21','UNKNOWN',' ',' ',IDUM,.FALSE.)
            CALL NEWLAB('AOU21INT',LUR12(I),LUPRI)
         END IF
      END IF
C
      DO 100 KR12 = 1, NOPP12 + 1
       IF (ITYPE .EQ. -1) THEN
         ICOUNT = 0
       ELSE IF (ITYPE .EQ. 0) THEN
         IF (INDA .NE. 0) WRITE (LUR12(KR12)) INDA
         WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),ICOUNT
         IF (KR12 .EQ. NOPP12 + 1) ICOUNT = 0
       ELSE
         IF (INDA .NE. 0) THEN
            WRITE (LUR12(KR12)) -INDA
            WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),ICOUNT
         ELSE
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),ICOUNT
            END IF
            WRITE (LUR12(KR12)) (BUF(L,KR12),L=1,LBUF),-1
         END IF
         CALL GPCLOSE(LUR12(KR12),'KEEP')
       END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck rn2out */
      SUBROUTINE RN2OUT(SO,NSOINT,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                  THRESH,NINDAB,NINDCD,IPRINT)
C
C     Version of UN2OUT for R12 integrals (WK/26-11-2002).
C
#include "implicit.h"
#include "priunit.h"
#include "r12int.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
      PARAMETER (LBUF_alloc = 600)
      LOGICAL DCMPAB, DCMPCD, DCMPAC, DRALTB, DRCLTD, FIRST, LAST,
     &        DRABAB, DCABAB, IAEQIC, IALTIC, IPNTLG(3,*), NOTEST,
     &        GTTHRS
      REAL*8  SO(NSOINT,*), BUF(LBUF_alloc,5),
     &        IPNTNO(4,*), IPNTRP(3,*),
     &        NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2)
#include "nuclei.h"
#include "eribuf.h"
#include "twocom.h"
#include "symmet.h"
      SAVE BUF, XBUF, ICOUNT
C
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine RN2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
      END IF
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) THEN
         LBUF = LBUF_alloc
         CALL R12WRT(BUF,LBUF,ICOUNT,-1,0,IPRINT)
      END IF
C
      ISOFF  = 0
      NBUFCL = 0
      NSTART = ICOUNT
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      DO 100 I = 1, NINTS
         NSTRNA = IPNTNO(1,I)
         NSTRNB = IPNTNO(2,I)
         NSTRNC = IPNTNO(3,I)
         NSTRND = IPNTNO(4,I)
         IREPA  = IPNTRP(1,I)
         IREPB  = IPNTRP(2,I)
         IREPC  = IPNTRP(3,I)
         IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
         IF (NOTEST) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     IC = KHKTC*(NINDCD(ICD,1) - 1)
                     ID = KHKTD*(NINDCD(ICD,2) - 1)
                     INDC = IPTSYM(NSTRNC + IC,IREPC)
                     INDD = IPTSYM(NSTRND + ID,IREPD)
                     CALL LAB64U(SO(ISOFF+INT,1),NSOINT,
     &                           INDA,INDB,INDC,INDD,XABCD,
     &                           THRESH,GTTHRS,IPRINT)
                     IF (GTTHRS) THEN
                        ICOUNT = ICOUNT + 1
                        BUF(ICOUNT,1) = XABCD
                        DO IOPP = 1, NOPP12
                           BUF(ICOUNT,IOPP+1) = SO(ISOFF+INT,IOPP)
                        END DO
                        IF (ICOUNT .EQ. LBUF) THEN
                           NBUFCL = NBUFCL + 1
                           CALL R12WRT(BUF,LBUF,ICOUNT,0,0,IPRINT)
                        END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
         ELSE
            DCMPAB = IPNTLG(1,I)
            DCMPCD = IPNTLG(2,I)
            DCABAB = IPNTLG(3,I)
            DRALTB = IREPA .LT. IREPB
            DRCLTD = IREPC .LT. IREPD
            DRABAB = DCABAB .AND. IREPA.EQ.IREPC .AND. IREPB.EQ.IREPD
            INT = 0
            DO 300 IAB = 1, NORBAB
               IA = KHKTA*(NINDAB(IAB,1) - 1)
               IB = KHKTB*(NINDAB(IAB,2) - 1)
               IF (DCMPAB) THEN
                  IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                     INT = INT + NORBCD
                     GO TO 300
                  END IF
               END IF
               INDA = IPTSYM(NSTRNA + IA,IREPA)
               INDB = IPTSYM(NSTRNB + IB,IREPB)
               DO 310 ICD = 1,NORBCD
                  IC = KHKTC*(NINDCD(ICD,1) - 1)
                  ID = KHKTD*(NINDCD(ICD,2) - 1)
                  INT = INT + 1
                  IF (DCMPCD ) THEN
                     IF (ID.GT.IC) GO TO 310
                     IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  END IF
                  IF (DRABAB) THEN
                     IF (IA.LT.IC.OR.(IA.EQ.IC.AND.IB.LT.ID)) GOTO 310
                  END IF
                  INDC = IPTSYM(NSTRNC + IC,IREPC)
                  INDD = IPTSYM(NSTRND + ID,IREPD)
                  CALL LAB64U(SO(ISOFF+INT,1),NSOINT,
     &                        INDA,INDB,INDC,INDD,XABCD,
     &                        THRESH,GTTHRS,IPRINT)
                  IF (GTTHRS) THEN
                     ICOUNT = ICOUNT + 1
                     BUF(ICOUNT,1) = XABCD
                     DO IOPP = 1, NOPP12
                        BUF(ICOUNT,IOPP+1) = SO(ISOFF+INT,IOPP)
                     END DO
                     IF (ICOUNT .EQ. LBUF) THEN
                        CALL R12WRT(BUF,LBUF,ICOUNT,0,0,IPRINT)
                        NBUFCL = NBUFCL + 1
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NGINT = LBUF*NBUFCL + ICOUNT - NSTART
      CALL DELSTA(0,NGINT)
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) CALL R12WRT(BUF,LBUF,ICOUNT,1,0,IPRINT)
      RETURN
      END
C  /* Deck us2out */
      SUBROUTINE US2OUT(SO,NSOINT,WRKBUF,
     &                  IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                  THRESH,NINDAB,NINDCD,IORBSH,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      LOGICAL FIRST, LAST, IPNTLG(*)
      DIMENSION SO(NSOINT,*), WRKBUF(*), IPNTNO(*), IPNTRP(*),
     &          NINDAB(*), NINDCD(*), IORBSH(*) 
#include "disbuf.h"
C
C---------------------------------
C     Call sort and write routine.
C---------------------------------
C
      CALL US2OU1(SO,NSOINT,WRKBUF(KDSBF),WRKBUF(KDUBF),
     &            WRKBUF(KDSIBF),WRKBUF(KDSNCT),
     &            WRKBUF(KDSORB),WRKBUF(KORBDS),IPNTNO,IPNTRP,IPNTLG,
     &            FIRST,LAST,THRESH,NINDAB,NINDCD,LDSBUF,NDIST,IORBSH,
     &            IPRINT)
C
      RETURN
      END
C  /* Deck us2ou1 */
      SUBROUTINE US2OU1(SO,NSOINT,BUF,CUF,
     &                  IBUF4,NCOUNT,IDSORB,IORBDS,IPNTNO,IPNTRP,
     &                  IPNTLG,FIRST,LAST,THRESH,NINDAB,NINDCD,LDSBUF,
     &                  NDIST,IORBSH,IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "eribuf.h"
#include "nuclei.h"
      LOGICAL DCMPCD, DRCLTD, FIRST, LAST, IPNTLG(3,*), NOTEST,
     &        DOINDX
      DIMENSION SO(NSOINT,*), BUF(LBUF,NDIST), CUF(LBUF,NDIST), 
     &          IPNTNO(4,*), IPNTRP(3,*), NCOUNT(NDIST),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          IDSORB(NDIST), IORBDS(NBASIS), IORBSH(*)
      INTEGER*4 IBUF4(LBUF*NIBUF,NDIST)
#include "twocom.h"
#include "symmet.h"
#include "drw2el.h"
#include "r12int.h"
C
      INTEGER*8 NWRIT
C

C
      IF (LBUF .NE. LDSBUF) THEN
        WRITE (LUPRI,*) 'LBUF .ne. LDSBUF :',LBUF,LDSBUF
        CALL QUIT('Error in US2OU1, LBUF .ne. LDSBUF')
      END IF
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine US2OUT',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
         WRITE (LUPRI,'(2X,A, I5)') 'NDIST ', NDIST
      END IF
C
      IF (NIBUF .EQ. 1) THEN
         NBITS = 8
         IBIT1 = 2**8  - 1
         IBIT2 = 2**16 - 1
      ELSE IF (NIBUF .EQ. 2) THEN
         NBITS = 16
         IBIT1 = 2**16 - 1
         IBIT2 = 0   ! not used when NIBUF .eq. 2
      ELSE
         CALL QUIT('ERROR US2OU1: NIBUF .ne. 1 .and. NIBUF .ne. 2')
      END IF
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) THEN
         CALL IZERO(NCOUNT,NDIST)
         CALL UN2WRU(BUF,CUF,IBUF4,ICOUNT,-1,0,IPRINT)
C        CALL UN2WRT(BUF,IBUF4,LBUF,NIBUF,ICOUNT,-1,NBITS,0,IPRINT)
         DOINDX = .TRUE.
         CALL AINDEX(ISHELA,NAINTS,IDSORB,DOINDX,IORBSH,IPRINT)
         DO 50 IDIST = 1, NDIST
             IORBDS(IDSORB(IDIST)) = IDIST
   50    CONTINUE
         NWRIT = 0
      END IF
C
      ISOFF  = 0
      NBUFCL = 0
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      DO 100 I = 1, NINTS
         NSTRNA = IPNTNO(1,I)
         NSTRNB = IPNTNO(2,I)
         NSTRNC = IPNTNO(3,I)
         NSTRND = IPNTNO(4,I)
         IREPA  = IPNTRP(1,I)
         IREPB  = IPNTRP(2,I)
         IREPC  = IPNTRP(3,I)
         IREPD  = IEOR(IEOR(IREPA,IREPB),IREPC)
         IF (NOTEST) THEN
            IF (NIBUF .EQ. 1) THEN
               INT = 0
               DO 200 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDAB = INDA*(IBIT1 + 1) + IPTSYM(NSTRNB + IB,IREPB)
                  IDIST = IORBDS(INDA)
                  DO 210 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT1 = SO(ISOFF+INT,1)
                     SOINT2 = SO(ISOFF+INT,2)
                     IF (ABS(SOINT1) .GT. THRESH .OR. 
     &                   ABS(SOINT2) .GT. THRESH) THEN
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) = - SOINT2
                       ELSE
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) =   SOINT2
                       END IF
                       IBUF4(ICOUNT,IDIST) = INDAB*(IBIT2 + 1) + INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRU(BUF(1,IDIST),CUF(1,IDIST),
     &                                IBUF4(1,IDIST),
     &                                ICOUNT,0,INDA,IPRINT)
C                         CALL UN2WRT(BUF(1,IDIST),IBUF4(1,IDIST),LBUF,
C    &                                NIBUF,ICOUNT,0,NBITS,INDA,IPRINT)
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
  210             CONTINUE
  200          CONTINUE
            ELSE
               INT = 0
               DO 205 IAB = 1, NORBAB
                  IA = KHKTA*(NINDAB(IAB,1) - 1)
                  IB = KHKTB*(NINDAB(IAB,2) - 1)
                  INDA = IPTSYM(NSTRNA + IA,IREPA)
                  INDB = IPTSYM(NSTRNB + IB,IREPB)
                  INDAB = INDA*(IBIT1 + 1) + INDB
                  IDIST = IORBDS(INDA)
                  DO 215 ICD = 1, NORBCD
                     INT = INT + 1
                     SOINT1 = SO(ISOFF+INT,1)
                     SOINT2 = SO(ISOFF+INT,2)
                     IF (ABS(SOINT1) .GT. THRESH .OR. 
     &                   ABS(SOINT2) .GT. THRESH) THEN
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       IC = KHKTC*(NINDCD(ICD,1) - 1)
                       ID = KHKTD*(NINDCD(ICD,2) - 1)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) = - SOINT2
                       ELSE
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) =   SOINT2
                       END IF
                       IBUF4(2*ICOUNT-1,IDIST) = INDAB
                       IBUF4(2*ICOUNT  ,IDIST) = INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          NBUFCL = NBUFCL + 1
                          CALL UN2WRU(BUF(1,IDIST),CUF(1,IDIST),
     &                                IBUF4(1,IDIST),
     &                                ICOUNT,0,INDA,IPRINT)
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
  215             CONTINUE
  205          CONTINUE
            END IF
         ELSE
            DCMPCD = IPNTLG(2,I)
            DRCLTD = IREPC .LT. IREPD
            INT = 0
            DO 300 IAB = 1, NORBAB
               IA = KHKTA*(NINDAB(IAB,1) - 1)
               IB = KHKTB*(NINDAB(IAB,2) - 1)
               INDA = IPTSYM(NSTRNA + IA,IREPA)
               INDB = IPTSYM(NSTRNB + IB,IREPB)
               INDAB = INDA*(IBIT1 + 1) + INDB
               IDIST = IORBDS(INDA)
               DO 310 ICD = 1,NORBCD
                  IC = KHKTC*(NINDCD(ICD,1) - 1)
                  ID = KHKTD*(NINDCD(ICD,2) - 1)
                  INT = INT + 1
                  IF (DCMPCD ) THEN
                     IF (ID.GT.IC) GO TO 310
                     IF (DRCLTD .AND. ID.EQ.IC) GO TO 310
                  END IF
                  SOINT1 = SO(ISOFF+INT,1)
                  SOINT2 = SO(ISOFF+INT,2)
                  IF (ABS(SOINT1) .GT. THRESH .OR.
     &                ABS(SOINT2) .GT. THRESH) THEN
                     IF (NIBUF .EQ. 1) THEN
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (BPH2OO.AND.INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) = - SOINT2
                       ELSE
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) =   SOINT2
                       END IF
                       IBUF4(ICOUNT,IDIST) = INDAB*(IBIT2 + 1) + INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRU(BUF(1,IDIST),CUF(1,IDIST),
     &                                IBUF4(1,IDIST),
     &                                ICOUNT,0,INDA,IPRINT)
                          NBUFCL = NBUFCL + 1
                          NCOUNT(IDIST) = 0
                       END IF
                     ELSE
                       NCOUNT(IDIST) = NCOUNT(IDIST) + 1
                       ICOUNT = NCOUNT(IDIST)
                       INDC = IPTSYM(NSTRNC + IC,IREPC)
                       INDD = IPTSYM(NSTRND + ID,IREPD)
                       INDCD  = MAX(INDC,INDD)*IBIT1 + INDC + INDD
                       IF (INDD.GT.INDC) THEN 
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) = - SOINT2
                       ELSE
                          BUF (ICOUNT,IDIST) =   SOINT1
                          CUF (ICOUNT,IDIST) =   SOINT2
                       END IF
                       IBUF4(2*ICOUNT-1,IDIST) = INDAB
                       IBUF4(2*ICOUNT  ,IDIST) = INDCD
                       IF (ICOUNT.EQ.LBUF) THEN
                          CALL UN2WRU(BUF(1,IDIST),CUF(1,IDIST),
     &                                IBUF4(1,IDIST),
     &                                ICOUNT,0,INDA,IPRINT)
                          NBUFCL = NBUFCL + 1
                          NCOUNT(IDIST) = 0
                       END IF
                     END IF
                  END IF
  310          CONTINUE
  300       CONTINUE
         END IF
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      NWRIT = NWRIT + LBUF*NBUFCL
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) THEN
         DO 400 IDIST = 1, NDIST
            NWRIT = NWRIT + NCOUNT(IDIST)
            CALL UN2WRU(BUF(1,IDIST),CUF(1,IDIST),IBUF4(1,IDIST),
     &                  NCOUNT(IDIST),1,IDSORB(IDIST),IPRINT)
  400    CONTINUE
         FNALL  = (NBASIS*(NBASIS + 1))/2
         FNALL  = FNALL*NBASIS
         FNALL  = FNALL*NDIST
         PERCNT = NWRIT
         PERCNT = 100.D0*PERCNT / FNALL
         IF (IPRINT.GT.0) WRITE (LUPRI,'(/1X,A,I10,A,F4.1,A)')
     &         'Number of two-electron integrals written:',NWRIT,
     &         ' (',PERCNT,'%)'
      END IF
C
      RETURN
      END
C  /* Deck un2wru */
      SUBROUTINE UN2WRU(BUF,CUF,IBUF4,ICOUNT,ITYPE,INDA,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "ibtpar.h"
#include "eritap.h"
#include "eribuf.h"
      DIMENSION BUF(LBUF), CUF(LBUF)
      INTEGER*4 IBUF4(LBUF*NIBUF), ICOUNT4
#include "drw2el.h"
#include "r12int.h"
#include "twosta.h"
#include "inftap.h"
#include "nuclei.h"
#include "symmet.h"

C
      NBUFX(0) = NBUFX(0) + 1
      ICOUNT4  = ICOUNT
C
      IF (ITYPE .EQ. -1) THEN
         REWIND LUINTR
         NBUFX(0) = NBUFX(0) - 1
c        CALL NEWLAB('BASINFO ',LUINTA,LUPRI)
c        WRITE (LUINTA) MAXREP+1,(NAOS(I),I=1,8),LBUF,NIBUF,NBITS
c        CALL NEWLAB('BASTWOEL',LUINTA,LUPRI)
c        ICOUNT = 0
c        NBUFX(0) = 0
      ELSE IF (ITYPE .EQ. 0) THEN
         IF (INDA .NE. 0) WRITE (LUINTR) INDA
         WRITE (LUAORC(0),REC=NBUFX(0)) BUF,IBUF4,ICOUNT4
         WRITE (LU21INT  ,REC=NBUFX(0)) CUF,IBUF4,ICOUNT4
C
         IF (IPRINT .GE. 6) THEN
            WRITE (LUPRI,'(2X,A,I5,A/)') 'UN2WRU '//
     &         'Integral buffer #',NBUFX(0),' has been written.'
            IBIT1 = 2**NBITS - 1
            DO 100 INT = 1, ICOUNT
               IF (NIBUF .EQ. 1) THEN
                  IJKL = IBUF4(INT) ! IJKL will always be standard integer
                  I = IAND(ISHFT(IJKL,-3*NBITS),IBIT1)
                  J = IAND(ISHFT(IJKL,-2*NBITS),IBIT1)
                  K = IAND(ISHFT(IJKL,  -NBITS),IBIT1)
                  L = IAND(       IJKL,         IBIT1)
               ELSE
                  IJ = IBUF4(2*INT-1)
                  KL = IBUF4(2*INT  )
                  I = IAND(ISHFT(IJ,-NBITS),IBIT1)
                  J = IAND(       IJ,       IBIT1)
                  K = IAND(ISHFT(KL,-NBITS),IBIT1)
                  L = IAND(       KL,       IBIT1)
               END IF
               WRITE (LUPRI,'(10X,A,2X,4I4,5X,2(1P,D16.8))')
     &                      ' ## ', I, J, K, L, BUF(INT), CUF(INT)
  100       CONTINUE
         END IF
         ICOUNT = 0
      ELSE
         IF (INDA .NE. 0) THEN
            WRITE (LUINTR) INDA
            WRITE (LUAORC(0),REC=NBUFX(0)) BUF,IBUF4,ICOUNT4
            WRITE (LU21INT  ,REC=NBUFX(0)) CUF,IBUF4,ICOUNT4
         ELSE
            CALL QUIT('Error in UN2WRU')
         END IF
C
         IF (IPRINT .GE. 6) THEN
            IF (ICOUNT .GT. 0) THEN
               WRITE (LUPRI,'(2X,A,I5,A,I5/)') 'UN2WRU '//
     &            'Integral buffer #',NBUFX(0),' has been written.'//
     &            '   INDA =',INDA
               IBIT1 = 2**NBITS - 1
               DO 200 INT = 1, ICOUNT
                  IF (NIBUF .EQ. 1) THEN
                     IJKL = IBUF4(INT) ! IJKL will always be standard integer
                     I = IAND(ISHFT(IJKL,-3*NBITS),IBIT1)
                     J = IAND(ISHFT(IJKL,-2*NBITS),IBIT1)
                     K = IAND(ISHFT(IJKL,  -NBITS),IBIT1)
                     L = IAND(       IJKL,         IBIT1)
                  ELSE
                     IJ = IBUF4(2*INT-1)
                     KL = IBUF4(2*INT  )
                     I = IAND(ISHFT(IJ,-NBITS),IBIT1)
                     J = IAND(       IJ,       IBIT1)
                     K = IAND(ISHFT(KL,-NBITS),IBIT1)
                     L = IAND(       KL,       IBIT1)
                  END IF
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,2(1P,D16.8))')
     &                         ' ## ', I, J, K, L, BUF(INT), CUF(INT)
  200          CONTINUE
            END IF
         END IF
C
C        Statistics
C
         IF (INDA .EQ. 0) THEN
            N2WRIT = LBUF*NBUFX(0) + ICOUNT
            IF (ICOUNT.GT.0 .AND. INDA.NE.0) THEN
               NBUFX(0) = NBUFX(0) + 2
            ELSE
               NBUFX(0) = NBUFX(0) + 1
            END IF
            IF (IRAT .EQ. 1) LWORD = 8
            IF (IRAT .EQ. 2) LWORD = 4
            FMBYTES = LWORD*(LBUF*IRAT + NIBUF*LBUF + 1)
            FMBYTES = NBUFX(0)*FMBYTES
            FMBYTES = FMBYTES / (1024.D0**2)
            FNALL  = (NBASIS*(NBASIS + 1))/2
            FNALL  = FNALL*(FNALL + 1.0D0)*0.5D0
            PERCNT = N2WRIT
            PERCNT = 100.D0*PERCNT / FNALL
            WRITE (LUPRI,'(/A,I10,A,F4.1,A/A,F10.3//)')
     &         ' Number of two-electron integrals written:',N2WRIT,
     &         ' (',PERCNT,'%)',
     &         ' Megabytes written:                       ',FMBYTES
         END IF
      END IF
      RETURN
      END
C end of FILE: abacus/aba2r12.F
