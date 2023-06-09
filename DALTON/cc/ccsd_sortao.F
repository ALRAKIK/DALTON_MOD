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
C  /* Deck ccsd_sortao */
      SUBROUTINE CCSD_SORTAO(WORK,LWORK)
C
C     Written by Henrik Koch 25-Sep-1993
C
C     Purpose: Read destribution of AO integrals.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
C
      LOGICAL FIRST, ENABLED
      DIMENSION WORK(LWORK)
C
      CHARACTER*8 NAME(8)
      CHARACTER*8 LBLSAV
C
C     SAVE FIRST
C
      DATA FIRST /.TRUE./ 
      DATA NAME  /'CCAOIN_1','CCAOIN_2','CCAOIN_3','CCAOIN_4',
     *            'CCAOIN_5','CCAOIN_6','CCAOIN_7','CCAOIN_8'/
      COMMON/SORTIO/LUAOIN(8)
#include "inftap.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccfop.h"
#include "ccsdinp.h"
#include "ccpack.h"
#include "r12int.h"
C
C--------------------
C     Only sort once.
C-------------------- 
C
C     IF (FIRST) THEN
C        FIRST = .FALSE.
C     ELSE
C        RETURN
C     ENDIF
C
C-------------------------
C     Open integral files.
C-------------------------
C
      IF (LUINTA.LE.0) THEN
        LUINTA = -1
        CALL MAKE_AOTWOINT(WORK,LWORK)
        CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
      END IF
C
      DO 50 ISYM = 1,NSYM
C
         NFILE = 0
         CALL WOPEN2(NFILE,NAME(ISYM),64,0)
         LUAOIN(ISYM) = NFILE
C
   50 CONTINUE
C
C------------------------------
C     Skip sorting if required.
C------------------------------
C     
      IF (NOSORT) THEN
C     
          DO ISYMD = 1,NSYM
             ISYDIS = MULD2H(ISYMD,ISYMOP)
             LENGTH = NDISAO(ISYDIS)
             IOFFID = 1
             DO ID = 1,NBAS(ISYMD)
                IDEL = ID + IBAS(ISYMD)
                IOFFINT(IDEL) = IOFFID
                IOFFID = IOFFID + LENGTH
             END DO
          END DO
C
          RETURN
C     
      END IF
C
C------------------------
C     Buffer information.
C------------------------
C
      IF (CCR12) THEN
        NALLBAS = 0
        DO I = 1, NSYM
          NALLBAS = NALLBAS + MBAS1(I) + MBAS2(I)
        END DO
      ELSE
        NALLBAS = NBAST 
      ENDIF
C
      LBUF = 600
      IF (NALLBAS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
         IBIT1 = 2**8  - 1
         IBIT2 = 2**16 - 1
      ELSE
         NIBUF = 2
         NBITS = 16
         IBIT1 = 2**16 - 1
         IBIT2 = 0   ! not used when NIBUF .eq. 2
      END IF
C
C-----------------------
C     Buffer allocation.
C-----------------------
C
      KRBUF = 1
      KIBUF = KRBUF + LBUF
      KAOAB = KIBUF + (NIBUF*LBUF + 1)/2 + 1 ! IBUF always integer*4
      KAOG  = KAOAB + (N2BASX     + 1)/IRAT + 1
      KEND1 = KAOG  + (NBAST*NSYM + 1)/IRAT + 1
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) CALL QUIT('Insufficient work space in CCRDAO')
C
C------------------------------------------------------
C     Calculate in the index arrays needed in the sort.
C------------------------------------------------------
C
      CALL CCSD_INIT2(WORK(KAOAB),WORK(KAOG))
C
C-------------------------------------------
C     set up table for packing of integrals:
C-------------------------------------------
C
      IF (LPACKINT) THEN
         DTIME = SECOND()
         CALL INITPCKR8(THRPCKINT,IPCKTABINT,ENABLED)

         IF (.NOT.ENABLED) THEN
           WRITE(LUPRI,'(A)') 
     &     'packing routines not enabled for this architecture...',
     &     '...the integral packing is switched off...'
           LPACKINT = .FALSE.
         END IF

         NTOTPCK  = 0
         NTOTINT  = 0
         PCKRATIO = 1.0D0
         PCKTIME   = SECOND() - DTIME
      END IF
C
C------------------------------------
C     Loop over batches of integrals.
C------------------------------------
C
C
      DO 100 ISYMD = 1,NSYM
C
         IOFF2 = 1
C
         ISYDIS = MULD2H(ISYMD,ISYMOP)
         NTOTD  = NBAS(ISYMD)
       IF (NTOTD .EQ. 0) GOTO 100
       LENGTH = NDISAO(ISYDIS)
C
         NUMBAT = MIN(NTOTD,LWRK1/LENGTH)
C
         IF (NUMBAT .EQ. 0) THEN
            WRITE(LUPRI,*) 'In CCSD_SORTAO NUMBAT is zero'
            CALL QUIT('Insufficient work space in CCRDAO')
         ENDIF
C
         ITOTBA = (NTOTD-1)/NUMBAT + 1
C
         ID1   = IBAS(ISYMD) + 1
         ID2   = IBAS(ISYMD)
         IOFF1 = IBAS(ISYMD)
C
         DO 200 I = 1,ITOTBA
C
            INUMBA = NUMBAT
            IF (NUMBAT*I .GT. NTOTD) THEN
               INUMBA = NTOTD - NUMBAT*(I-1)
            ENDIF
C
            ID2 = ID2 + INUMBA
C
            CALL DZERO(WORK(KEND1),LENGTH*INUMBA)
C
            CALL CCSD_SORT1(LUINTA,WORK(KEND1),WORK(KIBUF),WORK(KRBUF),
     *                      WORK(KAOAB),WORK(KAOG),ISYDIS,LENGTH,
     *                      IOFF1,ID1,ID2,NIBUF,LBUF,NBITS,IBIT1)
C
            CALL CCSD_SORT2(WORK(KEND1),IOFF2,INUMBA,LENGTH,
     *                      ID1,ID2,ISYMD)
C
            IF (IPRINT .GT. 50) THEN
               CALL AROUND('Integral distribution')
               IPRC = KEND1
               DO 210 IPRD = ID1,ID2
                  WRITE(LUPRI,*) 'D distribution',IPRD
                  DO 220 IPSYMG = 1,NSYM
                     WRITE(LUPRI,*) 'Gamma symmetry',IPSYMG
                     ISYMAB = MULD2H(IPSYMG,ISYDIS)
                     CALL OUTPUT(WORK(IPRC),1,NNBST(ISYMAB),1,
     *                           NBAS(IPSYMG),NNBST(ISYMAB),
     *                           NBAS(IPSYMG),1,LUPRI)
                     IPRC = IPRC + NNBST(ISYMAB)*NBAS(IPSYMG)
  220             CONTINUE
  210          CONTINUE
            END IF
C
            ID1   = ID1   + INUMBA
            IOFF1 = IOFF1 + INUMBA
C
  200    CONTINUE
C
  100 CONTINUE
C
C-------------------------------------
C     Print packing statistics:
C-------------------------------------
C
      IF (IPRINT.GT.0 .AND. LPACKINT) THEN
         WRITE (LUPRI,'(//10X,A,F9.2,A)')
     &        'Time needed to pack integrals:   ',
     &                   PCKTIME, ' seconds'
         WRITE (LUPRI,'(10X,A,G9.2)')
     &        'Threshold used for packing:      ',
     &                   THRPCKINT
         NTOTINT  = MAX(NTOTINT,1)
         PCKRATIO = DBLE(NTOTPCK)/DBLE(NTOTINT)
         WRITE (LUPRI,'(10X,A,F9.2,A)')
     &        'Reduction obtained by packing:   ',
     &      100.0D0*(1.0D0 - PCKRATIO),' %'
      END IF
C
C-------------------------------------
C     Close integral files and delete.
C-------------------------------------
C
      IF (KEEPAOTWO .EQ. 0) THEN
         CALL GPCLOSE(LUINTA,'DELETE')
      ELSE
         CALL GPCLOSE(LUINTA,'KEEP')
      ENDIF
C
      RETURN
      END
C  /* Deck ccsd_sort1 */
      SUBROUTINE CCSD_SORT1(LUINTA,XINT,IBUF4,RBUF,KAOAB,KAOG,ISYDIS,
     *                      LENGTH,IOFF,ID1,ID2,NIBUF,LBUF,NBITS,IBIT1)
C
C     Written by Henrik Koch 25-Sep-1993
C
#include "implicit.h"
#include "priunit.h"
#include "ibtpar.h"
#include "ccorb.h"
      REAL*8    XINT(LENGTH,*),RBUF(LBUF)
      INTEGER*4 IBUF4(NIBUF*LBUF), LENGTH4
      INTEGER   KAOAB(NBAST,NBAST),KAOG(NBAST,NSYM)
#include "ccsdsym.h"
#include "r12int.h"

C
C     INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      REWIND (LUINTA)
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
C
      IF (NIBUF .EQ. 1) THEN
C
   10    READ(LUINTA,ERR=2000) RBUF,IBUF4,LENGTH4
C
         IF (LENGTH4 .EQ. -1) GOTO 100
C
         DO I = 1,LENGTH4
C
            LABLE = IBUF4(I)
            VALUE = RBUF(I)
C
            IP = IAND(       LABLE         ,IBIT1)
            IQ = IAND(ISHFT(LABLE,  -NBITS),IBIT1)
            IR = IAND(ISHFT(LABLE,-2*NBITS),IBIT1)
            IS = IAND(ISHFT(LABLE,-3*NBITS),IBIT1)
            IF (NOAUXB) CALL IJKAUX(IP,IQ,IR,IS)
C
            IF ((IS .GE. ID1) .AND. (IS .LE. ID2)) THEN
               IADR = KAOG(IR,ISYDIS) + KAOAB(IP,IQ)
               XINT(IADR,IS-IOFF) = VALUE
            ENDIF
            IF ((IR .GE. ID1) .AND. (IR .LE. ID2)) THEN
               IADR = KAOG(IS,ISYDIS) + KAOAB(IP,IQ)
               XINT(IADR,IR-IOFF) = VALUE
            ENDIF
            IF ((IP .GE. ID1) .AND. (IP .LE. ID2)) THEN
               IADR = KAOG(IQ,ISYDIS) + KAOAB(IR,IS)
               XINT(IADR,IP-IOFF) = VALUE
            ENDIF
            IF ((IQ .GE. ID1) .AND. (IQ .LE. ID2)) THEN
               IADR = KAOG(IP,ISYDIS) + KAOAB(IR,IS)
               XINT(IADR,IQ-IOFF) = VALUE
            ENDIF
C
         END DO
C
         GOTO 10
  100    CONTINUE
C
      ELSE
C
   30    READ(LUINTA,ERR=2000) RBUF,IBUF4,LENGTH4
C
         IF (LENGTH4 .EQ. -1) GOTO 200
C
         DO 40 I = 1,LENGTH4
C
            LABLE1 = IBUF4(2*I-1)
            LABLE2 = IBUF4(2*I  )
            VALUE = RBUF(I)
C
            IP = IAND(       LABLE1       ,IBIT1)
            IQ = IAND(ISHFT(LABLE1,-NBITS),IBIT1)
            IR = IAND(       LABLE2,       IBIT1)
            IS = IAND(ISHFT(LABLE2,-NBITS),IBIT1)
            IF (NOAUXB) CALL IJKAUX(IP,IQ,IR,IS)
C
            IF ((IS .GE. ID1) .AND. (IS .LE. ID2)) THEN
               IADR = KAOG(IR,ISYDIS) + KAOAB(IP,IQ)
               XINT(IADR,IS-IOFF) = VALUE
            ENDIF
            IF ((IR .GE. ID1) .AND. (IR .LE. ID2)) THEN
               IADR = KAOG(IS,ISYDIS) + KAOAB(IP,IQ)
               XINT(IADR,IR-IOFF) = VALUE
            ENDIF
            IF ((IP .GE. ID1) .AND. (IP .LE. ID2)) THEN
               IADR = KAOG(IQ,ISYDIS) + KAOAB(IR,IS)
               XINT(IADR,IP-IOFF) = VALUE
            ENDIF
            IF ((IQ .GE. ID1) .AND. (IQ .LE. ID2)) THEN
               IADR = KAOG(IP,ISYDIS) + KAOAB(IR,IS)
               XINT(IADR,IQ-IOFF) = VALUE
            ENDIF
C
   40    CONTINUE
C
         GOTO 30
  200    CONTINUE
C
      ENDIF
C
      RETURN
 2000 CALL QUIT('Error reading AOTWOINT in CCSD_SORT1')
      END
C  /* Deck ccsd_sort2 */
      SUBROUTINE CCSD_SORT2(XINT,IOFF,INUMBA,LENGTH,ID1,ID2,ISYM)
C
C     Written by Henrik Koch 25-Sep-1993
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION XINT(*)
C
      CHARACTER*8 NAME(8)
C
      DATA NAME  /'CCAOIN_1','CCAOIN_2','CCAOIN_3','CCAOIN_4',
     *            'CCAOIN_5','CCAOIN_6','CCAOIN_7','CCAOIN_8'/

      COMMON/SORTIO/LUAOIN(8)
#include "ccorb.h"
#include "maxorb.h"
#include "ccpack.h"
C
      NFILE = LUAOIN(ISYM)
      KOFF  = 1
C
      DO IDEL = ID1, ID2
C
         NBYTE   = LENGTH*8
         NDWORDS = LENGTH
C
         IF (LPACKINT) THEN

            DTIME = SECOND()

            CALL PCKR8(XINT(KOFF),LENGTH,XINT(KOFF),NBYTE,
     &                 IPCKTABINT,LPACKINT)

            NDWORDS = (NBYTE+7)/8
            NTOTINT = NTOTINT + LENGTH  / NBAST
            NTOTPCK = NTOTPCK + NDWORDS / NBAST
            PCKTIME = PCKTIME  + SECOND() - DTIME

         END IF
C
         CALL PUTWA2(NFILE,NAME(ISYM),XINT(KOFF),IOFF,NDWORDS)
C
         IOFFINT(IDEL) = IOFF
         NPCKINT(IDEL) = NBYTE
C
         IOFF = IOFF + NDWORDS
         KOFF = KOFF + LENGTH
C
      END DO
C
      RETURN
      END
C  /* Deck ccsd_init2 */
      SUBROUTINE CCSD_INIT2(KAOAB,KAOG)
C
C     Henrik Koch and Alfredo Sanchez.       29-Jun-1994
C
C     Set up indexing arrays
C
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
      DIMENSION KAOAB(NBAST,NBAST),KAOG(NBAST,NSYM)
#include "ccsdsym.h"
C
C
      DO 100 ISYMAB = 1,NSYM
         ICOUNT = 0
         DO 110 ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYMAB)
            IF (ISYMB .GT. ISYMA) THEN
               DO 120 B = 1,NBAS(ISYMB)
                  IB = IBAS(ISYMB) + B
                  DO 130 A = 1,NBAS(ISYMA)
                     IA = IBAS(ISYMA) + A
                     ICOUNT = ICOUNT + 1
                     KAOAB(IA,IB) = ICOUNT
                     KAOAB(IB,IA) = ICOUNT
  130             CONTINUE
  120          CONTINUE
            ELSE IF (ISYMA .EQ. ISYMB) THEN
               DO 140 B = 1,NBAS(ISYMB)
                  IB = IBAS(ISYMB) + B
                  DO 150 A = 1,B
                     IA = IBAS(ISYMA) + A
                     ICOUNT = ICOUNT + 1
                     KAOAB(IA,IB) = ICOUNT
                     KAOAB(IB,IA) = ICOUNT
  150             CONTINUE
  140          CONTINUE
            END IF
  110    CONTINUE
  100 CONTINUE
C
      DO 200 ISYMD = 1,NSYM
         ISYDIS = MULD2H(ISYMD,ISYMOP)
         ICOUNT = 0
         DO 210 ISYMG = 1,NSYM
            ISYMAB = MULD2H(ISYMG,ISYDIS)
            DO 220 G = 1,NBAS(ISYMG)
               IG = IBAS(ISYMG) + G
               KAOG(IG,ISYMD) = ICOUNT
               ICOUNT = ICOUNT + NNBST(ISYMAB)
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
