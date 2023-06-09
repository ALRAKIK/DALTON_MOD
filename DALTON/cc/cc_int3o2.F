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
*----------------------------------------------------------------*
      SUBROUTINE CC_INT3O2(X3OINT,DSRHF,ISYDRHF,
     *                     XLAMDH,ISYMH,XLAMDP,ISYMP,
     *                     WORK,LWORK,IDEL,ISYMD,
     *                     FAC,LWRDSK,LUO3,FILO3,ITRAN,LSQRAB)
*----------------------------------------------------------------*
C   Purpose: To calculate (and write to disc) an integral batch with
C            three occupied indices for a given delta.
C
C    Written by Henrik Koch * Asger Halkier 27/7 - 1995
C
C    Modified by Asger Halkier to return integrals in X3OINT as
C    well as writing them to disc 28/10 - 1995.
C
C    Generalized to also calculated integrals (j del i k) where k is barred
C    thus lamda is not total symmetric.
C    Ove Christiansen 13-6-1996
C
C    Generalized to handle LambdaP and LambdaH matrices of different
C    symmetry. The symmetry of the distribution is passed from
C    outside. 
C    LWRDSK = .FALSE. the result is NOT written on file
C    FAC = 1.0  the integrals are ADDED to previous result of call
C    (second DGEMM).
C    ALL deltas are dumped on file for the given ITRAN
C    to allow later sequencial reading of the whole ITRAN batch
C    of deltas and transformation to occupied using DGEM
C    Generalized to handle a full (al bet| space of DSRHF (LSQRAB = true)
C    Sonia Coriani, 10-03-1999
C     
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION X3OINT(*),DSRHF(*),WORK(LWORK)
      DIMENSION XLAMDH(*),XLAMDP(*)
      CHARACTER FILO3*(*)
      LOGICAL LSQRAB,LWRDSK
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "cclr.h"
C
C-------------------------------------------
C     Work space allocation 1 * outer loops.
C-------------------------------------------
C
      N3ODMX = 0
      DO ISYM  = 1, NSYM
        N3ODMX = MAX(N3ODMX,N3ODEL(ISYM))
      END DO

      ISALBEJ = ISYDRHF
      ISYIKJ = MULD2H(ISALBEJ,MULD2H(ISYMP,ISYMH))
C
      DO 100 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISALBEJ,ISYMJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
C------------------------------------------------------------
C           Work space allocation 1 * unpacking of integrals
C           if LSQRAB = .FALSE.
C------------------------------------------------------------
C
            IF (.NOT.LSQRAB) THEN
              KUNINT = 1
              KEND1  = KUNINT + N2BST(ISALBE)
              LWRK1  = LWORK  - KEND1
C
              IF (LWRK1 .LT. 0) THEN
                 CALL QUIT(
     &              '1-Insufficient work space area in CC_INT3O2')
              ENDIF
C
              KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J-1)+1
C
              CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KUNINT))
            ELSE
              KUNINT = IDSRHFSQ(ISALBE,ISYMJ) + N2BST(ISALBE)*(J-1)+1
              KEND1 = 1
            END IF
C
            DO 120 ISYMK = 1,NSYM
C
C-----------------------------------------------------------------------
C              Transform remaining AO-indices of integrals to occ. space
C-----------------------------------------------------------------------
C
               ISYMBE = MULD2H(ISYMK,ISYMH)
               ISYMAL = MULD2H(ISALBE,ISYMBE)
               ISYMI  = MULD2H(ISYMAL,ISYMP)
               ISYMIK = MULD2H(ISYMI,ISYMK)
C
               KINMA1 = KEND1
               KEND2  = KINMA1 + NBAS(ISYMAL)*NRHF(ISYMK)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT(
     &              '2-Insufficient work space area in CCINT3O2')
               ENDIF
C
               KOFF1 = KUNINT + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMRH(ISYMBE,ISYMK) + 1
               KOFF3 = KINMA1
C
               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTB = MAX(NBAS(ISYMBE),1)
C
               IF (LSQRAB) THEN
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),
     &                       NBAS(ISYMBE),ONE,DSRHF(KOFF1),NTOTA,
     &                       XLAMDH(KOFF2),NTOTB,ZERO,WORK(KOFF3),NTOTA)
               ELSE
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),
     &                       NBAS(ISYMBE),ONE,WORK(KOFF1),NTOTA,
     &                       XLAMDH(KOFF2),NTOTB,ZERO,WORK(KOFF3),NTOTA)
               END IF
C
               KOFF1 = IGLMRH(ISYMAL,ISYMI) + 1
               KOFF2 = KINMA1
               KOFF3 = IMAIJK(ISYMIK,ISYMJ) + NMATIJ(ISYMIK)*(J - 1)
     *               + IMATIJ(ISYMI,ISYMK) + 1
C
               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTI = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMK),NBAS(ISYMAL),
     *                    ONE,XLAMDP(KOFF1),NTOTA,WORK(KOFF2),NTOTA,
     *                    FAC,X3OINT(KOFF3),NTOTI)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      IF ((.NOT. MP2).AND.(LWRDSK)) THEN
C
C--------------------------------
C        Write integrals on disc.
C--------------------------------
C
         D    = IDEL - IBAS(ISYMD)
C
         NTOT = NMAIJK(ISYIKJ)
         IOFF = N3ODMX*(ITRAN - 1) +        !skip previous ITRAN distrib.
     &          I3ODEL(ISYIKJ,ISYMD) + NTOT*(D - 1) + 1
C
         IF (NTOT .GT. 0) THEN
            CALL PUTWA2(LUO3,FILO3,X3OINT,IOFF,NTOT)
         ENDIF
      ENDIF
C
      RETURN
      END
