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
C  /* Deck cc_aodens3 */
      SUBROUTINE CC_AODENS3(XLAMDP1,ISYMP1,XLAMDH2,ISYMH2,
     &                      XLAMDP3,ISYMP3,XLAMDH4,ISYMH4,
     &                      DENS,ISYDEN,ICORE,IOPT,WORK,LWORK)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Calculate (special) AO-density matrix used in contructing
C     the (special) AO Fock matrix.
C
C     IOPT = 1 : XLAMDP1 x XLAMDH2 
C     IOPT = 2 : XLAMDP1 x XLAMDH2 + XLAMDP3 x XLAMDH4 
C
C     Sonia Coriani  10-Feb-1999, based on AODENS2 
C    
C     Careful: No special treatment of ICORE. 
C     ORDER IS IMPORTANT: always P1*H2 (+ P3*H4)
C     Debug 12.8.99 OK
C     Reinsert in AODENS2 at some point
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XLAMDP1(*), XLAMDH2(*), XLAMDP3(*), XLAMDH4(*)
      DIMENSION DENS(*), WORK(LWORK)
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "priunit.h"
#include "dummy.h"
C
*
* Symmetry test
*      
      ISYTOT1 = MULD2H(ISYMP1,ISYMH2) 
      IF (ISYTOT1.NE.ISYDEN) 
     *    CALL QUIT('Symmetry mismatch 1 in AODENS3' )
      IF (IOPT.EQ.2) THEN
        ISYTOT2 = MULD2H(ISYMP3,ISYMH4)
        IF (ISYTOT1.NE.ISYTOT2) 
     *      CALL QUIT('Symmetry mismatch 1 in AODENS3' )
      END IF
*
      IF (IOPT.GE.1) THEN
         CALL DZERO(DENS,N2BST(ISYDEN))
      END IF                                 
*
      DO ISYMK = 1,NSYM
C
        IF (IOPT.GE.1) THEN
C
            ISYMA = MULD2H(ISYMP1,ISYMK)      ! XLAMDP1
            ISYMB = MULD2H(ISYMH2,ISYMK)      ! XLAMDH2
C
            KOFF1 = 1 + IGLMRH(ISYMA,ISYMK)   ! XLAMDP1
            KOFF2 = 1 + IGLMRH(ISYMB,ISYMK)   ! XLAMDH2
            KOFF3 = 1 + IAODIS(ISYMA,ISYMB)
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),ONE,
     *                 XLAMDP1(KOFF1),NBASA,XLAMDH2(KOFF2),NBASB,ONE,
     *                 DENS(KOFF3),NBASA)
C
         END IF
C
         IF (IOPT.EQ.2) THEN                     !2nd contribution
C
            ISYMA = MULD2H(ISYMP3,ISYMK)         ! XLAMDP3
            ISYMB = MULD2H(ISYMH4,ISYMK)         ! XLAMDH4
C
            KOFF1 = 1 + IGLMRH(ISYMA,ISYMK)      ! XLAMDP3
            KOFF2 = 1 + IGLMRH(ISYMB,ISYMK)      ! XLAMDH4
            KOFF3 = 1 + IAODIS(ISYMA,ISYMB)
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),ONE,
     *                 XLAMDP3(KOFF1),NBASA,XLAMDH4(KOFF2),NBASB,ONE,
     *                 DENS(KOFF3),NBASA)
C
         END IF !result has been added to previous (beta=1)
C
      END DO
C
C
C-----------------------------
C     Include frozen orbitals.
C-----------------------------
C
      IF ( (FROIMP.OR.FROEXP) .AND. (ICORE .EQ. 1) ) THEN
C
         IF (IOPT.NE.0) THEN
           WRITE (LUPRI,*) 
     *              'CC_AODENS3: ICORE=1 not yet available for IOPT>0.'
           CALL QUIT(
     *              'CC_AODENS3: ICORE=1 not yet available for IOPT>0.')
         END IF
C
         IF (LWORK .LT. NLAMDS) THEN
            CALL QUIT('Insufficient space in CCSD_AODENS')
         ENDIF
C
C-------------------------------------------------
C        Read MO-coefficients from interface file.
C-------------------------------------------------
C
         LUSIFC = -1
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
C
         READ (LUSIFC)
         READ (LUSIFC) (WORK(I), I=1,NLAMDS)
C
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C-------------------------------------------------------
C        Add contribution from frozen occupied orbitals.
C-------------------------------------------------------
C
         KOFF1 = 0
         KOFF2 = 0
         DO ISYMK = 1,NSYM
C
CCH the following has been changed because ISYMP0 is not defined...?
CCH         ISYMA = MULD2H(ISYMP0,ISYMK)
CCH         ISYMB = MULD2H(ISYMP0,ISYMK)
CCH
CCH Sonia is this correct?
            ISYMA = ISYMK
            ISYMB = ISYMK
CCH
C
            DO II = 1,NRHFFR(ISYMK)
C
               K = KFRRHF(II,ISYMK)
C
               DO B = 1,NBAS(ISYMB)
                  DO A = 1,NBAS(ISYMA)
C
                     NAK = KOFF1 + NBAS(ISYMA)*(K - 1) + A
                     NBK = KOFF1 + NBAS(ISYMB)*(K - 1) + B
                     NAB = KOFF2 + NBAS(ISYMA)*(B - 1) + A
C
                     DENS(NAB) = DENS(NAB) + WORK(NAK)*WORK(NBK)
C
                  END DO
               END DO
C
            END DO
C
            KOFF1 = KOFF1 + NBAS(ISYMK)*NORBS(ISYMK)
            KOFF2 = KOFF2 + NBAS(ISYMA)*NBAS(ISYMB)
C
         END DO
C
      ENDIF
C
      END
