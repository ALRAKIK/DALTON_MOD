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
*=====================================================================*
      SUBROUTINE CC_1IDXTRAN(AMAT,ISYMA,BMAT,ISYMB,CMAT,ISYMC)
*---------------------------------------------------------------------*
*
*     Purpose: evaluate the 1-index transformation
*
*           C = A B + B A^T
*
*     Christof Haettig 7-2-1999
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include "ccorb.h"
#include "ccsdsym.h"
#include "priunit.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYMA, ISYMB, ISYMC

      DOUBLE PRECISION AMAT(*), BMAT(*), CMAT(*), ONE
      PARAMETER( ONE = 1.0D0 )

      INTEGER ISYMP,ISYMQ,ISYMR,NBASP,NBASQ,NBASR,KOFF1,KOFF2,KOFF3

*---------------------------------------------------------------------*
*     check symmetries and initialize output matrix:
*---------------------------------------------------------------------*
      IF (ISYMC .NE. MULD2H(ISYMA,ISYMB)) THEN
         CALL QUIT('Symmetry mismatch in CC_1IDXTRAN.')
      END IF

      CALL DZERO(CMAT,N2BST(ISYMC))

*---------------------------------------------------------------------*
*     Calculate A x B  and add to output matrix:
*---------------------------------------------------------------------*
      DO ISYMP = 1, NSYM

         ISYMQ = MULD2H(ISYMP,ISYMA)
         ISYMR = MULD2H(ISYMQ,ISYMB)

         KOFF1 = IAODIS(ISYMP,ISYMQ) + 1
         KOFF2 = IAODIS(ISYMQ,ISYMR) + 1
         KOFF3 = IAODIS(ISYMP,ISYMR) + 1

         NBASP = MAX(1,NBAS(ISYMP))
         NBASQ = MAX(1,NBAS(ISYMQ))

         CALL DGEMM('N','N',NBAS(ISYMP),NBAS(ISYMR),NBAS(ISYMQ),
     *              ONE,AMAT(KOFF1),NBASP,BMAT(KOFF2),NBASQ,
     *              ONE,CMAT(KOFF3),NBASP)
      END DO

 
*---------------------------------------------------------------------*
*     Calculate B x A^T  and add to output matrix:
*---------------------------------------------------------------------*
      DO ISYMP = 1, NSYM

         ISYMQ = MULD2H(ISYMP,ISYMB)
         ISYMR = MULD2H(ISYMQ,ISYMA)

         KOFF1 = IAODIS(ISYMP,ISYMQ) + 1
         KOFF2 = IAODIS(ISYMR,ISYMQ) + 1
         KOFF3 = IAODIS(ISYMP,ISYMR) + 1

         NBASP = MAX(1,NBAS(ISYMP))
         NBASR = MAX(1,NBAS(ISYMR))

         CALL DGEMM('N','T',NBAS(ISYMP),NBAS(ISYMR),NBAS(ISYMQ),
     *              ONE,BMAT(KOFF1),NBASP,AMAT(KOFF2),NBASR,
     *              ONE,CMAT(KOFF3),NBASP)
      END DO

 
*---------------------------------------------------------------------*
*     print to output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_1IDXTRAN> result of one-index '//
     &        'transformation:'
         CALL CC_PRONELAO(CMAT,ISYMC)
      END IF

      RETURN
      END
*======================================================================*
