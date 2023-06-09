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
      SUBROUTINE CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,ISIGN,INUM,
     &                         WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: find second-order operator from the labels of the
*              corresponding first-order perturbations and try to
*              read the operator from the AONEINT file
*              if the labels is identified and the operator is found
*              on file, it is added to the IROPER2 list
*              
*
*     Input:   LABELA,LABELB  --  labels of first-order perturbations
*
*     Output:  LABSOP    --  labels for second-order operator
*              ISYSOP    --  point group symmetry
*              INUM      --  index on IROPER2 list, -1 if not found
*              ISIGN     --  flag for sign, since for some operators
*                            the sign conventions differs from the
*                            one for the first-order operators
*              
*     Christof Haettig, March 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccroper.h"
#include "ccropr2.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*8 LABELA, LABELB, LABSOP
      LOGICAL FOUND, LOPNSAVE
      INTEGER LWORK, ISYSOP, INUM, IMATRIX, IERR, KPROPAO
      INTEGER KEND1, LWRK1, ISIGN
      INTEGER ISYMA, ISYMB, INUMA, INUMB, INUMAB

      REAL*8   WORK(LWORK)

* external functions:
      INTEGER IROPER2
      INTEGER IROPER

*---------------------------------------------------------------------*
* check if already known (if so, return without any further action):
*---------------------------------------------------------------------*
      LOPNSAVE = LOPR2OPN
      LOPR2OPN = .FALSE.
      LQUIET   = .TRUE.

      INUM = IROPER2(LABELA,LABELB,LABSOP,ISIGN,ISYSOP)

      LQUIET   = .FALSE.
      LOPR2OPN = LOPNSAVE

      IF (INUM.GT.0) RETURN

*---------------------------------------------------------------------*
* for some second-order operators the sign convention differs from the
* one for the first-order operators and we need to turn the sign
* (see below), but default is not to turn the sign.
*---------------------------------------------------------------------*
      ISIGN = +1

*---------------------------------------------------------------------*
* Hessian (geometric second derivatives): not yet available
*---------------------------------------------------------------------*
      IF    (LABELA(1:5).EQ.'1DHAM'.AND.LABELB(1:5).EQ.'1DHAM') THEN

        FOUND = .FALSE.

*---------------------------------------------------------------------*
* Dipole gradient:
*---------------------------------------------------------------------*
      ELSEIF(LABELA(2:8).EQ.'DIPLEN '.AND.LABELB(1:5).EQ.'1DHAM') THEN

        WRITE(LABSOP,'(A3,A4,A1)') LABELB(6:8), 'DPG ', LABELA(1:1)
        FOUND = .TRUE.
        ISIGN = -1

      ELSEIF(LABELB(2:8).EQ.'DIPLEN '.AND.LABELA(1:5).EQ.'1DHAM') THEN

        WRITE(LABSOP,'(A3,A4,A1)') LABELA(6:8), 'DPG ', LABELB(1:1)
        FOUND = .TRUE.
        ISIGN = -1

*---------------------------------------------------------------------*
* Second moment of charge:
*---------------------------------------------------------------------*
      ELSEIF(LABELA(3:8).EQ.'SECMOM'.AND.LABELB(1:5).EQ.'1DHAM') THEN

        WRITE(LABSOP,'(A3,A3,A2)') LABELB(6:8), 'QDG', LABELA(1:2)
        FOUND = .TRUE.
        ISIGN = -1

      ELSEIF(LABELB(3:8).EQ.'SECMOM'.AND.LABELA(1:5).EQ.'1DHAM') THEN

        WRITE(LABSOP,'(A3,A3,A2)') LABELA(6:8), 'QDG', LABELB(1:2)
        FOUND = .TRUE.
        ISIGN = -1

*---------------------------------------------------------------------*
* Third moment of charge:
*---------------------------------------------------------------------*
      ELSEIF(LABELA(5:8).EQ.'3MOM'.AND.LABELB(1:5).EQ.'1DHAM') THEN

        WRITE(LABSOP,'(A2,A3,A3)') LABELB(7:8), 'OCG', LABELA(1:3)
        FOUND = .TRUE.
        ISIGN = -1

      ELSEIF(LABELB(5:8).EQ.'3MOM'.AND.LABELA(1:5).EQ.'1DHAM') THEN

        WRITE(LABSOP,'(A2,A3,A3)') LABELA(7:8), 'OCG', LABELB(1:3)
        FOUND = .TRUE.
        ISIGN = -1

*---------------------------------------------------------------------*
* Quadrupole gradient: not yet available
*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
* Octupole gradient: not yet available
*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
* Magnetizabilities (magnetic second derivatives): not yet available
*---------------------------------------------------------------------*
      ELSEIF(LABELA(1:5).EQ.'dh/dB'.AND.LABELB(1:5).EQ.'dh/dB') THEN

        FOUND = .FALSE.

*---------------------------------------------------------------------*
* geometric derivatives of magnetic properties: 
*---------------------------------------------------------------------*
      ELSEIF(LABELA(1:5).EQ.'1DHAM'.AND.LABELB(1:5).EQ.'dh/dB') THEN

        FOUND = .FALSE.

      ELSEIF(LABELA(1:5).EQ.'dh/dB'.AND.LABELB(1:5).EQ.'1DHAM') THEN

        FOUND = .FALSE.

*---------------------------------------------------------------------*
* mixed electric/magnetic properties with London orbitals: 
*---------------------------------------------------------------------*
      ELSEIF(LABELA(2:7).EQ.'DIPLEN'.AND.LABELB(1:5).EQ.'dh/dB') THEN

        WRITE(LABSOP,'(A1,A5,A2)') LABELA(1:1), '-CM1 ', LABELB(6:7)
        FOUND = .TRUE.
        ISIGN = +1

      ELSEIF(LABELA(1:5).EQ.'dh/dB'.AND.LABELB(2:7).EQ.'DIPLEN') THEN

        WRITE(LABSOP,'(A1,A5,A2)') LABELB(1:1), '-CM1 ', LABELA(6:7)
        FOUND = .TRUE.
        ISIGN = +1

*---------------------------------------------------------------------*
* mixed electric/magnetic properties with usual (not London) orbitals: 
*---------------------------------------------------------------------*
      ELSEIF(LABELA(2:7).EQ.'DIPLEN'.AND.LABELB(2:7).EQ.'ANGMOM') THEN

        WRITE(LABSOP,'(A6,A1,A1)') '-> DxL',LABELA(1:1),LABELB(1:1)
        FOUND = .TRUE.
        ISIGN = 0

      ELSEIF(LABELA(2:7).EQ.'ANGMOM'.AND.LABELB(2:7).EQ.'DIPLEN') THEN

        WRITE(LABSOP,'(A6,A1,A1)') '-> DxL',LABELB(1:1),LABELA(1:1)
        FOUND = .TRUE.
        ISIGN = 0

*---------------------------------------------------------------------*
* two dipole operators --> second-order operator is zero 
*---------------------------------------------------------------------*
      ELSEIF(LABELA(2:7).EQ.'DIPLEN'.AND.LABELB(2:7).EQ.'DIPLEN') THEN

        WRITE(LABSOP,'(A6,A1,A1)') '-> DxD',LABELA(1:1),LABELB(1:1)
        FOUND = .TRUE.
        ISIGN = 0

*---------------------------------------------------------------------*
* nuclear shieldings:
*---------------------------------------------------------------------*
      ELSEIF(LABELA(1:5).EQ.'dh/dB'.AND.LABELB(1:4).EQ.'PSO ') THEN

        WRITE(LABSOP,'(A3,A4,A1)') LABELB(5:7),' NST',LABELA(6:6)
        FOUND = .TRUE.
        ISIGN = +1

      ELSEIF(LABELB(1:5).EQ.'dh/dB'.AND.LABELA(1:4).EQ.'PSO ') THEN

        WRITE(LABSOP,'(A3,A4,A1)') LABELA(5:7),' NST',LABELB(6:6)
        FOUND = .TRUE.
        ISIGN = +1

*---------------------------------------------------------------------*
* default: no second-order operator in the Hamiltonian
*---------------------------------------------------------------------*
      ELSE
        FOUND = .FALSE.
      END IF

*---------------------------------------------------------------------*
* check if the LABSOP is available on the AONEINT file:
*---------------------------------------------------------------------*
      IF (FOUND .AND. LABSOP(1:2).NE.'->') THEN

         KPROPAO = 1
         KEND1   = KPROPAO + N2BASX
         LWRK1   = LWORK   - KEND1

         CALL CCPRPAO(LABSOP,.TRUE.,WORK(KPROPAO),ISYSOP,IMATRIX,IERR,
     &                WORK(KEND1),LWRK1)

         IF (IERR.GT.0) FOUND = .FALSE.

         IF (IERR.LT.0) THEN
            INUMA   = IROPER(LABELA,ISYMA)
            INUMB   = IROPER(LABELB,ISYMB)
            ISYSOP  = MULD2H(ISYMA,ISYMB)
            IMATRIX = ISYMAT(INUMA) * ISYMAT(INUMB) 
            IF (LOCDBG) THEN
              WRITE(LUPRI,*) 'LABELA,ISYMA,INUMA:',LABELA,ISYMA,INUMA
              WRITE(LUPRI,*) 'LABELB,ISYMB,INUMB:',LABELB,ISYMB,INUMB
              WRITE(LUPRI,*) 'IMAT:',ISYMAT(INUMA),ISYMAT(INUMB)
            END IF
         END IF

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) '"',LABSOP,'" integrals found on file :',
     &           FOUND
         END IF

      END IF

      IF (FOUND .AND. LABSOP(1:2).EQ.'->') THEN
         INUMA   = IROPER(LABELA,ISYMA)
         INUMB   = IROPER(LABELB,ISYMB)
         ISYSOP  = MULD2H(ISYMA,ISYMB)
         IMATRIX = 0
      END IF

*---------------------------------------------------------------------*
* add operator to the IROPER2 list:
*---------------------------------------------------------------------*
      IF (FOUND) THEN
         LOPNSAVE = LOPR2OPN
         LOPR2OPN = .TRUE.

         INUM = IROPER2(LABELA,LABELB,LABSOP,ISIGN,ISYSOP)

         ! save symmetry of integral matrix on common blocks
         ISYMAT2(INUM) = IMATRIX

         ! get index of operator on IROPER list and
         ! set symmetry of integral matrix correct for this list:
         INUMAB = IROPER(LABSOP,ISYSOP)
         ISYMAT(INUMAB) = IMATRIX

         LOPR2OPN = LOPNSAVE
      ELSE
         INUM = -1
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_FIND_SO_OP> LABELA :',LABELA
        WRITE (LUPRI,*) 'CC_FIND_SO_OP> LABELB :',LABELB
        WRITE (LUPRI,*) 'CC_FIND_SO_OP> LABSOP :',LABSOP
        WRITE (LUPRI,*) 'CC_FIND_SO_OP> ISYSOP :',ISYSOP
        WRITE (LUPRI,*) 'CC_FIND_SO_OP> IMATRIX:',IMATRIX
      END IF

      RETURN
      END
*=====================================================================*
