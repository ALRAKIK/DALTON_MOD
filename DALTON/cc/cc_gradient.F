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
C /* deck cc_gradient */
C=====================================================================*
       SUBROUTINE CC_GRADIENT(WORK,LWORK)
C---------------------------------------------------------------------*
C
C    Purpose: calculation of the gradient
C
C             implemented wavefunctions models:
C                 MP2, CCS, CCD, CC2, CCSD, CCSD(T)
C                 RCCD (RPA, singlet only), DRCCD (DRPA) and SOSEX
C
C     Written by Asger Halkier and Christof Haettig summer 1996.
C     CCSD(T) added, Sonia Coriani spring 2002
C     RCCD (RPA, singlet only), DRCCD (DRPA) and SOSEX by Sonia, 2011
C
C=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "nuclei.h"
#include "energy.h"
#include "cch2d.h"
#include "second.h"

* local parameters:

      REAL*8 ZERO
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      INTEGER LWORK

      REAL*8 WORK(LWORK)
      REAL*8 TIM0, GRDNRM

      LOGICAL SYM1ONLY, LPTALSO
      INTEGER MAXDIF, IOPT

* external functions
      REAL*8 DDOT

! trkoor.h : NCOOR
#include "trkoor.h"
      REAL*8 ERGMOL
      REAL*8, allocatable :: GRDMOL(:), HESMOL(:,:)

*---------------------------------------------------------------------*
* print header for gradient section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*--------  OUTPUT FROM COUPLED CLUSTER ',
     &                                  'GRADIENT SECTION  ---------*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*

* test for wavefunction model: 
!
      IF (.NOT. (CCS. OR. MP2. OR. CC2. OR.
     &           CCD .OR. CCSD. OR. CCPT .OR.
     &           RCCD .OR. DRCCD) ) THEN
         CALL QUIT('CC_GRADIENT called for a model '//
     &             'not (yet) implemented.')
      END IF

* print available work space:
      IF (IPRINT .GT. 10) THEN
        WRITE(LUPRI,*) 'Available work space in CC_GRADIENT section:',
     &        LWORK
      END IF
  
* initialize timing:
      TIM0  = SECOND()

* initialize ABACUS common blocks DORPS and CBINUC: 
* set differentiation level to 1 and print level to IPRINT
* and variables on DORPS common block as needed gradient
C     MAXDIF = 1
C     CALL CCDORPSINI(MAXDIF,IPRINT,'GRADIENT')

      SYM1ONLY = .TRUE.
      CALL CC_SETDORPS('1DHAM   ',SYM1ONLY,IPRINT)

* initialize /CCH2D/ common block, controlling the calculation
* of derivative 2-el integrals for the CC program:
      TKTIME = .TRUE.
      RETUR  = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      NOCONT = .FALSE.
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0

*---------------------------------------------------------------------*
* call ABACUS to calculate nuclear repulsion contribution to gradient:
* result is returned on GRADNN on ABACUS /ENERGY/ common block...
*---------------------------------------------------------------------*
      CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),WORK(2*MXCOOR*MXCOOR+1))

*---------------------------------------------------------------------*
* call CC_GRAD0 to calculate:
*        GRADKE -- contrib. from derivative of one-electron hamilt.
*        GRADEE -- contrib. from derivative of two-electron integr.
*        GRADFS -- contrib. from overlap derivative integr.
* arrays are allocated on ABACUS /ENERGY/ common block, GRADNA not used
*---------------------------------------------------------------------*
C
C==================================================================
C     Asger Halkier January 99. Change to two-step algorithm. First
C     call calculates everything else than GRADEE, which is instead
C     calculated in the second call to the new routine CC_GRAD2E.
C==================================================================
C     IOPT = 1
C     CALL CC_GRAD(GRADKE,GRADFS,WORK,LWORK,IOPT)
C     IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'Home from CC_GRAD.'
C     CALL CC_GRAD2E(GRADEE,WORK,LWORK)
C     IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'Home from CC_GRAD2E.'
C     CALL FLSHFO(LUPRI)
C==================================================================
C     SC & CH: Change back to one-step algorithm. Keep CC_GRAD in
C              two steps since it is used for DPT
C==================================================================
C
!SONIA: CCSD(T) is a "fake" CCSD
!
      LPTALSO = .FALSE.
      IF (CCPT) THEN
         LPTALSO = .TRUE.
         CCPT    = .FALSE.
         CCSD    = .TRUE.
         CALL CCPT_GRAD0(GRADKE,GRADFS,GRADEE,WORK,LWORK,LPTALSO)
         CCPT = .TRUE.
         CCSD = .FALSE.
      ELSE IF ((RCCD).OR.(DRCCD)) THEN
         write(lupri,*)'Warning: you requested RPA/dRPA geo gradient'
         CALL RPA_GRAD(GRADKE,GRADFS,GRADEE,WORK,LWORK)
      ELSE
         LPTALSO = .FALSE.
         CALL CC_GRAD0(GRADKE,GRADFS,GRADEE,WORK,LWORK)
      END IF

*---------------------------------------------------------------------*
* zero nuclear attraction contribution since it is included in GRADKE
*---------------------------------------------------------------------*
      CALL DZERO(GRADNA,MXCOOR)

*---------------------------------------------------------------------*
* print gradient contributions (GRADNN is already printed by NUCREP)
*---------------------------------------------------------------------*
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('One-electron integral gradient',-1)
         CALL PRIGRD(GRADKE,WORK,WORK(MXCOOR*MXCOOR+1))
         CALL HEADER('Two-electron integral gradient',-1)
         CALL PRIGRD(GRADEE,WORK,WORK(MXCOOR*MXCOOR+1))
         CALL HEADER('Reorthonormalization gradient',-1)
         CALL PRIGRD(GRADFS,WORK,WORK(MXCOOR*MXCOOR+1))
      END IF

*---------------------------------------------------------------------*
* calculated total gradient using ABACUS routines:
* (uses ABACUS common block /TAYMOL/)
*---------------------------------------------------------------------*
      NCOOR = 3*NUCDEP ! define NCOOR in trkoor.h, used in ZERGRD, ADDGRD and ABAREAD_TAYMOL
      CALL ZERGRD
      CALL ADDGRD(GRADNN)
      CALL ADDGRD(GRADEE)
      CALL ADDGRD(GRADFS)
      CALL ADDGRD(GRADKE)

      allocate ( GRDMOL(NCOOR), HESMOL(NCOOR,NCOOR) )
      CALL HEADER('Molecular gradient',-1)
      CALL ABAREAD_TAYMOL(ERGMOL,GRDMOL,HESMOL,NCOOR)
      CALL PRIGRD(GRDMOL,WORK,WORK(MXCOOR*MXCOOR+1))

      GRDNRM = DDOT(NCOOR,GRDMOL,1,GRDMOL,1)

      WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *   'Molecular gradient norm:', GRDNRM

      deallocate ( GRDMOL, HESMOL )

*---------------------------------------------------------------------*
* print timing & return:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,A,F12.2," seconds.")') ' Total time used ',
     &  'in gradient section:', SECOND() - TIM0

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_GRADIENT                          *
*=====================================================================*
c /* deck ccdorpsini */
*=====================================================================*
       SUBROUTINE CCDORPSINI(MAXDIFCC,IPRINTCC,KEY)
*---------------------------------------------------------------------*
*
*  Purpose: initialize ABACUS common blocks /DOPRS/, /CBINUC/:
*
*---------------------------------------------------------------------*
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbinuc.h"
#include "dorps.h"
#include "abainf.h"

       CHARACTER KEY*(*)
       INTEGER MAXDIFCC, IPRINTCC

*---------------------------------------------------------------------*
* copy input variables to /CBINUC/ common block:
*---------------------------------------------------------------------*
       MAXDIF = MAXDIFCC
       IPRINT = IPRINTCC

*---------------------------------------------------------------------*
* set variables on /DORPS/ common block:
*---------------------------------------------------------------------*
       IF ( KEY(1:8) .EQ. 'GRADIENT' ) THEN

C        ------------------------------------------------------------
C        for gradient we need only derivatives for totally symmetric
C        symmetry coordinates. note that in eri irepps go from 0 to 7
C        ------------------------------------------------------------
         MAXREP   = 7

         DOSYM(1) = .TRUE.
         DO IREP = 2, 8
           DOSYM(IREP) = .TRUE.
         END DO

         DOREPS(0) = .TRUE.
         DO IREP = 1, 7
           DOREPS(IREP) = .TRUE.
         END DO

         WRITE (LUPRI,*) 'DOREPS set.'
         CALL FLSHFO(LUPRI)

         DO I = 1, 3*MXCENT
           DOPERT(I,1) = .TRUE.
           DOPERT(I,2) = .FALSE.
         END DO

         WRITE (LUPRI,*) 'DOPERT set.'
         CALL FLSHFO(LUPRI)

C        ICOOR = 0
         DO IATOM = 1, NUCIND
            DO IXYZ = 1, 3
               DOCOOR(IXYZ,IATOM) = .TRUE.
C              ICOOR = ICOOR + 1
C              DO IREP = 0, 7
C                IF (DOREPS(IREP).AND.(IPTCNT(ICOOR,IREP,1).GT.0)) THEN
C                  DOCOOR(IXYZ,IATOM) = .TRUE.
C                END IF
C              END DO
            END DO
         END DO

         WRITE (LUPRI,*) 'DOCOOR set.'
         CALL FLSHFO(LUPRI)

       ELSE
         CALL QUIT('Unknown KEY in CCNUCINI.')
       END IF

       RETURN
       END
*---------------------------------------------------------------------*
C  /* Deck cc_geoflag */
      SUBROUTINE CC_GEOFLAG
C
C     Written by Asger Halkier & Christof Haettig 19/8 - 1998
C
C     Purpose: To set logical flag for integral file in case of
C              geometry optimization.
C              set RSPIM to enforce calculation of response 
C              intermediates together with the energy
C
#include "implicit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      KEEPAOTWO = MAX(KEEPAOTWO,1)
      RSPIM     = .TRUE.
      KEEPAOIN  = .TRUE.
C
      RETURN
      END
*---------------------------------------------------------------------*
