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
C  /* Deck herdi1 */
      SUBROUTINE HERDI1(WORK,LWORK,IPRFCK)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
      LOGICAL   ABA, OLDDX
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include "ccinftap.h"
#include "inforb.h"
#include "cbitwo.h"
#include "nuclei.h"
#include "dorps.h"
#include "inftap.h"
#include "symmet.h"
#include "blocks.h"
#include "eribuf.h"
C
C
C     Inquire if this is an ABACUS run or not:
C
      CALL ABARUN(ABA)
C
C----------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C----------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
C
         IPRALL =  0
C
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
C
C
         IPRINT = IPRFCK
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
C     define NIBUF and NBITS in eribuf.h
C
      CALL ERIBUF_INI
C
      IF (LUINTR .LE. 0) CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUINTR
C
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KJORBS = KIORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KKORBS = KJORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KLAST  = KKORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('HERDI1','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
C
      RETURN
      END
C  /* Deck herdi2 */
      SUBROUTINE HERDI2(WORK,LWORK,INDEXA,ISHELA,NAINTS,IPRFCK)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"

CCOLD LOGICAL DOINDX, RELCAL, RETUR, TKTIME, NOCONT, NOPV, NODV
      LOGICAL DOINDX, RELCAL
      INTEGER NUMDIS
C
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include "inforb.h"
#include "cbitwo.h"
#include "dorps.h"
#include "inftap.h"
#include "nuclei.h"
#include "symmet.h"
#include "blocks.h"
#include "eritap.h"
#include "eribuf.h"
#include "r12int.h"
C
      CALL QENTER('HERDI2')
      IF (LUAORC(0) .LE. 0) THEN
         LBFINP = LBUF
C
         IF (NIBUF .ne. 1 .and. NIBUF .ne. 2) THEN
            CALL QUIT('HERDI2 error: NIBUF .ne 1 .and. NIBUF .ne. 2')
         END IF
#if defined (SYS_NEC)
         LRECL =   LBFINP + NIBUF*LBFINP/2 + 1   ! integer*8 units
#else
         LRECL = 2*LBFINP + NIBUF*LBFINP   + 1   ! integer*4 units
#endif
         CALL GPOPEN(LUAORC(0),'AO2DIS00','UNKNOWN','DIRECT',
     &        'UNFORMATTED',LRECL,OLDDX)
         IF (U21INT) THEN
            LU21INT = -1
            CALL GPOPEN(LU21INT,'AOTDIS00','UNKNOWN','DIRECT',
     &                  'UNFORMATTED',LRECL,OLDDX)
         END IF
      END IF
C
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KJORBS = KIORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KKORBS = KJORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KLAST  = KKORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('HERDI2','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      ITYPE  = 7
      MAXDIF = 0
      JPRINT = 0
      NBUFX(0) = 0
      IPRINT = 0
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      NDMAT  = 0
      IREPDM = 0
      IFCTYP = 0
      INDXAB = 0
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      I2TYP  = 0
      ICEDIF = 0
      IFTHRS = 0
      DOINDX = .TRUE.
      NEWDIS = .TRUE.
      RETUR  = .FALSE.
      TKTIME = .FALSE.
      NOCONT = .FALSE.
      NOPV   = .FALSE.
      NODV   = .FALSE.
      RELCAL = .FALSE.
C
      REWIND(LUINTR)
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
C
      CALL AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,WORK(KIORBS),JPRINT)
C
      NUMDIS = NAINTS
C
      IF (IPRFCK.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (LUINTA .LE. 0) CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      CALL TWOINT(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,IREPDM,
     &            IFCTYP,GMAT,INDXAB,NUMDIS,MAXDIS,
     &            ITYPE,MAXDIF,JATOM,NODV,
     &            NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &            IPRNTD,RETUR,ISHELA,I2TYP,
     &            WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,
     &            .false.)
      IF (IPRFCK.GT.0) CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
      IF (U21INT) CALL GPCLOSE(LU21INT,'KEEP')
C
      CALL QEXIT('HERDI2')
      RETURN
      END
C
C  /* Deck nerdi2 */
      SUBROUTINE NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHELA,ISHELB,NAINTS,
     &                  NBINTS,IPRFCK,SETUP,IOPT)
C
C     Option input :   IOPT = 1  :  Calculate ordenary integral distributions 
C                                   with two fixed indices.
C                           = 2  :  Calculate exchange type integrals entering
C                                   the diagonal elements.
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
#include "choles.h"

      LOGICAL DOINDX, RELCAL,SETUP
      INTEGER NUMDIS
C
      DIMENSION INDEXA(*)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,?
C
#include "inforb.h"
#include "cbitwo.h"
#include "dorps.h"
#include "inftap.h"
#include "nuclei.h"
#include "symmet.h"
#include "blocks.h"
#include "eritap.h"
#include "eribuf.h"
C
      IF (LUAORC(0) .LE. 0) THEN
         LBFINP = LBUF
C
         CALL ERIBUF_INI  ! set NIBUF, NBITS, IBIT1, IBIT2
#if defined (SYS_NEC)
         LRECL =   LBFINP + NIBUF*LBFINP/2 + 1   ! integer*8 units
#else
         LRECL = 2*LBFINP + NIBUF*LBFINP   + 1   ! integer*4 units
#endif
         CALL GPOPEN(LUAORC(0),'AO2DIS00','UNKNOWN','DIRECT',
     &        'UNFORMATTED',LRECL,OLDDX)
      END IF
C
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 - 1)/IRAT + 1
      KJORBS = KIORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KKORBS = KJORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      KLAST  = KKORBS + (MXSHEL*MXAOVC   - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('NERDI2','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      ITYPE  = 7
      MAXDIF = 0
      JPRINT = 0
      NBUFX(0) = 0
      IPRINT = 0
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      NDMAT  = 0
      IREPDM = 0
      IFCTYP = 0
      INDXAB = 0
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      I2TYP  = 0
      ICEDIF = 0
      IFTHRS = 0
      DOINDX = .TRUE.
      NEWDIS = .TRUE.
      RETUR  = .FALSE.
      TKTIME = .FALSE.
      NOCONT = .FALSE.
      NOPV   = .FALSE.
      NODV   = .FALSE.
      RELCAL = .FALSE.
C
      ISHLB = ISHELB
C
      IF     (IOPT .EQ. 1) THEN
         I2TYP = 4
      ELSEIF (IOPT .EQ. 2) THEN
         I2TYP = 5
c        LSAVE = LBUF
c        LBUF  = 50 000
      ELSE
          CALL QUIT('Wrong I2TYP in NERDI2')
      ENDIF
C
      REWIND(LUINTR)
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
C
      CALL AINDEX(ISHELB,NBINTS,INDEXB,DOINDX,WORK(KIORBS),JPRINT)
      CALL AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,WORK(KIORBS),JPRINT)
C
      IF (SETUP) RETURN
C
Casm  NUMDIS = NAINTS
C
      IF (IPRFCK.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (LUINTA .LE. 0) CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      CALL TWOINT(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,IREPDM,
     &            IFCTYP,GMAT,INDXAB,NUMDIS,MAXDIS,
     &            ITYPE,MAXDIF,JATOM,NODV,
     &            NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &            IPRNTD,RETUR,ISHELA,I2TYP,
     &            WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,
     &            .false.)
      IF (IPRFCK.GT.0) CALL TIMER('TWOINT',TIMSTR,TIMEND)
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
C
      RETURN
      END

