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
C  /* Deck anainp */
      SUBROUTINE ANAINP(WORD)
C
C  5-Jul-1985 Hans Jorgen Aa. Jensen
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 4)
      PARAMETER (MAXANG = 20)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "abainf.h"
      LOGICAL SKIP
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
      DATA TABLE /'.SKIP  ', '.XXXXXX', '.ANGLES', '.DIHEDR'/
      DATA MANG/0/, MDIHED/0/
C
      CALL ANAINI
C
      NEWDEF = (WORD .EQ. '*GEOANA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ANAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ANAINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) NANG
                  MANG = MIN(MAXANG,NANG)
                  DO 310 I = 1,MANG
                     READ(LUCMD,*) (IANG(J,I),J=1,3)
  310             CONTINUE
                  MANG = NANG - MANG
                  DO 320 I = 1,MANG
                     READ(LUCMD,'()')
  320             CONTINUE
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) NDIHED
                  MDIHED = MIN(MAXANG,NDIHED)
                  DO 410 I = 1,MDIHED
                     READ(LUCMD,*) (IDIHED(J,I),J=1,4)
  410             CONTINUE
                  MDIHED = NDIHED - MDIHED
                  DO 420 I = 1,MDIHED
                     READ(LUCMD,'()')
  420             CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ANAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ANAINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for GEOANA:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' GEOANA skipped in this run.'
         ELSE
            IF (NANG .GT. 0) THEN
               WRITE (LUPRI,'(/A/)')
     *            ' Following angles will be calculated:'
               DO 1310 I = 1,NANG
                  WRITE (LUPRI,'(I10,A,4I5)') I,' : ',(IANG(J,I),J=1,3)
 1310          CONTINUE
               IF (MANG .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3,A)') ' The last',MANG,
     *               ' angles specified go beyond current maximum',
     *               ' and will not be printed.'
               END IF
            END IF
            IF (NDIHED .GT. 0) THEN
               WRITE (LUPRI,'(/A/)')
     *            ' Following dihedral angles will be calculated:'
               DO 1410 I = 1,NDIHED
                  WRITE (LUPRI,'(I10,A,4I5)')I,' : ',(IDIHED(J,I),J=1,4)
 1410          CONTINUE
               IF (MDIHED .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3,A)') ' The last',MDIHED,
     *               ' dihedral angles specified go beyond current',
     *               ' maximum and will not be printed.'
               END IF
            END IF
         END IF
         WRITE (LUPRI,'(/)')
      END IF
      RETURN
      END
C  /* Deck anaini */
      SUBROUTINE ANAINI
C
C     Initialize /CBIANA/
C
#include "implicit.h"
      PARAMETER (MAXANG = 20)
      LOGICAL SKIP
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
      NANG   = 0
      NDIHED = 0
      SKIP   = .FALSE.
      RETURN
      END
C  /* Deck geoana */
      SUBROUTINE GEOANA(COORD,PRINT,DIF,NBONDS,LUPUNCH,WORK,LWORK)
C
C 30-Jun-1985 Hans Jorgen Aa. Jensen
C Modified for symmetry 25-Sep-1989 tuh
C Modified for differential geomtries 18-Oct-1989 tuh
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL   PRINT, DIF
      DIMENSION COORD(*),WORK(LWORK)
#include "nuclei.h"
C
      CALL QENTER('GEOANA')
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KVEC ,3*NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KBOND  ,NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KCHRG  ,NUCDEP,       WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPAIR,2*NUCDEP*NUCDEP,WORK,KFREE,LFREE)
C
      CALL GEOANA_1(COORD,PRINT,DIF,NBONDS,LUPUNCH,WORK(KVEC),
     &            WORK(KBOND),WORK(KCHRG),WORK(KPAIR))
C
      CALL MEMREL('GEOANA',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('GEOANA')
      RETURN
      END
C  /* Deck GEOANA_1 */
      SUBROUTINE GEOANA_1(COORD,PRINT,DIF,NBONDS,LUPUNCH,VEC,
     &                  BONDED,ICHARG,IPAIRS)
C
C Modified for more selective printing of bonded atoms,
C     Jan-1995 Hanne Heiberg
C     LUPUN .gt. 0: punching atom bonds for Gamess graphic output, K.Ruud-95
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
#include "facang.h"
      PARAMETER (MAXANG = 20)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D1P5 = 1.5D00)
      LOGICAL   SKIP, DIF, PRINT
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "nuclei.h"
#include "symmet.h"
#include "qm3.h"
C
      DIMENSION COORD(3,*), DIST(MXCENT*(MXCENT+1)/2),
     &          ANGLE(MAXANG), DIHED(MAXANG), ICHARG(NUCDEP)
      DIMENSION VEC(3,NUCDEP,NUCDEP), IPAIRS(2,NUCDEP*NUCDEP)
      LOGICAL   BONDED(NUCDEP,NUCDEP)
      CHARACTER*6 NUCNAM(4)
      SAVE DIST, ANGLE, DIHED
C     statement function:
      ARCCOS(ARG) = FACANG*ACOS(ARG)

C
#ifndef PRG_DIRAC
      IF (REDCNT) RETURN
!     ... REDCNT is a QM3 variable
#endif

      NBONDS = 0
      IF (NUCDEP .EQ. 1) RETURN
C
C     set up bond vectors in Angstrom in VEC
C     NUCIND is number of QM atoms
C     NCTOT .ge. NUCIND is number of QM atoms + number of MM atoms
C
      N_QMATOMS = 0
      IATOMA = 0
      DO 100 ICENTA = 1, NCTOT
         DO 110 IA = 0, MAXOPR
            IF (IAND(IA,ISTBNU(ICENTA)) .EQ. 0) THEN
               IATOMA = IATOMA + 1
               CXA = PT(IAND(ISYMAX(1,1),IA))*COORD(1,ICENTA)
               CYA = PT(IAND(ISYMAX(2,1),IA))*COORD(2,ICENTA)
               CZA = PT(IAND(ISYMAX(3,1),IA))*COORD(3,ICENTA)
C
               ICHARG(IATOMA) = IZATOM(ICENTA)
C
               IATOMB = 0
               DO 200 ICENTB = 1, NCTOT
                  DO 210 IB = 0, MAXOPR
                     IF (IAND(IB,ISTBNU(ICENTB)) .EQ. 0) THEN
                        IATOMB = IATOMB + 1
                        IF (IATOMB .GT. IATOMA) GO TO 110
C                       ... next IATOMA, only IATOMB .le. IATOMA needed
                          CXB=PT(IAND(ISYMAX(1,1),IB))*COORD(1,ICENTB)
                          CYB=PT(IAND(ISYMAX(2,1),IB))*COORD(2,ICENTB)
                          CZB=PT(IAND(ISYMAX(3,1),IB))*COORD(3,ICENTB)
                          VEC(1,IATOMB,IATOMA) = XTANG*(CXA - CXB)
                          VEC(2,IATOMB,IATOMA) = XTANG*(CYA - CYB)
                          VEC(3,IATOMB,IATOMA) = XTANG*(CZA - CZB)
                          VEC(1,IATOMA,IATOMB) = -VEC(1,IATOMB,IATOMA)
                          VEC(2,IATOMA,IATOMB) = -VEC(2,IATOMB,IATOMA)
                          VEC(3,IATOMA,IATOMB) = -VEC(3,IATOMB,IATOMA)
                     END IF
  210             CONTINUE
  200          CONTINUE
            END IF
  110    CONTINUE
         IF (ICENTA .EQ. NUCIND) N_QMATOMS = IATOMA
  100 CONTINUE
      IF (IATOMA .NE. NUCDEP) THEN
        WRITE(LUPRI,*) 'GEOANA error, IATOMA .ne. NUCDEP:',IATOMA,NUCDEP
        CALL QUIT('NCTOT and NUCDEP inconsistent')
      END IF
C
C
C     Set up distance matrix in Angstrom
C
      DIST_MAX   = D0
      ABSDIST_MAX= D0
      IDIST_MAX  = 0
      JDIST_MAX  = 0
      ADIST_MAX  = D0
      IADIST_MAX = 0
      JADIST_MAX = 0

      N_SHORT_HX = 0
      N_SHORT_YX = 0
      ADISTHX_MIN= 1.D200
      ADISTYX_MIN= 1.D200

      IJ = 0
      DO 400 I = 1,NUCDEP
         DO 300 J = 1,I
            IJ = IJ + 1
            DISTAN = VEC(1,J,I)*VEC(1,J,I) + VEC(2,J,I)*VEC(2,J,I)
     *             + VEC(3,J,I)*VEC(3,J,I)
            DISTAN = SQRT(DISTAN)
            IF (DIF) THEN
               DISTAN = DISTAN - DIST(IJ)
            ELSE IF (I .LE. N_QMATOMS .AND. I.NE.J) THEN
               IF (ICHARG(I) .EQ. 0 .OR. ICHARG(J) .EQ. 0) THEN
                  ! do nothing - we do not want to include floating orbitals in minimum atom distance
               ELSE IF (ICHARG(I) .NE. 1 .AND. ICHARG(J) .NE. 1) THEN
                  ADISTYX_MIN = MIN(ADISTYX_MIN,DISTAN)
                  IF (DISTAN .LE. 1.0D0) N_SHORT_YX = N_SHORT_YX + 1 ! R(Y-X) .lt. 1.0 Angstrom is usually an error
               ELSE
                  ADISTHX_MIN = MIN(ADISTHX_MIN,DISTAN)
                  IF (DISTAN .LE. 0.7D0) N_SHORT_HX = N_SHORT_HX + 1 ! R(H-X) .lt. 0.7 Angstrom is usually an error
               END IF
            END IF
            DIST(IJ)   = DISTAN
            IF (ABS(DISTAN) .GT. ABSDIST_MAX) THEN
               DIST_MAX  = DISTAN
               ABSDIST_MAX = ABS(DIST_MAX)
               IDIST_MAX = I
               JDIST_MAX = J
               IF (I .LE. N_QMATOMS) THEN
C                 ... these are QM atoms
                  ADIST_MAX  = ABS(DISTAN)
                  IADIST_MAX = I
                  JADIST_MAX = J
               END IF
            END IF
  300    CONTINUE
  400 CONTINUE
C
      IF (PRINT) THEN
         JPRIDIS = MAX(1,IPRUSR)
         IF (DIST_MAX .NE. D0) THEN ! exclude atoms (DIST_MAX .eq. 0.D0)

           IF (NUCDEP .LE. 16*JPRIDIS) THEN
C          hjaaj Oct 2003: this output is only useful for small molecules ...
            IF (DIF) THEN
              CALL HEADER
     *        ('Differential interatomic separations (in Angstrom):',2)
            ELSE
              CALL HEADER('Interatomic separations (in Angstrom):',2)
            END IF
            CALL PRIDIS(NAMDEP,DIST,NUCDEP)
           END IF
C
           IF (DIF) THEN
            WRITE(LUPRI,'(/A,F12.6,A/A,I5,A,I5,5A)')
     &      '  Max differential change in interatomic separation is',
     &      ADIST_MAX,' Angstrom',
     &      '  between atoms',IADIST_MAX,' and',JADIST_MAX,
     &      ', "',NAMDEP(IADIST_MAX),'" and "',NAMDEP(JADIST_MAX),'".'
           ELSE
            WRITE(LUPRI,'(/A,2(F10.4,A)/A,I5,A,I5,5A)')
     &      '  Max    interatomic separation is',
     &       ADIST_MAX,' Angstrom (',ADIST_MAX/XTANG,' Bohr)',
     &      '  between atoms',IADIST_MAX,' and',JADIST_MAX,
     &      ', "',NAMDEP(IADIST_MAX),'" and "',NAMDEP(JADIST_MAX),'".'
            IF (ADISTHX_MIN .LT. 1.D10) WRITE(LUPRI,'(/A,2(F10.4,A))')
     &      '  Min HX interatomic separation is',
     &       ADISTHX_MIN,' Angstrom (',ADISTHX_MIN/XTANG,' Bohr)'
            IF (ADISTYX_MIN .LT. 1.D10) WRITE(LUPRI,'(/A,2(F10.4,A))')
     &      '  Min YX interatomic separation is',
     &       ADISTYX_MIN,' Angstrom (',ADISTYX_MIN/XTANG,' Bohr)'

            IF (N_SHORT_YX .GT. 0 .OR. N_SHORT_HX .GT. 0) THEN
              NWARN = NWARN + 1
              WRITE(LUPRI,'(/A,2I5/A/A)')
     &          '@ WARNING: Number of short HX and YX bond lengths:',
     &            N_SHORT_HX,N_SHORT_YX,
     &          '@ WARNING: If not intentional, maybe your coordinates'
     &            //' were in Angstrom,',
     &          '@ WARNING: but "Angstrom" was not specified'
     &            //' in .mol file'
            END IF

           END IF
           IF (N_QMATOMS .NE. NUCDEP) THEN
            IF (DIF) THEN
             WRITE(LUPRI,'(/A,F12.6,A/A,I5,A,I5,5A)')
     &         '  Max differential change in QM+MM interatomic '//
     &           'separation is',DIST_MAX,' Angstrom',
     &         '  between the QM+MM centers',IDIST_MAX,' and',JDIST_MAX,
     &         ', "',NAMDEP(IDIST_MAX),'" and "',NAMDEP(JDIST_MAX),'".'
            ELSE
             WRITE(LUPRI,'(/A,2(F10.4,A)/A,I5,A,I5,5A)')
     &         '  Max QM+MM interatomic separation is',
     &         DIST_MAX,' Angstrom (',DIST_MAX/XTANG,' Bohr)',
     &         '  between the QM+MM centers',IDIST_MAX,' and',JDIST_MAX,
     &         ', "',NAMDEP(IDIST_MAX),'" and "',NAMDEP(JDIST_MAX),'".'
            END IF
           END IF
         END IF
      END IF
C
      IF (PRINT .AND. .NOT. DIF) THEN

        IJ = 0
        DO 10 J= 1,N_QMATOMS
C     ... only QM atoms, not MM atoms
          RADJ = RADIUS(ICHARG(J))
          IF (RADJ .LT. 0.0D0 .AND. ICHARG(J) .GT. 0) THEN
C             do not print if cavity center or floating orbital /hjaaj
              WRITE(LUPRI,*)
     &        'INFO: RADIUS FOR ATOM WITH ATOMIC NUMBER ',
     &        ICHARG(J),' IS UNAVAILABLE, USING 2.0 AA'
            RADJ = 2.0D0
          END IF
          DO 20 I= 1, J-1
            IJ = IJ + 1
            RADI = RADIUS(ICHARG(I))
            IF (RADI .LT. 0.0D0 .AND. ICHARG(I) .GT. 0) THEN
              RADI = 2.0D0
            END IF
            IF (RADI .LT. 0 .OR. RADJ .LT. 0) THEN
C             not bonded if cavity center or floating orbital /hjaaj
              BONDED(I,J) = .FALSE.
              BONDED(J,I) = .FALSE.
            ELSE IF (DIST(IJ) .LT. (1.2D0 * (RADI + RADJ))) THEN
              NBONDS = NBONDS + 1
              IPAIRS(1,NBONDS) = I
              IPAIRS(2,NBONDS) = J
              BONDED(I,J) = .TRUE.
              BONDED(J,I) = .TRUE.
            ELSE
              BONDED(I,J) = .FALSE.
              BONDED(J,I) = .FALSE.
            END IF
  20     CONTINUE
          IJ = IJ + 1
         BONDED(J,J) = .FALSE.
  10  CONTINUE
C
      IF (.NOT. QM3) THEN
          CALL HEADER('Bond distances (Angstrom):',1)
          WRITE (LUPRI,'(18X,A/18X,A)')
     &          'atom 1     atom 2       distance',
     &          '------     ------       --------'
          IJ = 0
          DO I = 1, N_QMATOMS
            DO J = 1, I-1
              IJ = IJ + 1
              IF (BONDED(I,J)) THEN
                NUCNAM(1) = NAMDEP(I)
                NUCNAM(2) = NAMDEP(J)
                WRITE(LUPRI,'(A,2X,A6,5X,A6,F15.6)')
     &               '  bond distance:',
     &               NUCNAM(1), NUCNAM(2), DIST(IJ)
              END IF
            END DO
            IJ = IJ + 1
          END DO
C
          IF (N_QMATOMS .GT. 2 .AND. NANG .LE. 0) THEN
            CALL HEADER('Bond angles (degrees):',1)
            WRITE (LUPRI,'(18X,A/18X,A)')
     $            'atom 1     atom 2     atom 3         angle',
     $            '------     ------     ------         -----'
C
            IJK = 0
            DO 40, I= 1,N_QMATOMS
              DO 50, J= 1, N_QMATOMS - 1
                DO 60, K= J + 1, N_QMATOMS
                  IF (BONDED(I,J) .AND. BONDED(I,K)) THEN
                    IJK = IJK + 1
                    NUCNAM(1) = NAMDEP(J)
                    NUCNAM(2) = NAMDEP(I)
                    NUCNAM(3) = NAMDEP(K)
                    ANG = FACANG*VECANG(VEC(1,I,J),VEC(1,I,K))
                    WRITE(LUPRI,'(A,5X,A6,5X,A6,5X,A6,F14.3)')
     &                 '  bond angle:',NUCNAM(1),NUCNAM(2),NUCNAM(3),ANG
                  END IF
 60             CONTINUE
 50           CONTINUE
 40         CONTINUE
            IF (IJK .EQ. 0) WRITE(LUPRI,'(5X,A)') 'No angles found'
          END IF
        END IF ! not QM3
      END IF
C
C     Punch bonding information in Gamess output format on unit LUPUNCH
C
      IF (LUPUNCH .GT. 0) THEN
         IF(NBONDS.LE.6) THEN
            WRITE(LUPUNCH,8010) (IPAIRS(1,I),IPAIRS(2,I),I=1,NBONDS)
         ELSE
            WRITE(LUPUNCH,8020) (IPAIRS(1,I),IPAIRS(2,I),I=1,6)
            WRITE(LUPUNCH,8030) (IPAIRS(1,I),IPAIRS(2,I),I=7,NBONDS)
         END IF
      END IF
C         
      IF (NANG .GT. 0) THEN
         IF (PRINT) THEN
            CALL HEADER('Angles according to input list:',2)
            WRITE (LUPRI,'(A/A)')
     *       '    atom 1     atom 2     atom 3         angle (degrees)',
     *       '    ------     ------     ------         ---------------'
         END IF
         DO 1000 I = 1,NANG
            I1 = IANG(1,I)
            I2 = IANG(2,I)
            I3 = IANG(3,I)
            IMX = MAX(I1,I2,I3)
            IF (IMX .GT. N_QMATOMS) THEN
               IF (PRINT) WRITE (LUPRI,'(/A/)')
     &            ' *GEOANA input error for .ANGLES: non-existent atom'
               GO TO 1000
            END IF
            NUCNAM(1) = NAMDEP(I1)
            NUCNAM(2) = NAMDEP(I2)
            NUCNAM(3) = NAMDEP(I3)
            IF (I1 .NE. I2 .AND. I2 .NE. I3) THEN
               ANG = FACANG*VECANG(VEC(1,I2,I1),VEC(1,I2,I3))
               IF (.NOT.DIF) THEN
                  ANGLE(I) = ANG
               ELSE
                  ANGLE(I) = ANG - ANGLE(I)
               END IF
               IF (PRINT) WRITE (LUPRI,'(4X,A6,5X,A6,5X,A6,F20.3)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),ANGLE(I)
            ELSE
               IF (PRINT) WRITE (LUPRI,'(4X,A6,5X,A6,5X,A6,10X,A)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),'undefined'
            END IF
 1000    CONTINUE
      END IF
C
      IF (NDIHED .GT. 0) THEN
          IF (PRINT) WRITE (LUPRI,'(//A/A)')
     *       '    atom 1     atom 2     atom 3     atom 4'//
     *       '    dihedral angle (degrees)',
     *       '    ------     ------     ------     ------'//
     *       '    ------------------------'
         DO 2000 I = 1,NDIHED
            I1 = IDIHED(1,I)
            I2 = IDIHED(2,I)
            I3 = IDIHED(3,I)
            I4 = IDIHED(4,I)
            IMX = MAX(I1,I2,I3,I4)
            IF (IMX .GT. N_QMATOMS) THEN
               IF (PRINT) WRITE (LUPRI,'(/A/)')
     &            ' *GEOANA input error for .DIHEDR: non-existent atom'
               GO TO 2000
            END IF
            NUCNAM(1) = NAMDEP(I1)
            NUCNAM(2) = NAMDEP(I2)
            NUCNAM(3) = NAMDEP(I3)
            NUCNAM(4) = NAMDEP(I4)
            X1 = VEC(2,I2,I1)*VEC(3,I2,I3) - VEC(2,I2,I3)*VEC(3,I2,I1)
            X2 = VEC(3,I2,I1)*VEC(1,I2,I3) - VEC(3,I2,I3)*VEC(1,I2,I1)
            X3 = VEC(1,I2,I1)*VEC(2,I2,I3) - VEC(1,I2,I3)*VEC(2,I2,I1)
            Y1 = VEC(2,I3,I2)*VEC(3,I3,I4) - VEC(2,I3,I4)*VEC(3,I3,I2)
            Y2 = VEC(3,I3,I2)*VEC(1,I3,I4) - VEC(3,I3,I4)*VEC(1,I3,I2)
            Y3 = VEC(1,I3,I2)*VEC(2,I3,I4) - VEC(1,I3,I4)*VEC(2,I3,I2)
            Z1 = X2*Y3 - X3*Y2
            Z2 = X3*Y1 - X1*Y3
            Z3 = X1*Y2 - X2*Y1
            SENSE = Z1*VEC(1,I2,I3) + Z2*VEC(2,I2,I3) + Z3*VEC(3,I2,I3)
            SENSE = SIGN(D1,SENSE)
            ANG = X1*Y1 + X2*Y2 + X3*Y3
            DDD = (X1*X1 + X2*X2 + X3*X3) * (Y1*Y1 + Y2*Y2 + Y3*Y3)
            IF (DDD .GT. 1.D-10) THEN
               ANG = ANG / SQRT(DDD)
               IF (ABS(ANG) .GT. D1) ANG = SIGN(D1,ANG)
               ANG = SENSE*ARCCOS(ANG)
               IF (.NOT.DIF) THEN
                  DIHED(I) = ANG
               ELSE
                  DIHED(I) = ANG - DIHED(I)
               END IF
               IF (PRINT) WRITE(LUPRI,'(4X,A6,5X,A6,5X,A6,5X,A6,F20.3)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),NUCNAM(4),DIHED(I)
            ELSE
               IF (PRINT) WRITE(LUPRI,'(4X,A6,5X,A6,5X,A6,5X,A6,10X,A)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),NUCNAM(4),'undefined'
            END IF
 2000    CONTINUE
      END IF
C
      IF (PRINT) WRITE (LUPRI,'(/)')
      RETURN
 8010 FORMAT('BONDATOMS ',6(I4,I4,2X))
 8020 FORMAT('BONDATOMS ',6(I4,I4,2X),' >')
 8030 FORMAT(7(I4,I4,2X),:,' >')
      END
C  /* Deck pridis */
      SUBROUTINE PRIDIS (NAMDEP,DISMAT,NROW)
C
C 30-Jun-1985 Hans Jorgen Aa. Jensen
C (based on OUTPAK by Nelson H.F. Beebe)
C
C Print bond distance matrix (or other matrix over atoms)
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (KCOL=6)
      CHARACTER*6 NAMDEP(*)
      DIMENSION DISMAT(*)
      INTEGER BEGIN
C
      LAST = MIN(NROW,KCOL)
      BEGIN = 1
 1050 NCOL = 1
      WRITE (LUPRI,1000) (NAMDEP(I),I = BEGIN,LAST)
      WRITE (LUPRI,1000) ('------' ,I = BEGIN,LAST)
      DO 40 K = BEGIN,NROW
         KTOTAL = (K*(K-1))/2 + BEGIN - 1
         WRITE (LUPRI,2000) NAMDEP(K),
     *      (DISMAT(KTOTAL+J),J = 1,NCOL)
         IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40 CONTINUE
      WRITE (LUPRI,'()')
      LAST = MIN(LAST+KCOL,NROW)
      BEGIN = BEGIN+NCOL
      IF (BEGIN.LE.NROW) GO TO 1050
      RETURN
 1000 FORMAT (8X,6(4X,A6,2X))
 2000 FORMAT (1X,A6,':',6F12.6)
      END
C  /* Deck radius */
      FUNCTION RADIUS(NCHARGE)
#include "implicit.h"
#include "priunit.h"
C
C     Based on covalent radii and metallic radii in Angstrom.
C     Returns -1 where data is inavailable
C     Oct 2006 hjaaj: changed Hydrogen from 30 to 40 pm,
C              such that H2 is printed as bonded ;-) .
C
      DIMENSION RAD(100)
      DATA (RAD(I), I = 1, 100)/
     &        40.,  155.,  160.,  110.,
     & 90.,   80.,   70.,   68.,   65.,
     &154.,  190.,  160.,  140.,  110.,
     &110.,  105.,  105.,  190.,  238.,
     &200.,  165.,  145.,  135.,  130.,
     &125.,  125.,  125.,  125.,  125.,
     &140.,  140.,  130.,  120.,  120.,
     &120.,  200.,  255.,  215.,  180.,
     &160.,  145.,  140.,  135.,  130.,
     &130.,  135.,  140.,  155.,  160.,
     &160.,  140.,  140.,  140.,  220.,
     &270.,  220.,  185.,  180.,  180.,
     &180.,  180.,  180.,  200.,  180.,
     &175.,  175.,  175.,  175.,  170.,
     &170.,  170.,  155.,  145.,  140.,
     &135.,  135.,  135.,  135.,  145.,
     &155.,  170.,  175.,  170.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100./
C
      IF (NCHARGE .LE. 0) THEN
Chj      =  0: solvent cavity center or floating orbital
Chj      = -Z: multiple basis for r12 methods for nuclear charge Z
Chj      = -1234567890: point charges
C
         RADIUS = -1.0D0
      ELSE IF (NCHARGE .LT. 1 .OR. NCHARGE .GT. 100) THEN
         WRITE (LUPRI,*) 
     &          'ERROR, RADIUS called with CHARGE =',NCHARGE
         CALL QUIT('RADIUS called with unvalid CHARGE')
      ELSE
         RADIUS = 0.01D0 * RAD(NCHARGE)
      END IF
      RETURN
      END
C  /* Deck vdwrad */
      FUNCTION VDWRAD(NCHARGE)
#include "implicit.h"
#include "priunit.h"
C     Based on van der Waals radii in Angstrom.
C     Returns -1 where data is inavailable
      DIMENSION RAD(100)
      DATA (RAD(I), I = 1, 100)/
     &       110.,  220.,  122.,   63.,
     &155.,  155.,  140.,  135.,  130.,
     &154.,  190.,  160.,  140.,  110.,
     &202.,  220.,  150.,  150.,  220.,
     &188.,  181.,  175.,  277.,  239.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100.,   -100.,   -100.,   -100.,   -100.,
     & -100./
C
      IF (NCHARGE .LE. 0) THEN
Chj      =  0: solvent cavity center or floating orbital
Chj      = -Z: multiple basis for r12 methods for nuclear charge Z
Chj      = -1234567890: point charges
C
         VDWRAD = -1.0D0
      ELSE IF (NCHARGE .GT. 100) THEN
         WRITE (LUPRI,*) 'ERROR, VDWRAD called with CHARGE =',NCHARGE
         CALL QUIT('VDWRAD called with illegal CHARGE')
      ELSE
         VDWRAD = 0.01D0 * RAD(NCHARGE)
         IF (VDWRAD .LT. 0.0D0) THEN
C           if no table value, use covalent radius plus 0.65 AA
C           (as for B-F above) /Dec. 2006, hjaaj
            WDWRAD = RADIUS(NCHARGE) 
            IF (VDWRAD .GT. 0.0D0) VDWRAD = VDWRAD + 0.65D0
         END IF
      END IF
      RETURN
      END
