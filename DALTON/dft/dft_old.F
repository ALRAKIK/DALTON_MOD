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
#ifdef INCLUDE_DFT_OLD
      SUBROUTINE DFT_COMP_KSM(DMAT,KSM,WRK,LWORK,IPRFCK)
#include "implicit.h"
#include "inforb.h"
      DIMENSION DMAT(NBAST,NBAST), KSM(NBAST,NBAST), WRK(LWORK)
      KFREE = 1+NBAST*NBAST
      IF(KFREE.GT.LWORK) STOP "no mem in DFT_COMP_KSM"
      LF = LWORK-KFREE+1
      CALL DZERO(WRK,NBAST*NBAST)
      CALL DFTEXC(DMAT,1,WRK,1,.FALSE.,.TRUE.,.FALSE.,.FALSE.,
     &     .FALSE.,.FALSE.,DUMMY,.FALSE.,WRK(KFREE),LF,IPRFCK)
      CALL DAXPY(NBAST*NBAST,1D0,WRK,1,KSM,1)
      END
C /* Deck dftexc */
      SUBROUTINE DFTEXC(DMAT,NDMAT,EXCMAT,NXCMAT,MAKDEN,DOERG,DOGRD,
     &                  DOLND,DOATR,TRPLET,HES,DOHES,WORK,LWORK,KPRINT)
C
C     T. Helgaker sep 99 / feb 01
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
C 
      PARAMETER (FACTHR = 1.0D-8) 
      LOGICAL MAKDEN, DOERG, DOGRD, DOLND, DOATR, DOHES, TRPLET, DOGGA,
     &        DODRC, DOVWN, DOBCK, DOLYP
      DIMENSION DMAT(NBAST,NBAST,NDMAT), EXCMAT(NBAST,NBAST,NXCMAT),
     &          HES(*), WORK(LWORK)
C
#include "inforb.h"
#include "nuclei.h"
#include "dftcom.h"
#include "dftinf.h"
#include "dfterg.h"
#include "functionals.h"
C
      IPRINT = 5
      NBUF   = 1000000
C
C     Square DMAT if necessary
C 
      IF (MAKDEN) CALL DFTDNS(DMAT,WORK,LWORK,IPRINT)
C
C     Calculate grid
C
      IF (.NOT.DFTGRID_DONE_OLD) THEN
         CALL MAKE_DFTGRID(WORK,LWORK,1,.FALSE.)
         DFTGRID_DONE_OLD = .TRUE.
      END IF
C
      NC0  = 0
      IF (DOHES) NC0  = NORBT
C
      FACDRC = WDFTX
      FACVWN = WDFTC
      FACBCK = WDFTB
      FACLYP = WDFTL
C
      DODRC = DABS(FACDRC) .GT. FACTHR
      DOVWN = DABS(FACVWN) .GT. FACTHR 
      DOBCK = DABS(FACBCK) .GT. FACTHR 
      DOLYP = DABS(FACLYP) .GT. FACTHR 
      DOGGA = DFT_ISGGA()
C
C     Number of AOs and their addresses
C
      CALL SETUPSOS(DOGRD,DOGGA,DOLND,DFTPOT)
C
C     Allocations 
C
      KX    = 1
      KY    = KX    + NBUF
      KZ    = KY    + NBUF
      KW    = KZ    + NBUF
      KGSO  = KW    + NBUF
      KCNT  = KGSO  + NTYPSO*NBAST 
      KDST  = KCNT  + NBAST
      KC0   = KDST  + NATOMS 
      KC1   = KC0   + NC0
      KC2   = KC1   + 3*NC0
      KDGA  = KC2   + NC0
      KLST  = KDGA  + NBAST
      LWRK  = LWORK - KLST + 1
      IF (KLST.GT.LWORK) CALL STOPIT('DFTEXC','DFTDRV',KLST,LWORK)
      CALL DFTDRV(DMAT(1,1,1),DMAT(1,1,2),EXCMAT,NXCMAT,DOERG,DOGRD,
     &            DOLND,DOATR,WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &            WORK(KGSO),WORK(KCNT),WORK(KDST),WORK(KC0),WORK(KC1),
     &            WORK(KC2),WORK(KDGA),HES,DOHES,TRPLET,FACDRC,FACVWN,
     &            FACBCK,FACLYP,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &            WORK(KLST),LWRK,IPRINT)
      END
C /* Deck dftdrv */
      SUBROUTINE DFTDRV(DMAT,DTRMAT,EXCMAT,NXCMAT,DOERG,DOGRD,DOLND,
     &                  DOATR,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,NCNT,DST,
     &                  C0,C1,C2,DMAGAO,HES,DOHES,TRPLET,FACDRC,FACVWN,
     &                  FACBCK,FACLYP,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &                  WORK,LWORK,IPRINT)
C
C     T. Helgaker sep 99 / feb 01
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, DP5 = 0.5D0, DP3 = 1D0/3D0)
c      PARAMETER (DFTHR0 = 1.0D-11, DFTHRL = 1.0D-14, DFTHRI = 1.0D-15)
C
#include "inforb.h"
#include "nuclei.h"
#include "dfterg.h"
#include "energy.h"
#include "dftcom.h"
#include "dftinf.h"
#include "orgcom.h"
#include "functionals.h"
#include "maxorb.h"
#include "infpar.h"
#include "expopt.h"
C
      INTEGER A, B
      LOGICAL DOERG, DOGRD, DODRC, DOVWN, DOBCK, DOLYP, 
     &        DOLND, DOATR, FROMVX, DOGGA, TRPLET, DOHES 
      DIMENSION CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF), 
     &          DMAT(NBAST,NBAST), DTRMAT(NBAST,NBAST),
     &          GSO(NBAST*NTYPSO), EXCMAT(NBAST,NBAST,NXCMAT), 
     &          NCNT(NBAST), DST(NATOMS), C0(NORBT), C1(NORBT,3), 
     &          C2(NORBT),DMAGAO(NBAST), HES(NOCCT,NVIRT,NOCCT,NVIRT), 
     &          RHG(3), WORK(LWORK), VX(DERVS1), COR(3)
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C
      IF (IPRINT .GT. 100) THEN
         CALL HEADER('DMAT in DFTDRV ',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IF (DOATR) THEN
            CALL HEADER('DTRMAT in DFTDRV ',-1)
            CALL OUTPUT(DTRMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
C
      EXCTRO = 0D0
      IF (DOERG) EXCTRO = DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
      IF (DOGRD) CALL DZERO(GRADFT,MXCOOR)
C
      FROMVX = DOGGA .AND. DFTPOT
C
      ELCTRN = D0
      ENERGY = D0
C
      VD  = D0
      VV  = D0
      VB1 = D0 
      VB2 = D0 
      VL1 = D0 
      VL2 = D0 
C
      LUQUAD = -1
      CALL GPOPEN(LUQUAD,'DALTON.QUAD','OLD','SEQUENTIAL',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      R2TOT = D0
      NPNTS = 0
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
         DO 300 IPNT = 1, NPOINT
            IF (IPRINT .GT. 100) THEN
               WRITE (LUPRI,'(2X,I6,4F12.6)') 
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF
C
            WGHT = WEIGHT(IPNT)
c     2002.02.25, pawsa: The grid weights are sometimes wrong on
c     IEEE-conformant architectures (weight<=0) - skip them. New grid
c     generator will probably resolve this.
            IF(WGHT.LE.0D0) GO TO 300
C
C           AOs
C           ===
C
            THRINT = DFTHRI/WGHT
            COR(1) = CORX(IPNT)
            COR(2) = CORY(IPNT)
            COR(3) = CORZ(IPNT)
            CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,
     &                  NBAST,DOLND,DOGGA,THRINT,IPRINT)
C
C           Density
C           =======
C
            CALL GETRHO(DMAT,GSO,RHO,DMAGAO,THRINT)
            if(RHO.LT.0D0) PRINT *, "RHO=", RHO, WGHT
C
            IF (RHO.GT.DFTHR0) THEN 
               RHO13 = RHO**(DP3)
C
C              Gradient of density
C              ===================
C
               IF (DOGGA) THEN
                  CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAO,1,
     &                       D0,RHG,1)
                  RHOGRD = DSQRT(RHG(1)**2 + RHG(2)**2 + RHG(3)**2)
               END IF

#ifdef DEBUG_DFT_NUMERICALLY
               CALL DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
#endif
C              Hessian of density
C              ===================
C
               IF (FROMVX) THEN
                  CALL DFTRHH(DMAT,DMAGAO,GSO(KSO0),GSO(KSO1),GSO(KSO2),
     &                        RHG,RHOLAP,RHOGHG)
               END IF
C
C              Number of electrons
C              ===================
C
               ELCTRN = ELCTRN + WGHT*RHO
C
C              Energy
C              ======
C
               IF (DOERG) THEN
                  ENERGY = ENERGY + DFTENERGY(RHO,RHOGRD)*WGHT
               END IF
C
C              Exchange-correlation potential 
C              ==============================
C
               IF (DOERG.OR.DOGRD.OR.DOLND) THEN                 
                  CALL DFTPOT0(VX, WGHT, RHO, RHOGRD)
                  IF(DOGGA) VX(FZ0) = 2D0*VX(FZ0)/RHOGRD
               END IF
C
C              Exchange-correlation contribution to Kohn-Sham matrix
C              =====================================================
C
               IF (DOERG) THEN
                  CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG,VX(FR0),
     &                        VX(FZ0), DOGGA,FROMVX,DFTHRL)
               END IF
C
C              Exchange-correlation contribution to molecular gradient
C              =======================================================
C
               IF (DOGRD) THEN
                 IF (EXPGRA) THEN
                      CALL DFTEXP(DMAT,GSO(KSO0),GSO(KSO1),GSO(KSO2),
     &                            VX(FR0),VX(FZ0),RHG(1),RHG(2),RHG(3),
     &                            CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                            DOGGA)
                   ELSE
                      CALL DFTFRC(DMAT,GSO(KSO0),GSO(KSO1),GSO(KSO2),
     &                            VX(FR0),VX(FZ0),RHG,DOGGA)
                   END IF
               END IF
C
C              London contributions to Kohn-Sham matrix 
C              ========================================
C
               IF (DOLND) THEN
                  CALL DFTMAG(EXCMAT,
     &                        CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                        GSO(KSO0),GSO(KSO1),GSO(KSOB),GSO(KSOB1),
     &                        VX(FR0),VX(FZ0),RHG,DOGGA,
     &                        FROMVX)
               END IF
C
C              Hessian transformation
C              ======================
C
               IF (DOATR) THEN
                  CALL DFTLTR(DTRMAT,EXCMAT,GSO(KSO0),GSO(KSO1),
     &                        C0,C1,C2,HES,RHO,RHOGRD,RHG,WGHT,
     &                        DOGGA,TRPLET,DOHES,DMAGAO)
               END IF
            END IF
C
  300    CONTINUE
C
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      END IF
  200 CONTINUE
C
      CALL GPCLOSE(LUQUAD,'KEEP')
C
      IF (DOERG) THEN
         EDFTY = ENERGY+EXCTRO-DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
      END IF
      IF (DOERG) THEN
         DO I = 1, NBAST
         DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,1) + EXCMAT(J,I,1))
            EXCMAT(I,J,1) = AVERAG 
            EXCMAT(J,I,1) = AVERAG 
         END DO
         END DO
      END IF 
      IF (DOATR) THEN
         IF (DOHES) THEN
            DO B = 1, NVIRT
            DO J = 1, NOCCT 
            DO A = B + 1, NVIRT
            DO I = 1, NOCCT 
               HES(I,A,J,B) = HES(J,B,I,A)
            END DO
            END DO
            END DO
            END DO
         ELSE
            DO I = 1, NBAST
            DO J = 1, I - 1 
               AVERAG = EXCMAT(I,J,1) + EXCMAT(J,I,1)
               EXCMAT(I,J,1) = AVERAG 
               EXCMAT(J,I,1) = AVERAG 
            END DO
            END DO
         END IF
      END IF 
      IF (DOLND) THEN
         DO K = 1, 3
         DO I = 1, NBAST
         DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,K) - EXCMAT(J,I,K))
            EXCMAT(I,J,K) =   AVERAG 
            EXCMAT(J,I,K) = - AVERAG 
         END DO
         END DO
         END DO
      END IF 
C
C     Print section
C
      IF (IPRINT .GT. 2) THEN 
         WRITE (LUPRI,'(/2X,A,F14.7,6X,D8.2,I14)')
     &      ' Number of electrons/abscissas:  ', 
     &        ELCTRN,ELCTRN-NINT(ELCTRN),NPNTS
         IF (DOERG) THEN
c            WRITE (LUPRI,'(2X,A,3F14.6)') 
c     &         ' Dirac and Becke exchange energy:',       
c     &           ERGDRC, ERGBCK, ERGDRC+ERGBCK
c             WRITE (LUPRI,'(2X,A,3F14.6)') 
c     &          ' VWN and LYP correlation energy :', 
c     &          ERGVWN, ERGLYP, ERGVWN + ERGLYP
         WRITE (LUPRI,'(2X,A,28X,F14.6)') 
     &          ' DFT exchange-correlation energy:', 
     &          ENERGY
         END IF
      END IF
      IF (IPRINT .GT. 100) THEN
         IF (DOERG) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DOLND) THEN
            CALL HEADER('EXCMAT in DFTDRV (x direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,1),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (y direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,2),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (z direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,3),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DOGRD) THEN
            CALL HEADER('GRADFT in DFTDRV ',-1)
            CALL OUTPUT(GRADFT,1,3,1,NATOMS,3,NATOMS,1,LUPRI)
         END IF
         IF (DOATR) THEN
            CALL HEADER('EXCMAT after DFTDRV ',-1)
            CALL OUTPUT(EXCMAT(1,1,1),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      CALL TIMER('DFTDRV',TIMSTR,TIMEND)
      RETURN 
      END
C  /* Deck dftatr */
      SUBROUTINE DFTATR(FMAT,CMO,ZYMAT,TRPLET,KSYMOP,WORK,LWORK)
C
C     T. Helgaker Sep 99 / Feb 01
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
      PARAMETER (D2 = 2.0D0)
C
      DIMENSION CMO(*),ZYMAT(NORBT,NORBT),FMAT(NORBT,NORBT),WORK(LWORK)
      LOGICAL DOPRNT, TRPLET
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
C
      DOPRNT = .FALSE.
      IPRINT = 2
C
      KDAO = 1
      KDA1 = KDAO + N2BASX
      KFAO = KDA1 + N2BASX
      KLST = KFAO + N2BASX
      LWRK = LWORK - KLST + 1
      IF (KLST .GT. LWORK) CALL STOPIT('DFTATR','LWORK',KLST,LWORK)
C
C     AO density matrix
C
      JKEEP = JWOPSY
      JWOPSY = 1
      CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDAO),DUMMY,CMO,DUMMY,
     &            WORK(KLST),LWRK)
      JWOPSY = JKEEP
      IF (DOPRNT) THEN
         CALL HEADER('AO density matrix in DFTATR ',-1)
         CALL OUTPUT(WORK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     AO transformation matrix
C
      IF (DOPRNT) THEN
         CALL HEADER('ZYMAT in DFTATR ',-1)
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL DEQ27(CMO,ZYMAT,DUMMY,WORK(KDA1),DUMMY,WORK(KLST),LWRK)
      CALL DSCAL(N2BASX,0.5D0,WORK(KDA1),1)
      IF (DOPRNT) THEN
         CALL HEADER('AO transformation matrix in DFTATR ',-1)
         CALL OUTPUT(WORK(KDA1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     AO Fock matrix contribution
C
      CALL DZERO(WORK(KFAO),N2BASX)
      CALL DFTEXC(WORK(KDAO),2,WORK(KFAO),1,
     &            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.,
     &            TRPLET,DUMMY,.FALSE.,WORK(KLST),LWRK,IPRINT)
      IF (DOPRNT) THEN
         CALL HEADER('AO Fock matrix contribution in DFTATR ',-1)
         CALL OUTPUT(WORK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     MO Fock matrix contribution 
C
      CALL DZERO(WORK(KDA1),N2ORBX)
      DO ISYM = 1, NSYM
         JSYM  = MULD2H(ISYM,KSYMOP)
         NORBI = NORB(ISYM)
         NORBJ = NORB(JSYM)
         IF (NORBI.GT.0 .AND. NORBJ.GT.0) THEN 
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &                 WORK(KFAO),NBAS,NBAST,WORK(KDA1),NORB,
     &                 NORBT,WORK(KLST),LWRK)
         END IF
      END DO
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix contribution in DFTATR ',-1)
         CALL OUTPUT(WORK(KDA1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Add to Fock matrix
C
      IF (DOPRNT) THEN
         CALL HEADER('Original MO Fock matrix in DFTATR ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL DAXPY(NORBT*NORBT,D2,WORK(KDA1),1,FMAT,1) 
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix in DFTATR ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck dftatx */
      SUBROUTINE DFTATX(FMAT,CMO,ZYMAT,TRPLET,WORK,LWORK)
C
C     T. Helgaker Sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
      PARAMETER (NOC = 7, NVT = 1)
C
      INTEGER A, B
      DIMENSION CMO(*),ZYMAT(NORBT,NORBT),FMAT(NORBT,NORBT),WORK(LWORK)
      DIMENSION HES(NOC*NVT*NOC*NVT) 
      LOGICAL DOPRNT,FIRST,TRPLET
      SAVE FIRST, HES
      DATA FIRST/.TRUE./
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
C
      DOPRNT = .FALSE.
      IPRINT = 2
C
      NDIM1 = (NOCCT*NVIRT)**2
      NDIM2 = (NOC*NVT)**2
      IF (NDIM1 .GT. NDIM2) CALL STOPIT('DFTATX','HES',NDIM1,NDIM2)
C
      KDAO = 1
      KDA1 = KDAO + N2BASX
      KFAO = KDA1 + N2BASX
      KLST = KFAO + N2BASX
      IF (KLST .GT. LWORK) CALL STOPIT('DFTATX','LWORK',KLST,LWORK)
      LWRK = LWORK - KLST + 1
C
C     AO density matrix
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDAO),DUMMY,CMO,DUMMY,
     &               WORK(KLST),LWRK)
         CALL DZERO(HES,NOCCT*NVIRT*NOCCT*NVIRT)
         CALL DCOPY(NBAST*NORBT,CMO,1,WORK(KDA1),1)
         IF (DOPRNT) THEN
            CALL HEADER('AO density matrix in DFTATX ',-1)
            CALL OUTPUT(WORK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('MO coefficient matrix in DFTATX ',-1)
            CALL OUTPUT(WORK(KDA1),1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
         END IF
         CALL DFTEXC(WORK(KDAO),2,WORK(KFAO),1,
     &               .FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.,
     &               TRPLET,HES,.TRUE.,WORK(KLST),LWRK,IPRINT)
         IF (DOPRNT) THEN
            NPAIR = NOCCT*NVIRT
            CALL HEADER('Hessian matrix in DFTATX ',-1)
            CALL OUTPUT(HES,1,NPAIR,1,NPAIR,NPAIR,NPAIR,1,LUPRI)
         END IF
      END IF
      IF (DOPRNT) THEN
         CALL HEADER('ZYMAT in DFTATX ',-1)
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IHES = 0
      DO A = 1, NVIRT
      DO I = 1, NOCCT
         SUM = D0 
         DO B = 1, NVIRT
         DO J = 1, NOCCT
            IHES = IHES + 1
            SUM = SUM + HES(IHES)*(ZYMAT(J,B+NOCCT) - ZYMAT(B+NOCCT,J))
         END DO
         END DO
         WORK(KDA1 + NORBT*(NOCCT + A - 1) + I - 1) = SUM
         WORK(KDA1 + NORBT*(I - 1) + NOCCT + A - 1) = SUM 
      END DO
      END DO
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix contribution in DFTATX ',-1)
         CALL OUTPUT(WORK(KDA1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Add to Fock matrix
C
      IF (DOPRNT) THEN
         CALL HEADER('Original MO Fock matrix in DFTATX ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL DAXPY(NORBT*NORBT,D2,WORK(KDA1),1,FMAT,1) 
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix in DFTATX ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
c /* Deck dftlnd */
      SUBROUTINE DFTLND(FX,FY,FZ,WORK,LWORK,IPRINT)
C
C     T. Helgaker oct 99
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      DIMENSION FX(N2BASX), FY(N2BASX), FZ(N2BASX), WORK(LWORK)
#include "inforb.h"
C
      KXCMAT = 1
      KDAO   = KXCMAT + 3*N2BASX
      KLAST  = KDAO   +   N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTLND',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL DZERO(WORK(KXCMAT),3*N2BASX)
      CALL DFTEXC(WORK(KDAO),1,WORK(KXCMAT),3,.TRUE.,.FALSE.,.FALSE.,
     &            .TRUE.,.FALSE.,.FALSE.,
     &            DUMMY,.FALSE.,WORK(KLAST),LWRK,IPRINT)
C
      KXC = KXCMAT - 1
      KYC = KXCMAT - 1 +   N2BASX
      KZC = KXCMAT - 1 + 2*N2BASX
      DO I = 1, N2BASX
         FX(I) = FX(I) + WORK(KXC + I)
         FY(I) = FY(I) + WORK(KYC + I)
         FZ(I) = FZ(I) + WORK(KZC + I)
      END DO
C
      RETURN
      END
C  /* Deck dftsol */
      SUBROUTINE DFTSOL(WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION WORK(*)
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"
C
      KVEC = 1
      KFCK = KVEC + NWOPPT
      KLST = KFCK + NNORBT - 1
      IF (KLST .GT. LWORK) CALL STOPIT('DFTSOL',' ',KLST,LWORK)
      CALL DFTSO1(WORK(KVEC),WORK(KFCK))
      RETURN
      END
C  /* Deck dftso1 */
      SUBROUTINE DFTSO1(VEC,FCK)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
      PARAMETER (FAC = -1.0D0/4.0D0)
      DIMENSION VEC(*), FCK(*)
#include "inftap.h"
#include "infvar.h"
#include "inflin.h"
#include "inforb.h"
      IA(I) = JWOP(1,I)
      IB(I) = JWOP(2,I)
      IAA(I) = IA(I)*(IA(I)+1)/2
      IBB(I) = IB(I)*(IB(I)+1)/2
C
      REWIND (LUSIFC)
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) 
      READ (LUSIFC) 
      READ (LUSIFC)
      READ (LUSIFC) 
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NNORBT,FCK)
      DO 100 I = 1,NWOPPT
         VEC(I) = FAC*VEC(I)/(FCK(IAA(I)) - FCK(IBB(I)))
  100 CONTINUE
      IF (.FALSE.) THEN
         WRITE(LUPRI,'(//,A)')' Orbital part of solution'
         PRFAC = 0.1D0
         CALL PRKAP(NWOPPT,VEC,PRFAC,LUPRI)
      END IF
      RETURN
      END
C
C /* Deck dftltr */
      SUBROUTINE DFTLTR(DTRMAT,EXCMAT,GAO,GAO1,C0,C1,C2,HES,RHO,
     &                  RHOGRD,RHG,WGHT,DOGGA,TRPLET,DOHES,DTGAO)
C
C     T. Helgaker sep 99/oct 00
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
C
#include "inforb.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
#include "wrkrsp.h"
#include "functionals.h"
C
      INTEGER A, B
      LOGICAL DOGGA, TRPLET, DOHES
      DIMENSION DTRMAT(NBAST,NBAST),
     &          GAO(NBAST), GAO1(NBAST,3), 
     &          EXCMAT(NBAST,NBAST), 
     &          C0(NORBT), C1(NORBT,3), C2(NORBT),
     &          HES(NOCCT,NVIRT,NOCCT,NVIRT),RHG(3),DTGAO(NBAST)
      DIMENSION B3(3), VXC(DERVS2)
C
      IF (DOHES .AND.NSYM.GT.1) THEN
         WRITE (LUPRI,'(2X,A,/A)') 
     &      ' Symmetry not implemented for explicit Hessian in DFTLTR.',
     &      ' Program aborted.'
         CALL QUIT('Symmetry not implementd for DOHES in DFTLTR')
      END IF
      IF (DOHES) THEN
         B0 = D1
      ELSE
         CALL DGEMV('N',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
         B0 = DDOT(NBAST,DTGAO,1,GAO,1)
      END IF
C
C     ***************
C     ***** LDA *****
C     ***************
C
      IF (.NOT.DOGGA) THEN
         IF (DABS(B0).GT.DFTHRL) THEN 
            CALL DFTPOT1(VXC, WGHT,RHO, 0D0, TRPLET)
            VT = VXC(FRR)*B0
C
C           Linear transformation
C
            IF (.NOT.DOHES) THEN
               IF (NSYM.EQ.1) THEN
                  DO I = 1, NBAST
                     GVI = VT*GAO(I)
                     DO J = 1, I
                        EXCMAT(J,I) = EXCMAT(J,I) + GVI*GAO(J)
                     END DO
                  END DO
               ELSE
                  DO ISYM = 1, NSYM
                     JSYM = MULD2H(ISYM,KSYMOP) 
                     IF (ISYM.GE.JSYM) THEN
                        ISTR = IBAS(ISYM) + 1
                        IEND = IBAS(ISYM) + NBAS(ISYM)
                        JSTR = IBAS(JSYM) + 1
                        JEND = IBAS(JSYM) + NBAS(JSYM)
                        DO I = ISTR, IEND
                           GVI = VT*GAO(I)
                           DO J = JSTR, MIN(I,JEND) 
                              EXCMAT(J,I) = EXCMAT(J,I) + GVI*GAO(J)
                           END DO
                        END DO
                     END IF
                  END DO
               END IF
C
C           Explicit Hessian
C
            ELSE
               CALL DGEMM('T','N',NORBT,1,NBAST,1.D0,
     &                    DTRMAT,NBAST,
     &                    GAO,NBAST,0.D0,
     &                    C0,NORBT)
               DO 300 B = NOCCT + 1, NORBT
               DO 300 J = 1, NOCCT
               DO 300 A = NOCCT + 1, NORBT
               DO 300 I = 1, NOCCT
                  IA = A - NOCCT
                  IB = B - NOCCT
                  HES(I,IA,J,IB) = HES(I,IA,J,IB) 
     &                           + VT*C0(A)*C0(I)*C0(B)*C0(J)
  300          CONTINUE
            END IF
         END IF
      ELSE
C
C        ***************
C        ***** GGA *****
C        ***************
C
C        B0, BX=B3(1), BY=B3(2), BZ=B3(3), BMAX
C
         IF (DOHES) THEN
            BMAX = D1
         ELSE
C           B3 = GAO1'*DTGAO
            CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D0,B3,1)
C           DTGAO= DTRMAT'*GAO
            CALL DGEMV('T',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
C           B3 = B3 + GAO1'*DTGAO
            CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D1,B3,1)
            BMAX = MAX(DABS(B0),DABS(B3(1)),DABS(B3(2)),DABS(B3(3)))
         END IF
C
         IF (BMAX.GT.DFTHRL) THEN
C
C           ZNV, FZ0, FRR, FRZ, FZZ
C
            CALL DFTPOT1(VXC, WGHT,RHO, RHOGRD, TRPLET)
            ZNV = D1/RHOGRD
            VXC(FZ0) = znv*VXC(FZ0) ! FZ0 appears always with this prefactor
            RX = ZNV*RHG(1) 
            RY = ZNV*RHG(2)
            RZ = ZNV*RHG(3)
C
C           Linear transformation
C
            IF (.NOT.DOHES) THEN
               BR = B3(1)*RX + B3(2)*RY + B3(3)*RZ
               FAC0 = VXC(FRR)*B0 + VXC(FRZ)*BR
               FACR = VXC(FRZ)*B0 + VXC(FZZ)*BR
               IF (NSYM.EQ.1) THEN
                  DO I = 1, NBAST
                     G0 = GAO(I)
                     GX = GAO1(I,1)
                     GY = GAO1(I,2)
                     GZ = GAO1(I,3)
                     DO J = 1, I 
                        A0 = G0*GAO(J)
                        AX = GX*GAO(J) + G0*GAO1(J,1)
                        AY = GY*GAO(J) + G0*GAO1(J,2)
                        AZ = GZ*GAO(J) + G0*GAO1(J,3)
                        AR = AX*RX + AY*RY + AZ*RZ
                        AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) - AR*BR
                        EXCMAT(J,I) = EXCMAT(J,I)+FAC0*A0+FACR*AR
     $                              +VXC(FZ0)*AB
                     END DO
                  END DO
               ELSE
                  DO ISYM = 1, NSYM
                     ISTR = IBAS(ISYM) + 1
                     IEND = IBAS(ISYM) + NBAS(ISYM)
                     JSYM = MULD2H(ISYM,KSYMOP)
                     IF (ISYM.GE.JSYM) THEN
                        JSTR = IBAS(JSYM) + 1
                        JEND = IBAS(JSYM) + NBAS(JSYM)
                        DO I = ISTR, IEND
                           G0 = GAO(I)
                           GX = GAO1(I,1)
                           GY = GAO1(I,2)
                           GZ = GAO1(I,3)
                           DO J = JSTR, MIN(I,JEND) 
                              A0 = G0*GAO(J)
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR = AX*RX + AY*RY + AZ*RZ
                              AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) -AR*BR
                              EXCMAT(J,I) = EXCMAT(J,I) + FAC0*A0 
     &                                    + FACR*AR + VXC(FZ0)*AB
                           END DO
                        END DO
                     END IF
                  END DO
               END IF
C
C           Explicit Hessian
C
            ELSE
               CALL DGEMM('T','N',NORBT,1,NBAST,1.D0,
     &                    DTRMAT,NBAST,
     &                    GAO,NBAST,0.D0,
     &                    C0,NORBT)
               CALL DGEMM('T','N',NORBT,3,NBAST,1.D0,
     &                    DTRMAT,NBAST,
     &                    GAO1,NBAST,0.D0,
     &                    C1,NORBT)
C
               DO I = 1, NORBT
                 C2(I) = RX*C1(I,1)+RY*C1(I,2)+RZ*C1(I,3)
               END DO
C
               DO B = NOCCT + 1, NORBT 
               DO J = 1, NOCCT 
                 IB  = B - NOCCT
                 GBJ = C0(B)*C0(J)
                 PBJ = C2(B)*C0(J)   + C2(J)*C0(B)
                 CBX = C0(B)*C1(J,1) + C1(B,1)*C0(J)
                 CBY = C0(B)*C1(J,2) + C1(B,2)*C0(J) 
                 CBZ = C0(B)*C1(J,3) + C1(B,3)*C0(J)
                 DO A = NOCCT + 1, B
                 DO I = 1, NOCCT
                   IA  = A - NOCCT
                   GAI = C0(A)*C0(I)
                   PAI = C2(A)*C0(I) + C2(I)*C0(A)
                   CAB = CBX*(C0(A)*C1(I,1) + C1(A,1)*C0(I))
     &                 + CBY*(C0(A)*C1(I,2) + C1(A,2)*C0(I))
     &                 + CBZ*(C0(A)*C1(I,3) + C1(A,3)*C0(I))
                   HES(I,IA,J,IB) = HES(I,IA,J,IB) 
     &                            + VXC(FRR)*GAI*GBJ
     &                            + VXC(FRZ)*(PAI*GBJ + GAI*PBJ)
     &                            + VXC(FZZ)*PAI*PBJ
     &                            + VXC(FZ0)*(CAB - PAI*PBJ)
                  END DO
                  END DO
               END DO
               END DO
            END IF
         END IF
      END IF
      RETURN
      END
C
C /* Deck dftest */
      SUBROUTINE DFTTST(WORK,LWORK,NBUF,LUQUAD,IPRINT)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C 
      DIMENSION WORK(LWORK)
C
#include "inforb.h"
#include "nuclei.h"
C
      KX    = 1
      KY    = KX   + NBUF
      KZ    = KY   + NBUF
      KW    = KZ   + NBUF
      KGAO  = KW   + NBUF
      KGA1  = KW   + NBUF
      KDEN  = KGA1 + 3*NBAST
      KCNT  = KDEN + NBAST*NBAST
      KLST  = KCNT + NBAST
      LWRK  = LWORK - KLST + 1
      CALL DFTES1(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &            WORK(KDEN),WORK(KGAO),WORK(KGA1),WORK(KCNT),
     &            WORK(KLST),LWRK,IPRINT,LUQUAD)
      RETURN
      END
C /* Deck dftest */
      SUBROUTINE DFTES1(CORX,CORY,CORZ,WEIGHT,NBUF,DMAT,GAO,GAO1,NCNT,
     &                  WORK,LWORK,IPRINT,LUQUAD)
C
C     T. Helgaker sep 99
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
C
#include "inforb.h"
#include "dftinf.h"
C
      DIMENSION CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF),
     &          DMAT(NBAST,NBAST), GAO(NBAST), GAO1(NBAST,3),
     &          NCNT(NBAST), WORK(LWORK)
C
      LOGICAL NODV, NODC, DOLND
C
      DOLND = .FALSE.
C
C     Calculate density matrix
C     
      NODC = .FALSE.
      NODV = .FALSE. 
      CALL GETDEN(DMAT,WORK,LWORK,NODC,NODV,IPRINT)
      ELCTRN = D0
C
      NDER = 0
C
C     Integrate
C
      REWIND LUQUAD
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
         DO 300 IPNT = 1, NPOINT
            CALL DFTAOS(GAO,GAO1,NCNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                  NBAST,DOLND,IPRINT)
            DO 400 I = 1, NBAST
            DO 400 J = 1, NBAST
               ELCTRN = ELCTRN + WEIGHT(IPNT)*DMAT(I,J)*GAO(I)*GAO(J)
  400       CONTINUE 
  300    CONTINUE
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      ELSE
         GO TO 200
      END IF
  200 CONTINUE
C
      WRITE (LUPRI,'(2X,A,F20.12)') ' Number of electrons ',ELCTRN
C
      RETURN 
      END
C  /* Deck dftnum */
      SUBROUTINE DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, 
     &           DP3 = D1/D3, ADD = 1.0D-06)
      LOGICAL DODRC,DOVWN,DOLYP,DOBCK
C
      ADD2 = D2*ADD
C
      R0 = RHO
      RP = RHO + ADD
      RM = RHO - ADD
      Z0 = RHOGRD
      ZP = RHOGRD + ADD 
      ZM = RHOGRD - ADD
      R013 = R0**DP3
      RP13 = RP**DP3
      RM13 = RM**DP3
C
      IF (DODRC) THEN
         CALL VDRC(G10,R013)
         CALL V1DRC(F10,F20,R0,R013)
C
         CALL EDRC(EP,RP,RP13)
         CALL VDRC(V1P,RP13)
         CALL EDRC(EM,RM,RM13)
         CALL VDRC(V1M,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' DR010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' DRC20 test',F20,H20,F20-H20
      END IF
      IF (DOVWN) THEN
         CALL VVWN(G10,R0,R013)
         CALL V1VWN(F10,F20,R0,R013)
C
         CALL EVWN(EP,RP,RP13)
         CALL VVWN(V1P,RP,RP13)
         CALL EVWN(EM,RM,RM13)
         CALL VVWN(V1M,RM,RM13)
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' VW010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' VWN20 test',F20,H20,F20-H20
      END IF
C
      IF (DOBCK) THEN
         CALL GBCK(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1BCK(F10,F01,F20,F11,F02,R0,Z0)
C
         CALL EBCK(EP,RP,RP13,Z0)
         CALL GBCK(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL EBCK(EM,RM,RM13,Z0)
         CALL GBCK(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
C
         CALL EBCK(EP,R0,R013,ZP)
         CALL GBCK(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL EBCK(EM,R0,R013,ZM)
         CALL GBCK(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' BC010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BC001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' BCK01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' BCK20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' BCK11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' BCK02 test',F02,H02,F02-H02
      END IF
      IF (DOLYP) THEN
         CALL GLYP(G10,G01,R0,R013,Z0)
         G01 = DP5*G01*Z0 
         CALL V1LYP(F10,F01,F20,F11,F02,R0,Z0)
C
         CALL ELYP(EP,RP,RP13,Z0)
         CALL GLYP(V1P,V2P,RP,RP13,Z0)
         V2P = DP5*Z0*V2P
         CALL ELYP(EM,RM,RM13,Z0)
         CALL GLYP(V1M,V2M,RM,RM13,Z0)
         V2M = DP5*Z0*V2M
         H10 = (EP - EM)/ADD2
         H20 = (V1P - V1M)/ADD2
         H11 = (V2P - V2M)/ADD2
C
         CALL ELYP(EP,R0,R013,ZP)
         CALL GLYP(V1P,V2P,R0,R013,ZP)
         V2P = DP5*ZP*V2P
         CALL ELYP(EM,R0,R013,ZM)
         CALL GLYP(V1M,V2M,R0,R013,ZM)
         V2M = DP5*ZM*V2M
         H01 = (EP - EM)/ADD2
         H02 = (V2P - V2M)/ADD2
C
         WRITE(LUPRI,'(A,3F15.10)')' LY010 test',G10,H10,G10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LY001 test',G01,H01,G01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP10 test',F10,H10,F10-H10
         WRITE(LUPRI,'(A,3F15.10)')' LYP01 test',F01,H01,F01-H01
         WRITE(LUPRI,'(A,3F15.10)')' LYP20 test',F20,H20,F20-H20
         WRITE(LUPRI,'(A,3F15.10)')' LYP11 test',F11,H11,F11-H11
         WRITE(LUPRI,'(A,3F15.10)')' LYP02 test',F02,H02,F02-H02
      END IF
      RETURN
      END
#endif   /* INCLUDE_DFT_OLD */
