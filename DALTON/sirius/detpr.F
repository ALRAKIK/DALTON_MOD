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
C  /* Deck pr4inm */
      SUBROUTINE PR4INM(RHO,NDIM,IPACK,LUPRI)
C
C PRINT (PACKED) 4 INDEX MATRIX
C
C IPACK :
C         1 : MATRIX IS PACKED USING 8-FOLD PERMUTATIONAL SYMMETRY
C
#include "implicit.h"
      DIMENSION RHO(*)
      PARAMETER ( THRPRI = 1.0D-8 , THRZER = 1.0D-15 )
C     ... 880420: limit chosen based on the F13.8 format below
C
      IF (IPACK .EQ. 1 ) THEN
      WRITE(LUPRI,'(/A,1P,E10.2,A/)')
     *   ' ( only elements with absolut value greater than',THRPRI,' )'
C
      IPRPRI = 0
      IPRSML = 0
      IPRZER = 0
      IJKL = 0
      DO 100 I = 1, NDIM
        DO 90 J = 1, I
          DO 80 K = 1, I
            IF( K .EQ. I ) THEN
              LMAX = J
            ELSE
              LMAX = K
            END IF
C
            DO 70 L = 1, LMAX
              IJKL = IJKL + 1
              IF (ABS(RHO(IJKL)) .GT. THRPRI) THEN
                 WRITE(LUPRI,'(A,4I3,A,F13.8)')
     &           '  ',I,J,K,L,'...',RHO(IJKL)
                 IPRPRI = IPRPRI + 1
              ELSE IF (ABS(RHO(IJKL)) .GT. THRZER) THEN
                 IPRSML = IPRSML + 1
              ELSE
                 IPRZER = IPRZER + 1
              END IF
   70       CONTINUE
            WRITE(LUPRI,'()')
C
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
         WRITE(LUPRI,'(/I10,A,1P,E10.2,A/I10,A,E10.2,A/I10,A/)')
     *      IPRPRI,' elements larger than',THRPRI,' printed,',
     *      IPRSML,' small non-zero (',THRZER,') not printed, and',
     *      IPRZER,' zero elements not printed.'
C
      ELSE
         WRITE(LUPRI,'(/A,I4,A/)') ' PR4INM:  IPACK =',IPACK,
     *      ' is not implemented, nothing printed.'
      END IF
C
      RETURN
      END
C  /* Deck prwf */
      SUBROUTINE PRWF(C,XNDXCI,IOUT,THPWF,WORK,LFREE)
C
C  880528 hjaaj
C
!     module dependencies
      use lucita_mcscf_ci_cfg
      use lucita_mcscf_ci_interface_procedures
      use mcscf_or_gasci_2_define_cfg, only :
     &    define_lucita_cfg_dynamic

#include "implicit.h"
      DIMENSION C(*), XNDXCI(*), WORK(*)
C
C Used from common blocks:
C  INFINP : LSYM
C  INFORB : MULD2H
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infpri.h"
#include "ciinfo.h"
#include "priunit.h"
#include "dummy.h"
C
      integer :: xctype1, xctype2
      MAXPWF = 100
C
      IF ( NASHT .GT. 1 ) THEN
         IF (FLAG(27)) THEN
C           ... use determinant expansion
            ICSF = 0
         ELSE
            ICSF = 1
#if defined (VAR_OLDCODE)
#include "cinfo.h"
            WRITE (IOUT,'(/A/A)')
     *         ' Print of vector in CSF expansion not'//
     *         ' fully implemented yet.',
     *         ' (No decomposition in spin-couplings)'
            NCSF = NCSASM(LSYM)
            CALL PRVEC(NCSF,C,1,THPWF,MAXPWF,IOUT)
            CALL PRMGN(NCSF,C,1,12,IOUT)
#endif
         END IF
         if(ci_program .eq. 'LUCITA   ')then
           if(icsf .ne. 0) 
     &     call quit('*** prwf: lucita cannot handle csfs yet ***')
           kfree = 1
           lwrk  = kfree-1 + lfree

           xctype1 =  1
           xctype2 = -1

           call define_lucita_cfg_dynamic(lsym,        ! rhs symmetry
     &                                    lsym,        ! lhs symmetry
     &                                     0,          ! no istaci
     &                                     0,           ! spin for operator
     &                                     0,           ! spin for operator
     &                                     1,           ! one root
     &                                    ispin,       ! singlet, doublet, triplet, ...
     &                                    -1,          ! # of davidson iterations
     &                                    -1,          ! 1/2-particle density calcs
     &                                    iprsir,      ! print level      
     &                                    1.0d-10,     ! screening factor
     &                                    0.0d0,       ! davidson ci convergence
     &                                    0.0d0,       ! davidson ci convergence (auxiliary)
     &                                    .false.,     ! do not calculate nat. orb. occ. num. (internally inside LUCITA)
     &                                    .true.,      ! wave function analysis
     &                                    .false.,     ! no 2-el operators active
     &                                    .false.,     ! integrals provided by / return density matrices to the MCSCF environment
     &                                    .false.,     ! calculate 2-particle density matrix
     &                                    .false.,     ! calculate transition density matrix
     &                                    xctype1,     ! vector exchange type1
     &                                    xctype2,     ! vector exchange type2
     &                                    .false.,     ! vector exchange active in parallel runs in mc2lu interface (both mc-->lu and lu-->mc)
     &                                    .false.)     ! io2io vector exchange between mcscf and lucita: default mem2io/io2mem


           call mcscf_lucita_interface(c,vdummy,vdummy,vdummy,
     &                                 'analyze Cvec',
     &                                 work(kfree),lwrk,iprsir)
         else if(ci_program .eq. 'SIRIUS-CI')then
           CALL ANACIN(C,LSYM,THPWF,MAXPWF,XNDXCI,MULD2H,IOUT,ICSF)
C          CALL ANACIN(CIVEC,ICSYM,THRES,MAXTRM,XNDXCI,SYMPRO,IOUT,INCSFB)
         end if
      ELSE
         WRITE(IOUT,50)C(1)
   50    FORMAT(/'  Configuration no. 1 out of 1, coefficient',F10.6)
      END IF
C
      RETURN
      END
C  /* Deck anacin */
      SUBROUTINE ANACIN(CIVEC,ICSYM,THRES,MAXTRM,XNDXCI,SYMPRO,IOUT,
     &                 INCSFB)
C
C OUTER ROUTINE FOR CI ANALYZER
C
C Revision 910812-hjaaj: transfer ILTSOB to ANACI2 and ANACSF
C
#include "implicit.h"
C
      DIMENSION CIVEC(*),XNDXCI(*)
      INTEGER   SYMPRO(8,8)
C
#include "mxpdim.h"
#include "ciinfo.h"
#include "detbas.h"
#include "csfbas.h"
#include "strnum.h"
C
c remco
      real(8), allocatable :: detvec(:), cvec(:)
#include "maxorb.h"
#include "infinp.h"
#include "infvar.h"
#include "priunit.h"
c einde remco
      IF(INCSFB.EQ.0) THEN
C. SD basis
        IF( ICOMBI .EQ. 0 ) THEN
          NCIVAR = NDTASM(ICSYM)
        ELSE
          NCIVAR = NCMASM(ICSYM)
        END IF
          CALL ANACI2(CIVEC,XNDXCI(KNSSOA),XNDXCI(KNSSOB),NOCTPA,NOCTPB,
     &                MAXSYM,SYMPRO,THRES,MAXTRM,XNDXCI(KISSOA),
     &                XNDXCI(KISSOB),NAEL,NBEL,XNDXCI(KIOCOC),ICSYM,
     &                NCIVAR,ICOMBI,XNDXCI(KIASTR),XNDXCI(KIBSTR),
     &                XNDXCI(KLTSOB),IOUT)
      ELSE ! CSF basis (INCSFB .ne. 0)

        CALL ANACSF(CIVEC,XNDXCI(KICONF(1)),XNDXCI(KCFTP),ICSYM,
     &              NAEL+NBEL,THRES,MAXTRM,XNDXCI(KLTSOB),IOUT)

        if (DONEVPT) then
           write(lupri,100)
           allocate(detvec(ncdets))
           allocate(cvec(ncdets))
           call dcopy(nconf,civec,1,cvec,1)
           call csdtvc(cvec,detvec,1,XNDXCI(KDTOC),
     *                XNDXCI(KICTS(1)),lsym,0,0)
           IF( ICOMBI .EQ. 0 ) THEN
              NCIVAR = NDTASM(ICSYM)
           ELSE
              NCIVAR = NCMASM(ICSYM)
           END IF
           CALL ANACI2(detvec,XNDXCI(KNSSOA),XNDXCI(KNSSOB),NOCTPA,
     &              NOCTPB,MAXSYM,SYMPRO,THRES,MAXTRM,XNDXCI(KISSOA),
     &              XNDXCI(KISSOB),NAEL,NBEL,XNDXCI(KIOCOC),ICSYM,
     &              NCIVAR,ICOMBI,XNDXCI(KIASTR),XNDXCI(KIBSTR),
     &              XNDXCI(KLTSOB),IOUT)
           deallocate (detvec)
           deallocate (cvec)
        endif
      END IF
C
      RETURN
 100  format(//'   *** Calculation was done in CSF basis, but results ar
     *e printed in determinant basis for NEVPT2 calculation ***',//)
      END
C  /* Deck anacsf */
      SUBROUTINE ANACSF(CIVEC,ICONF,IPROCS,IREFSM,NEL,THRES,MAXTRM,
     &                  ILTSOB,IOUT)
C
C Analyze CI vector in CSF basis
C
! hjaaj Nov 07: TODO, make output something like:
! (is configuration number not irrelevant?)
!
!   Symmetry      1---------- 2-------
!   Orbital       1   2   3   4   5   6
! -0.06584956     2   2       2                   <- spin coupling
! -0.06242366     2       2       2
! -0.07355010     2   2   1  -1   1  -1
C
#include "implicit.h"
#include "mxpdim.h"
#include "strnum.h"
#include "spinfo.h"
#include "ciinfo.h"
      DIMENSION ICONF(*),IPROCS(*),ILTSOB(*)
      DIMENSION CIVEC(*)
C
      NCIVAR = NCSASM(IREFSM)
      VNRM = DNRM2(NCIVAR,CIVEC,1)
      ITRM = 0
      ILOOP = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
C
2001  CONTINUE
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = VNRM
        XMIN = VNRM/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
C
      WRITE(IOUT,'(/A,1P,E12.4,A,E12.4/A)')
     &' Printout of coefficients in interval',XMIN,' to',XMAX,
     &' ==============================================================='
C. Loop over configurations and CSF's for given configuration
      ICNF = 0
      ICSF = 0
      INRANG = 0
      NOPRT = 0
      ICNBS0 = 1
      IPBAS = 1
      DO 1000 ITYP = 1, NTYP
        IOPEN = ITYP + MINOP - 1
        ICL = (NEL - IOPEN) / 2
        IOCC = IOPEN + ICL
        IF( ITYP .GT. 1 ) THEN
          ICNBS0 = ICNBS0 + NCNFTP(ITYP-1,IREFSM)*(NEL+IOPEN-1)/2
          IPBAS = IPBAS + NCSFTP(ITYP-1)*(IOPEN-1)
        END IF
C. Configurations of this type
        NNCNF = NCNFTP(ITYP,IREFSM)
        NNCSF = NCSFTP(ITYP)
        DO 900  IC = 1, NNCNF
          ICNF = ICNF + 1
          ICNBS = ICNBS0 + (IC-1)*(IOPEN+ICL)
C. CSF's in this configuration
          DO 800 IICSF = 1, NNCSF
            ICSF = ICSF+1
            IF( XMAX .GE. ABS(CIVEC(ICSF)) .AND.
     &          XMIN .LT. ABS(CIVEC(ICSF)) ) THEN
              ITRM  = ITRM + 1
              INRANG = INRANG + 1
              IF( ITRM .LE. MAXTRM ) THEN
                CNORM = CNORM + CIVEC(ICSF) ** 2
                WRITE(IOUT,'(/A,I10,A,F16.8,1P,E18.8)')
     &           '  Coefficient of CSF no.',ICSF,' is',
     &           CIVEC(ICSF),CIVEC(ICSF)
                WRITE(IOUT,'(A,40I4)') '  Orbital       ',
     &           (ILTSOB(ICONF(ICNBS-1+II)),II = 1,ICL+IOPEN)
                WRITE(IOUT,'(A,40I4)') '  Spin coupling ',
     &           ( 2, II = 1, ICL),
     &           ((2*IPROCS(IPBAS+(IICSF-1)*IOPEN-1+II)-1),II=1,IOPEN)
              ELSE
                NOPRT = NOPRT + 1
              END IF
            END IF
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
      IF(INRANG .EQ. 0 ) WRITE(IOUT,'(A)') '   ( no coefficients )'
      IF(NOPRT.NE.0) WRITE(IOUT,'(A,I10)')
     & ' Number of coefficients not printed', NOPRT
      IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
       WRITE(IOUT,'(/A,F15.8)')
     & ' Norm of printed CI vector .. ', CNORM
C
       WRITE(IOUT,'(/A/A/A,1P,E10.2,A/)')
     &    ' Magnitude of CI coefficients',
     &    ' ============================',
     &    ' ( Ranges are relative to norm of vector :',VNRM,' )'
C
       CNORM = 0.0D0
       ISUM = 0
       XMIN = VNRM + 0.1D0
       DO 200 IPOT = 0, 10
         CLNORM = 0.0D0
         INRANG = 0
         XMAX = XMIN
         XMIN = VNRM * (0.1D0 ** (IPOT+1))
C
         DO 180 IDET = 1, ICSF
           IF( ABS(CIVEC(IDET)) .LE. XMAX  .AND.
     &         ABS(CIVEC(IDET)) .GT. XMIN ) THEN
                 INRANG = INRANG + 1
                 CLNORM = CLNORM + CIVEC(IDET) ** 2
           END IF
  180    CONTINUE
         CNORM = CNORM + CLNORM
C
         IF (INRANG .GT. 0)
     &     WRITE(IOUT,'(A,I2,A,I2,I10,E18.8,F18.10)')
     &     '  10-',IPOT+1,' to 10-',IPOT,INRANG,CLNORM,CNORM
C
         ISUM = ISUM + INRANG
  200 CONTINUE
C
      WRITE(IOUT,*) ' Number of coefficients less than 10^-11 times ',
     & 'norm is ',ICSF - ISUM
C
      RETURN
      END
C  /* Deck anaci2 */
      SUBROUTINE ANACI2(CIVEC,NSSOA,NSSOB,NOCTPA,NOCTPB,MAXSYM,SYMPRO,
     &                 THRES,MAXTRM,ISSOA,ISSOB,NAEL,NBEL,IOCOC,
     &                 ICSM,NCIVAR,ICOMBI,IASTR,IBSTR,ILTSOB,IOUT)
C
C ANALYZE CI VECTOR :
C       PRINT atmost MAXTRM DETERMINANTS OR COMBINATIONS
C       with COEFFICIENTS LARGER THAN THRES.
C
C       NUMBER OF COEFFICIENTS OF EACH MAGNITUDE
C
C.. RAS VERSION
C
C JEPPE OLSEN JANUARY 1989
c
c Modified Spring 1992 to allow vector with general norm
C
#include "implicit.h"
#include "dummy.h"
#include "maxorb.h"
#include "infinp.h"
#include "infpt2.h"
      INTEGER SYMPRO(8,8)
      DIMENSION CIVEC(*)
      DIMENSION NSSOA(NOCTPA,MAXSYM),NSSOB(NOCTPB,MAXSYM)
      DIMENSION ISSOA(NOCTPA,MAXSYM),ISSOB(NOCTPB,MAXSYM)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*),ILTSOB(*)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      PARAMETER (LEN_LINE = 300)
      CHARACTER*(LEN_LINE) ALINE, BLINE
      data itimes/0/
      save itimes
      logical owrite
C
      itimes=itimes+1
      owrite=.false.
      if (doci) then
          if (itimes.eq.istnevci) owrite=.true.
      else
          owrite=.true.
      endif
      LUN98=-13300
      if (owrite) CALL GPOPEN(
     *   LUN98,'DETSTR','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)


      VNRM = DNRM2(NCIVAR,CIVEC,1)
C
C     Warn if M_S = 0 and triplet combination of determinants
C     (In determinant basis, triplet states for M_S = 0 have coefficients with
C      same absolute value but opposite signs when switching alpha- and beta-strings.
C      For singlet states they have same sign.)
C
      IF (NAEL .EQ. NBEL .AND. ICOMBI .EQ. 0) THEN
      IF (ABS(1.0D0-VNRM) .LT. 1.D-14) THEN ! only print for normalized CI vectors
         XMIN = 0.0D0
         XMAX = 0.0D0
         DO I = 1, NCIVAR
            XMAX = MAX(XMAX, CIVEC(I))
            XMIN = MIN(XMIN, CIVEC(I))
         END DO
         IF (ABS(XMAX+XMIN) .LT. 1.D-14) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(//A/)')
     &      '@ WARNING : This M_S=0 CI determinant vector is '//
     &         'not a singlet state (S=0)',
     &      '@ WARNING : It is probably a triplet state (S=1) '//
     &         'but it might also have S=3, or S=5, or ...',
     &      '@ WARNING : If you wanted a singlet state, '//
     &         'then use CSFs or specify .PLUS COMBINATIONS'
         END IF
      END IF
      END IF
c
      ITRM = 0
      IDET = 0
      IIDET = 0
      ILOOP = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
2001  CONTINUE
      IIDET = 0
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = VNRM
        XMIN = VNRM/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IDET = 0
C
      WRITE(IOUT,'(//A,E11.4,A,E11.4/A)')
     &'  Printout of coefficients in interval  ',XMIN,' to ',XMAX,
     &'  =============================================================='
      NOPRT = 0
      DO 100 IASM = 1, MAXSYM
        IBSM = SYMPRO(ICSM,IASM)
        DO 95 IATP = 1, NOCTPA
        DO 94 IBTP = 1, NOCTPB
        IF( IOCOC(IATP,IBTP) .LE. 0 ) GOTO 94
C
        IABAS = ISSOA(IATP,IASM)
        IBBAS = ISSOB(IBTP,IBSM)
C
        NIA = NSSOA(IATP,IASM)
        NIB = NSSOB(IBTP,IBSM)
C
        DO 90 IB = IBBAS,IBBAS+NIB-1
          DO 80 IA = IABAS,IABAS+NIA-1
            IF(IB.GT.IA .AND. ICOMBI.NE.0)GOTO 80
C
            IDET = IDET + 1
c--renzo begin
c            print*,'idet,civec',idet,civec(idet)
            if(iloop.eq.1.and.owrite)then
               write(LUN98) civec(idet),
     &              (ILTSOB(IASTR(IEL,IA)), IEL=1,NAEL),
     &              (ILTSOB(IBSTR(IEL,IB)), IEL=1,NBEL)
            endif
c--renzo end
            IF( XMAX .GE. ABS(CIVEC(IDET)) .AND.
     &          XMIN .LT. ABS(CIVEC(IDET))) THEN
              ITRM = ITRM + 1
              IIDET = IIDET + 1
              IF( ITRM .LE. MAXTRM ) THEN
                CNORM = CNORM + CIVEC(IDET) ** 2
C
                IF(ICOMBI .EQ. 0 )THEN
                  WRITE(IOUT,'(/A,I10,A,F16.8,1P,E16.8)')
     &            ' Coefficient of determinant',IDET,' is',
     &            CIVEC(IDET),CIVEC(IDET)
                ELSE
                  WRITE(IOUT,'(/A,I10,A,F16.8,1P,E16.8)')
     &            ' Coefficient of combination',IDET,' is',
     &            CIVEC(IDET),CIVEC(IDET)
                END IF
C
                aline = ' '
                do iel = 1, nael
                   iael = iltsob(iastr(iel,ia))
                   if (iael .le. (LEN_LINE/3) ) then
                      jael = 3*iael - 2  ! from: jael = 1 + 3*(iael-1)
                      write(aline(jael:),'(I3)') iael
                   else
                      write (iout,*) 'Overflow alpha orbital:', iael
                   end if
                end do
                bline = ' '
                do iel = 1, nbel
                   ibel = iltsob(ibstr(iel,ib))
                   if (ibel .le. (LEN_LINE/3) ) then
                      jbel = 3*ibel - 2  ! from: jbel = 1 + 3*(ibel-1)
                      write(bline(jbel:),'(I3)') ibel
                   else
                      write (iout,*) 'Overflow  beta orbital:', ibel
                   end if
                end do
                iend = LNBLNK(aline)
                write(iout,'(2A)') ' alpha-string:',aline(1:iend)
                iend = LNBLNK(bline)
                write(iout,'(2A)') '  beta-string:',bline(1:iend)
c               WRITE(IOUT,'(A,20I4:/,(14X,20I4))')
c    &          ' alpha-string:',(ILTSOB(IASTR(IEL,IA)),IEL = 1, NAEL )
c               WRITE(IOUT,'(A,20I4:/,(14X,20I4))')
c    &          '  beta-string:',(ILTSOB(IBSTR(IEL,IB)),IEL = 1, NBEL )
              ELSE
                NOPRT = NOPRT + 1
              END IF
            END IF
   80     CONTINUE
   90   CONTINUE
   94   CONTINUE
   95  CONTINUE
  100 CONTINUE
       IF(IIDET .EQ. 0 ) WRITE(IOUT,'(A)') '   ( no coefficients )'
       IF(NOPRT .NE. 0 ) WRITE(IOUT,'(A,I10)')
     & ' Number of coefficients not printed', NOPRT
       IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
C
       WRITE(IOUT,'(/A,F15.8)')
     & ' Norm of printed part of CI vector .. ', CNORM
C
C.. SIZE OF CI-COEFFICIENTS
C
       WRITE(IOUT,'(/A)') '  Magnitude of CI coefficients'
       WRITE(IOUT,'(A/)') '  ============================'
       WRITE(IOUT,'(A,1P,E10.2,A/)')
     & '  ( Ranges are relative to norm of vector :',VNRM,' )'
C
       CNORM = 0.0D0
       ISUM = 0
       XMIN = VNRM + 0.1D0
C      ... add 0.1 to eliminate round-off problems if only one element
       DO 200 IPOT = 0, 10
         CLNORM = 0.0D0
         INRANG = 0
         XMAX = XMIN
         XMIN = VNRM * (0.1D0 ** (IPOT+1))
C
         DO 180 IDET = 1, NCIVAR
           IF( ABS(CIVEC(IDET)) .LE. XMAX  .AND.
     &         ABS(CIVEC(IDET)) .GT. XMIN ) THEN
                 INRANG = INRANG + 1
                 CLNORM = CLNORM + CIVEC(IDET) ** 2
           END IF
  180    CONTINUE
         CNORM = CNORM + CLNORM
C
         IF (INRANG .GT. 0)
     &     WRITE(IOUT,'(A,I2,A,I2,3X,I7,3X,E15.8,3X,E15.8)')
     &     '  10-',IPOT+1,' to 10-',IPOT,INRANG,CLNORM,CNORM
C
         ISUM = ISUM + INRANG
  200 CONTINUE
C
      WRITE(IOUT,*) ' Number of coefficients less than 10^-11 times ',
     & 'norm is  ',NCIVAR - ISUM
C
      IF (owrite) CALL GPCLOSE(LUN98,'KEEP')
      RETURN
      END
