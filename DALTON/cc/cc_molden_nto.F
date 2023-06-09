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
C  /* Deck cc_molden */
      subroutine cc_molden_nto(work,lwork)
C
C     R. Faber et al., 2018
C
C     based on CC_MOLDEN by
C     Rolf H. Myhre and H. Koch
C
C     Purpose: Calculate natural transition orbitals and write them in
C     molden format to file
C
      implicit none

      character(len=*), parameter :: myname = 'CC_MOLDEN_NTO'
!
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccfop.h"
#include "ccsections.h"
#include "ccexgr.h"
#include "ccexci.h"
#include "ccexcinf.h"
#include "inftap.h"

!
      integer :: kcmot, knto, kucmo, koccu, krbtr, korbv
      integer :: kvec1, kumat, kvmat, ksval, kpos
      integer :: isyma, isymi, isym
      integer :: excisym, iexci
      integer :: mindim
      integer :: nto_count(8), ncount, nfound, ncount2
      integer :: i_info, idummy, ilist, ilistoff, iopt
      integer ludena
      integer :: icmo(8), iorbs(8)
      integer :: kend1, kend2, kend3
      integer :: lwrk1, lwrk2, lwrk3
      integer, intent(in) :: lwork
!
      character(len=20):: filename
      character(len=10):: model
      character(len=3), parameter :: list_type = 'RE '!'LE '
!
      double precision, intent(out):: work(lwork)
      double precision, parameter :: one = 1.d0, half = 0.5d0,
     &                               zero = 0d0,
     &                               thresh_nto = 1.0D-2
      double precision :: ddot, DNRM2
      double precision trace, dummy
!
      model ='CCSD      '
!
      nto_count = 0
!
      ncount  = 0
      ncount2 = 0
      do isym = 1, nsym
         icmo(isym) = ncount
         iorbs(isym) = ncount2
         !ncount = ncount + nbas(isym)*norb(isym)
         ncount = ncount + nbas(isym)*norbs(isym)
         ncount2 = ncount2 + norbs(isym)
      end do
!
!     Dynamic allocation
!     ------------------
!
      kcmot = 1
      knto = kcmot + nlamds
      kucmo = knto + nlamds
      krbtr = kucmo + nbast*norbts
      korbv = krbtr + nbast*nbast
      koccu = korbv + norbts
      kend1 = koccu + norbts
      lwrk1 = lwork - kend1 + 1
!      write(lupri,'(6I5)') kcmot, kucmo, krbtr, korbv, koccu, kend1
!
      if (lwrk1 .lt. 0) then
         stop 'Insufficient work space in '//myname
      endif

!
!     Close molden_MOS, we want new files
!     -----------------------------------
      call gpclose(LUMOLDEN_MOS,'KEEP')
!
!
!     Read in MO coefficients
!     -------------------------
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUERR)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMOT+I-1),I=1,nlamds)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C     Reorder and delete frozen orbitals
      CALL CMO_REORDER(WORK(KCMOT), WORK(KEND1), LWRK1)
!
!     Natural transition orbitals
!     ---------------------
!
      do excisym = 1, nsym

         ! symmetry specific entries
         kvec1 = kend1
         kend2 = kvec1 + nt1am(excisym)
         lwrk2 = lwork - kend2

         ilistoff = ISYOFE(excisym)

         do iexci = 1, nccexci(excisym,1) + nccexci(excisym,3)
            iopt = 1 ! Read only singles part of vector
            ilist = ilistoff + iexci
            ! Read in actual vector
            call cc_rdrsp(list_type, ilist, excisym, iopt,
     &                     model, work(kvec1), dummy)

            ! Zero the nto memory
            call dzero(work(knto), nlamds)
            call dzero(work(koccu), norbts)

            ! Loop over symmetry blocks of the vector
            do isymi = 1, nsym
               isyma = muld2h(isymi,excisym)

               mindim = min(nvir(isyma),nrhf(isymi))
               ! Skip if there is no work to do
               if (mindim.lt.1) cycle

               ! Allocate memory
               kumat = kend2
               kvmat = kumat + nvir(isyma)**2
               ksval = kvmat + nrhf(isymi)**2
               kend3 = ksval + mindim
               lwrk3 = lwork - kend3

               kpos = kvec1 + it1am(isyma,isymi)

               ! Do the SVD... Maybe we can use 'S' instead of 'A'
               CALL DGESVD(
     &            'A','A',nvir(isyma),nrhf(isymi),
     &            WORK(kpos),
     &            nvir(isyma), work(ksval),
     &            work(kumat), nvir(isyma),
     &            work(kvmat), nrhf(isymi),
     &            work(kend3), lwrk3, i_info
     &         )
               if (i_info .ne. 0) then
                  write(LUPRI,'(A,i5)') 'ERROR in DGSVD ', i_info
                  cycle
               endif


               ! work(kumat) now has the left/virtual eigenvectors on
               ! the columns
               ! work(kvmat) now has the right/occupied eigevectors on
               ! the rows

               ! Find the number of NTOs to write
               nfound = 0
               do n = 1, mindim
                  if (work(ksval-1+n) .lt. thresh_nto) exit
                  nfound = nfound + 1
               end do
               ! Ensure that we don't write too many, i.e. more than
               ! half the number ofbitals
               nfound = min(nfound, norbs(isymi)/2, norbs(isyma)/2)

               ! Transform the hole NTOs to AO basis
               call dgemm('N','T', nbas(isymi), nfound, nrhf(isymi),
     &                    one,
     &                    work(kcmot+ilmrhf(isymi)), nbas(isymi), ! AO x occ-MO matrix
     &                    work(kvmat), nrhf(isymi), ! NTO x occ-MO matrix
     &                    zero, work(knto+icmo(isymi)), nbas(isymi))

               ! Transform the particle NTOs to AO basis
               kpos = knto + icmo(isyma)
     &              + nbas(isyma)*(norbs(isyma)-nfound)
               call dgemm('N','N', nbas(isyma), nfound, nvir(isyma),
     &                    one,
     &                    work(kcmot+ilmvir(isyma)), nbas(isyma),
     &                    work(kumat), nvir(isyma),
     &                    zero, work(kpos), nbas(isyma))

               ! Put the hole NTO eigenvalues in the occupation vector
               kpos = koccu + iorbs(isymi)
               work(kpos:kpos+nfound-1) = - work(ksval:ksval+nfound-1)
               kpos = koccu + iorbs(isyma) + (norbs(isyma) - nfound)
               work(kpos:kpos+nfound-1) = work(ksval:ksval+nfound-1)

            end do
!
!           Create filename
!           ---------------
!
            write(filename,"(A5,I1,A,I0,A7)")
     &         "exci_", excisym, '_', iexci, ".molden"
!
!           Call molden_mos to write MOs in molden format
!           ---------------------------------------------
            call gpopen(LUMOLDEN_MOS,trim(filename),'NEW',' ',
     &                  'FORMATTED',idummy,.false.)
            rewind lumolden_mos
            call molden_mos(1, work(knto), work(koccu), work(krbtr),
     &                      work(kucmo),work(korbv))
            call gpclose(LUMOLDEN_MOS,'KEEP')
         end do
      end do
!
!     open original molden_mos.tmp
!     ----------------------------
      call gpopen(LUMOLDEN_MOS,'molden_MOS.tmp','UNKNOWN',' ',
     &            'FORMATTED',idummy,.false.)
!
      return
      end
