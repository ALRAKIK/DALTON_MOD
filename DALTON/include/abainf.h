!
!     File: abainf.h
!     Purpose: Control of what to do in ABACUS module
!
!     NOTE: DOSYM(NSYML) is sent to slaves in a parallel calculation and
!           should not be moved from this common block !
!
      INTEGER NSYML, IPRDEF, NWNABA, IRVIBG
      PARAMETER (NSYML = 8)
      LOGICAL ABAHF,  ABACI
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,                   &
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,                     &
     &        H2MO,   DOSYM(NSYML),   DOEXCI, SHIELD,                   &
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,                      &
     &        NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,                   &
     &        SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,                   &
     &        DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,                   &
     &        RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,                   &
     &        MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,                   &
     &        ABA_ALPHA,      EXPGRD, CTOCD,  NUMHES, DOD2DQ2,          &
     &        OECD,   MVEOR,  QPGRAD, SECNDM, THIRDM, VIB_G,            &
     &        GSDIP,  GSQUAD, GSOCT,  GSDIDI, GSDIQU, GSQUQU,           &
     &        VIBGIR, RNLRSCSHI, TDA_SINGLET, TDA_TRIPLET, RNLRSCEFG
      COMMON /ABAINF/ IPRDEF, NWNABA, IRVIBG,                           &
     &        ABAHF,  ABACI,                                            &
     &        MOLGRD, MOLHES, DIPDER, POLAR,  TSTINP,                   &
     &        VIB,    RESTAR, DOWALK, GDALL,  CCSD,                     &
     &        H2MO,   DOSYM,          DOEXCI, SHIELD,                   &
     &        SPNSPN, MAGSUS, VCD,    NACME,  AAT,                      &
     &        NOLOND, FCKDDR, ECD,    NODIFC, DODRCT,                   &
     &        SUPMAT, MOLGFA, OPTROT, SPINRO, MASSVE,                   &
     &        DARWIN, ABALNR, VROA,   NOCMC,  EXPFCK,                   &
     &        RAMAN,  QUADRU, NQCC,   HYPER,  VERDET,                   &
     &        MCD,    HELFEY, LINCPL, ABASOP, SKIPAB,                   &
     &        ABA_ALPHA,      EXPGRD, CTOCD,  NUMHES, DOD2DQ2,          &
     &        OECD,   MVEOR,  QPGRAD, SECNDM, THIRDM, VIB_G,            &
     &        GSDIP,  GSQUAD, GSOCT,  GSDIDI, GSDIQU, GSQUQU,           &
     &        VIBGIR, RNLRSCSHI, TDA_SINGLET, TDA_TRIPLET, RNLRSCEFG
! -- end of abainf.h --
