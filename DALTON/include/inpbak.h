! File:  inpbak.h
!
!     The following common blocks are used for backup of
!     sirius variables in infinp.h
!
      COMMON /INTBAK/ ISPINB,ISTATB,LSYMB,NACTLB, MCTYPB,
     &                LSLMXB,NLMSLB,NLMN1B,NLMX1B,NLMN3B,NLMX3B,
     &                LROTSB,NROTSB,IROOTB(MAXRTS),
     &                NOROTB(MXCORB),IMORDB(MXCORB),
     &                IORTOB,ICI0B ,KDELB ,ICHCKB,NTITB,
     &                MXMACB,MXMICB,MAXJTB,MXCITB,MXUITB,MXAPMB,
     &                MXABSB,ITRLVB,ITRFNB,JCHSMB,JCHRBB,
     &                NROCIB,ISTCIB,NFILDB,MXCMAB,ICCNOB,IMCNOB,
     &                NMCAVB
      LOGICAL  FLAGB(NFLAG), DOSCFB,DOMP2B,DOCNOB,DOCIB, DOMCB,
     &         DORSPB,FCVRBB,LNROTB,LMORDB,DRFCKB,CORHLB,CRRLXB,
     &         RSPHPB,JOLSNB,       INRSIB,INRSFB,SUPSMB,PRSIRB,
     &         PCMB
      COMMON /LOGBAK/  FLAGB,DOSCFB,DOMP2B,DOCNOB,DOCIB, DOMCB,
     &         DORSPB,FCVRBB,LNROTB,LMORDB,DRFCKB,CORHLB,CRRLXB,
     &         RSPHPB,JOLSNB,       INRSIB,INRSFB,SUPSMB,PRSIRB,
     &         PCMB
      COMMON /RELBAK/  SPINB,POTNCB,EPSOLB,EPSTTB,RSOLB(3),
     &                TRGRDB,TRPWFB,TRCIB,TRMCB,TRCGRB,
     &                EFILDB(MXFELT), CMXMOB, TROVLB,
     &                TRSSYB
      CHARACTER*60  TITLEB(6)
      CHARACTER*200 TITMLB(2)
      CHARACTER*4   CENTB, TYPEB
      CHARACTER*8   LFILDB
      COMMON /CHRBAK/ TITLEB, TITMLB,
     &                CENTB(MXCORB), TYPEB(MXCORB), LFILDB(MXFELT)
! -- end of inpbak.h --
