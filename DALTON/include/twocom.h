CTSAUE Modififed 240996
      LOGICAL         TCONAB,TCONCD,TPRIAB,TPRICD,
     &                DIAGAB,DIAGCD,DIACAB,DIACCD,
     &                SHAEQB,SHCEQD,SHABAB,GENAB, GENCD,
     &                RPRIAB, RPRICD, RCNTAB, RCNTCD,
     &                SPHRA, SPHRB, SPHRC, SPHRD,
     &                SPHRAB, SPHRCD, DIRAC, DOSCRN
      COMMON /TWOCOM/ CORAX0,CORAY0,CORAZ0,
     &                CORBX0,CORBY0,CORBZ0,
     &                CORCX0,CORCY0,CORCZ0,
     &                CORDX0,CORDY0,CORDZ0,
     &                ISHELA,ISHELB,ISHELC,ISHELD,
     &                ISYBLA,ISYBLB,ISYBLC,ISYBLD,
     &                NHKTA, NHKTB, NHKTC, NHKTD,
     &                KHKTA, KHKTB, KHKTC, KHKTD, KHKTAB, KHKTCD,
     &                KCKTA, KCKTB, KCKTC, KCKTD, KCKTAB, KCKTCD,KCKMAX,
     &                NCENTA,NCENTB,NCENTC,NCENTD,
     &                MULA,  MULB,  MULC,  MULD,
     &                MULTA, MULTB, MULTC, MULTD,
     &                NUCA,  NUCB,  NUCC,  NUCD,
     &                NUCAB, NUCCD, MXUCAB,MXUCCD,
     &                NORBA, NORBB, NORBC, NORBD, NORBAB,NORBCD,NOABCD,
     &                NSETA, NSETB, NSETC, NSETD,
     &                NSTRA, NSTRB, NSTRC, NSTRD,
     &                JMAXA, JMAXB, JMAXC, JMAXD,
     &                TCONAB,TCONCD,TPRIAB,TPRICD,
     &                DIAGAB,DIAGCD,DIACAB,DIACCD,
     &                GENAB, GENCD, IGENAB,IGENCD,
     &                SHAEQB,SHCEQD,SHABAB,
     &                MAXAB, MAXCD,
     &                IAB0X, IAB0Y, IAB0Z,
     &                ICD0X, ICD0Y, ICD0Z,
     &                NODCAB,NODCCD,
     &                NOPREP, IPTREP(0:8,2), NINTS, NINTSR(8), NINTMX,
     &                RPRIAB, RPRICD, RCNTAB, RCNTCD,
     &                SPHRA, SPHRB, SPHRC, SPHRD, SPHRAB, SPHRCD,
     &                IASMAX,IBSMAX,ICSMAX,IDSMAX,
     &                IASTRT,IBSTRT,ICSTRT,IDSTRT,
     &                DIRAC
      COMMON/SCREEN/SCRTHR,DNSMAX,ICEFLG,NEM,NCM,DOSCRN
      EQUIVALENCE (FCKTHR, SCRTHR)
