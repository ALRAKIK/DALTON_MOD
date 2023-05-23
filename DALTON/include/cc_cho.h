C Stuff for Cholesky CCSD(T)
      LOGICAL RSTVIR, SKIVI1,SKIVI2,
     &        UKNE4V, UKNE5V, UKNE4O, UKNE5O,
     &        RSTH, RSTH1, RSTF1, RSTC1, RSTC2
      COMMON /RSTCHT/ OLD4V, OLD5V, OLD4O, OLD5O,
     &                LUCHOU, IFVISY, IFVIOR, RSTVIR, SKIVI1,SKIVI2,
     &                UKNE4V, UKNE5V, UKNE4O, UKNE5O,
     &                RSTH, RSTH1, RSTF1, RSTC1, RSTC2
C
      COMMON /CHOXEN/ XENA,XENB,XENC,XEND,XENE,XENF,XENG,XENH,XENI,XENJ,
     &                XEN4,XEN5
      COMMON /CHOLE/ THRCHO, MXCHVE
      PARAMETER (MAXCHO = 100)
      COMMON /TIMINGS/ SCNDSA(MAXCHO),SCNDSB(MAXCHO),SCNDSC(MAXCHO),
     &                 SCNDSD(MAXCHO),SCNDSE(MAXCHO),SCNDSF(MAXCHO),
     &                 SCNDSG(MAXCHO),SCNDSH(MAXCHO),SCNDSI(MAXCHO),
     &                 SCNDSJ(MAXCHO),SCNDSX(MAXCHO)
      COMMON /ENRGCHO/ ENERGA(MAXCHO),ENERGB(MAXCHO),ENERGC(MAXCHO),
     &                 ENERGD(MAXCHO),ENERGE(MAXCHO),ENERGF(MAXCHO),
     &                 ENERGG(MAXCHO),ENERGH(MAXCHO),ENERGI(MAXCHO),
     &                 ENERGJ(MAXCHO),ENERGX(MAXCHO)

