      LOGICAL ESRHFC, HFCFC, HFCSD, HFCSO
      LOGICAL HSOBRT, HSOEFF, UNGAUSS
      INTEGER IPRHFC
      DOUBLE PRECISION AFCEXP, AFCPOL, ASDEXP, ASDPOL, ASOVAL
      COMMON /RSPHFC/ ASDEXP(MXCENT,3,3), ASDPOL(MXCENT,3,3),
     &       ASOVAL(MXCENT,3,3), AFCEXP(MXCENT), AFCPOL(MXCENT),
     &       IPRHFC, ESRHFC, HFCFC, HFCSD, HFCSO, HSOBRT, HSOEFF,
     &       UNGAUSS
