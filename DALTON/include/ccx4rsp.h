      INTEGER MAXX4LBL
      PARAMETER ( MAXX4LBL = 1 )

      LOGICAL LX4OPN
      INTEGER NX4LBL, ISYOFX4

      INTEGER ISYX4(MAXX4LBL,4)

      CHARACTER*8 LBLX4(MAXX4LBL,4)

      LOGICAL LORXX4(MAXX4LBL,4)

#if defined (SYS_CRAY)
      REAL FRQX4(MAXX4LBL,4)
#else
      DOUBLE PRECISION FRQX4(MAXX4LBL,4)
#endif

      COMMON/IX4RSP/ ISYX4, NX4LBL, ISYOFX4(8), LX4OPN, LORXX4
      COMMON/CX4RSP/ LBLX4
      COMMON/RX4RSP/ FRQX4

