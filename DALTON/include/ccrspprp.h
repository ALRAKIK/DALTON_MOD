! FILE: ccrspprp.h
      INTEGER MAXLBL_CC
      PARAMETER ( MAXLBL_CC = 1000 )
      INTEGER NPRLBL_CC
      CHARACTER*8 PRPLBL_CC
      COMMON/RSPPRPCC/ NPRLBL_CC
      COMMON/RSPLBLCC/ PRPLBL_CC(MAXLBL_CC)
! -- end of ccrspprp.h --