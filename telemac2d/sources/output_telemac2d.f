      SUBROUTINE OUTPUT_TELEMAC2D(TIME)
C-----------------------------------------------------------------------
C This file is part of TELEMAC-2D v5p7 
C-----------------------------------------------------------------------
C Written by JP Renaud
C
C This subroutine is just a "wrapper" for DESIMP so that outputs
C can be done from within ESTEL-3D when using the coupled model rather
C than relying on DESIMP (and its funny reliance on LT) directly.
C-----------------------------------------------------------------------
C
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C-----------------------------------------------------------------------
C  TIME: TIME TO PRINTOUT IN RECORD
C-----------------------------------------------------------------------
C ARGUMENTS
      DOUBLE PRECISION, INTENT(IN) :: TIME
C-----------------------------------------------------------------------
C LOCAL VARIABLES (required by DESIMP, no idea why...)
      DOUBLE PRECISION :: HIST(1)
      DATA HIST /9999.D0/
C-----------------------------------------------------------------------
C     PREPARATION OF THE RESULTS
      CALL PRERES_TELEMAC2D
C
C     OUTPUT A STANDARD TIME STEP
      CALL BIEF_DESIMP(T2D_FILES(T2DRES)%FMT,VARSOR,
     *                 HIST,0,NPOIN,T2D_FILES(T2DRES)%LU,'STD',
     *                 TIME,1,1,1,
     *                 SORLEO,SORIMP,MAXVAR,TEXTE,0,0)
C
C-----------------------------------------------------------------------
C
      RETURN
      END SUBROUTINE OUTPUT_TELEMAC2D
