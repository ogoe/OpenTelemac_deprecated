C                       ********************
                        SUBROUTINE OIL_SPILL
C                       ********************
C
C
C***********************************************************************
C TELEMAC 2D VERSION 6.0         20/04/2010         CEDRIC GOEURY (LHSV)
C***********************************************************************
C
C  FUNCTION  : OIL SPILL MODEL
C
C              CALLED IF KEYWORD 'OIL SPILL MODEL' SET TO YES
C              AND USES 'MIGRHYCAR STEERING FILE'
C
C              LOGICAL UNIT OF MIGRHYCAR STEERING FILE IS:
C              T2D_FILES(T2DMIG)%LU
C
C 
C-----------------------------------------------------------------------
C  EXAMPLES OF DATA IN DECLARATIONS_TELEMAC2D 
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |      AT        | -->| CURRENT TIME
C |      LT        | -->| ITERATION NUMBER
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT :
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
!     DOUBLE PRECISION 
!     LOGICAL 
!     INTEGER
!     SAVE 
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
      RETURN
      END
