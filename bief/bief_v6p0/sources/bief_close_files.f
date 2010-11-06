C                       ***************************
                        SUBROUTINE BIEF_CLOSE_FILES
C                       ***************************
C
     &(CODE,FILES,NFILES,PEXIT)
C
C***********************************************************************
C BIEF VERSION 6.0     01/04/2009     J-M HERVOUET (LNHE) 01 30 87 80 18
C 
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                | -->|
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : HOMERE
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      USE BIEF, EX_BIEF_CLOSE_FILES => BIEF_CLOSE_FILES
C
      USE DECLARATIONS_TELEMAC
      USE M_MED
C
      IMPLICIT NONE
      INTEGER     LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          , INTENT(IN)     :: NFILES
      CHARACTER(LEN=24), INTENT(IN)     :: CODE
      LOGICAL, INTENT(IN)               :: PEXIT
      TYPE(BIEF_FILE)   , INTENT(INOUT) :: FILES(NFILES)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      DO I=1,NFILES
C
        IF(FILES(I)%NAME(1:1).NE.' ') THEN
C 
C         FERMETURE DU FICHIER
C
          IF(FILES(I)%FMT.EQ.'MED     ') THEN
            CALL CLOSE_FILE_MED(FILES(I)%LU)
          ELSE
            CLOSE(FILES(I)%LU)
          ENDIF
C
        ENDIF
C
      ENDDO
C
C-----------------------------------------------------------------------
C
C     PARALLELISM: STOPPING IF PEXIT
C
      IF(NCSIZE.GT.0.AND.PEXIT) CALL P_EXIT
C
C-----------------------------------------------------------------------
C
      RETURN
      END
