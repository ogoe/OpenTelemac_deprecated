C                       ********************
                        LOGICAL FUNCTION EOF
C                       ********************
C
     *(LUNIT)
C
C***********************************************************************
C BIEF VERSION 5.1             17/08/94
C                              ORIGINE : ANTOINE YESSAYAN (MERCI TONIO)
C***********************************************************************
C
C     FONCTION  : DETECTION D'UNE FIN DE FICHIER
C
C                 SI EOF = .TRUE.  : FIN DE FICHIER
C                 SI EOF = .FALSE. : ON PEUT CONTINUER
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      LUNIT     | -->| UNITE LOGIQUE DU FICHIER QUE L'ON LIT
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: LUNIT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      EOF = .TRUE.
C
      READ ( UNIT=LUNIT , ERR=100 , END=100 )
C
      EOF = .FALSE.
C
100   CONTINUE
C
      BACKSPACE ( UNIT = LUNIT )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
