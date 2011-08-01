C                       ********************************
                        DOUBLE PRECISION FUNCTION P_DSUM
C                       ********************************
C
     *(MYPART)
C
C***********************************************************************
C  PARA VERSION 5.6         08/01/97        J-M HERVOUET (LNH)
C
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |  MYPART        | -->| CONTRIBUTION DU PROCESSEUR APPELANT.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      DOUBLE PRECISION MYPART
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_DSUM VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_DSUM IN ITS VOID VERSION'
C
      P_DSUM=MYPART
C
C-----------------------------------------------------------------------
C
      STOP
      END
 
