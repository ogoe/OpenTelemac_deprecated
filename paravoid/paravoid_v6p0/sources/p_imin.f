C                       ***********************
                        INTEGER FUNCTION P_IMIN
C                       ***********************
C
     *(MYPART)
C
C***********************************************************************
C  PARA VERSION 5.1         08/01/97        J-M HERVOUET (LNH)
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
      INTEGER MYPART
C
C----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_IMAX VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_IMAX IN ITS VOID VERSION'
C
C----------------------------------------------------------------------
C
      P_IMIN=MYPART
C
      STOP
      END
 
