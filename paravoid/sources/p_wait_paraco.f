C                       ************************
                        SUBROUTINE P_WAIT_PARACO
C                       ************************
C
     *(IBUF,NB)
C
C***********************************************************************
C  PARA       VERSION 5.9         23/06/2008        PASCAL VEZOLLE (IBM)
C  
C***********************************************************************
C
C      FONCTIONS: ATTENTE A LA FIN DE PARACO
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |                | -->|
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
      INTEGER IBUF(*), NB, IER
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_WAIT_PARACO VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_WAIT_PARACO IN VOID VERSION'
C
C----------------------------------------------------------------------
C
      RETURN
      END
