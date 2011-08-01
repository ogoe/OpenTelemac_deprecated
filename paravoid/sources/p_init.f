C                       *****************
                        SUBROUTINE P_INIT
C                       *****************
C
     *(CHAINE,NCAR,IPID,NCSIZE)
C
C***********************************************************************
C  PARA VERSION 5.9       /06/96           HANS HERRMANN (HANOVRE)
C             MODIFIED        08/06/96     REINHARD HINKELMANN (HANOVRE)
C             MODIFIED        17/12/96            J-M HERVOUET (LNH)
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
C APPELE PAR : PRINCI
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER IPID,NCSIZE
C
      CHARACTER*144 CHAINE
      INTEGER NCAR
C
      NCAR = 0
      CHAINE =' '
      IPID=0
C
C     IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_INIT VERSION VIDE'
C     IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_INIT IN ITS VOID VERSION'
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
