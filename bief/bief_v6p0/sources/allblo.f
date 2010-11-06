C                       *****************
                        SUBROUTINE ALLBLO
C                       *****************
C
     *( BLO , NOM )
C
C***********************************************************************
C BIEF VERSION 5.5           10/01/95    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION  : ALLOCATION EN MEMOIRE D'UNE STRUCTURE DE BLOC
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   BLO          | -->| BLOC RESULTAT
C |   NOM          | -->| NOM FORTRAN DU TABLEAU
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_ALLBLO => ALLBLO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: BLO
      CHARACTER(LEN=6), INTENT(IN)    :: NOM
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ERR
C
C-----------------------------------------------------------------------
C  PARTIE COMMUNE A TOUS LES OBJETS
C-----------------------------------------------------------------------
C
C     MARQUAGE DE L'OBJET
C
      BLO%KEY = 123456
C
C     TYPE DE L'OBJET
C
      BLO%TYPE = 4
C
C     NOM DE L'OBJET
C
      BLO%NAME = NOM
C
C-----------------------------------------------------------------------
C  PARTIE SPECIFIQUE AUX BLOCS
C-----------------------------------------------------------------------
C
C     NOMBRE D'OBJETS DANS LE BLOC
C
      BLO%N = 0
C
C     ALLOCATION DU TABLEAU DE POINTEURS ADR
C
      BLO%MAXBLOCK = 128
      ALLOCATE(BLO%ADR(BLO%MAXBLOCK),STAT=ERR)
C
C-----------------------------------------------------------------------
C
      IF(ERR.EQ.0) THEN
C       IF(LNG.EQ.1) WRITE(LU,*) 'BLOC : ',NOM,' ALLOUE'
C       IF(LNG.EQ.2) WRITE(LU,*) 'BLOCK: ',NOM,' ALLOCATED'
      ELSE
        IF(LNG.EQ.1) WRITE(LU,10) NOM,ERR
        IF(LNG.EQ.2) WRITE(LU,20) NOM,ERR
10      FORMAT(1X,'ERREUR A L''ALLOCATION DU BLOC : ',A6,/,1X,
     *            'CODE D''ERREUR : ',1I6)
20      FORMAT(1X,'ERROR DURING ALLOCATION OF BLOCK: ',A6,/,1X,
     *            'ERROR CODE: ',1I6)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
