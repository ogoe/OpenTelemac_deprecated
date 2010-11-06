C                       ****************
                        SUBROUTINE ECRI2
C                       ****************
C
     *(X , I , C , NVAL , TYPE , CANAL , STD , ISTAT)
C
C***********************************************************************
C BIEF VERSION 5.1          17/08/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C FORMER SUBROUTINE ECRIT, NAME CHANGED BECAUSE THIS NAME EXISTS
C                          IN THE CALCIUM LIBRARY.
C
C FONCTION  :  ECRITURE DE VALEURS SUIVANT DIFFERENTS STANDARDS
C
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      X         | -->| TABLEAU A ECRIRE S'IL EST REEL
C |      I         | -->| TABLEAU A ECRIRE S'IL EST ENTIER
C |      C         | -->| CHAINE DE CARACTERES A ECRIRE
C |      NVAL      | -->| NOMBRE DE VALEURS DANS LE TABLEAU
C |                |    | OU NOMBRE DE CARACTERES DE LA CHAINE
C |      TYPE      | -->| TYPE DES DONNEES A ECRIRE :
C |                |    | 'I' , 'CH' , 'R4' , 'R8'
C |      CANAL     | -->| UNITE LOGIQUE POUR L'ECRITURE
C |      STD       | -->| STANDARD D'ECRITURE : STD , IBM OU I3E
C |      ISTAT     |<-- | ENTIER EN CAS D'ERREUR
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PRECAUTIONS D'EMPLOI : LES SOUS-PROGRAMMES ECRIBM ET ECRI3E SONT
C                         DEPENDANTS DE CHAQUE MACHINE UTILISEE. PAR
C                         EXEMPLE ECRIBM FIGURE DANS LA BIBLIOTHEQUE
C                         GENERALE DE IMA (EDF). SUR STATION, ECRIBM
C                         EST UN SOUS-PROGRAMME VIDE.
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NVAL,CANAL
      DOUBLE PRECISION, INTENT(IN) :: X(NVAL)
      INTEGER, INTENT(IN) :: I(NVAL)
      CHARACTER*(*), INTENT(IN) :: TYPE,STD,C
      INTEGER, INTENT(OUT) :: ISTAT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER J
C
      INTRINSIC REAL
C
C-----------------------------------------------------------------------
C
      ISTAT = 0
C
C-----------------------------------------------------------------------
C
      IF(STD(1:3).EQ.'STD') THEN
C
         IF (TYPE(1:2).EQ.'R4') THEN
            WRITE(CANAL)(REAL(X(J)),J=1,NVAL)
         ELSEIF (TYPE(1:2).EQ.'R8') THEN
            WRITE(CANAL)(X(J),J=1,NVAL)
         ELSEIF (TYPE(1:1).EQ.'I') THEN
            WRITE(CANAL)(I(J),J=1,NVAL)
         ELSEIF (TYPE(1:2).EQ.'CH') THEN
            WRITE(CANAL) C(1:NVAL)
         ELSE
            IF(LNG.EQ.1) THEN
              WRITE(LU,20) TYPE
20            FORMAT(1X,'ECRI2 : TYPE INCONNU :',A2)
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,21) TYPE
21            FORMAT(1X,'ECRI2 : UNKNOWN TYPE:',A2)
            ENDIF
            CALL PLANTE(0)
            STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
C     ELSEIF(STD(1:3).EQ.'IBM') THEN
C
C ATTENTION : ICI LA DOUBLE PRECISION CRAY N'EST PAS PREVUE
C
C        IF (TYPE(1:2).EQ.'R4') THEN
C           CALL ECRIBM( X , NVAL , TYPE , CANAL )
C        ELSEIF (TYPE(1:2).EQ.'R8') THEN
C           CALL ECRIBM( X , NVAL , TYPE , CANAL )
C        ELSEIF (TYPE(1:1).EQ.'I') THEN
C           CALL ECRIBM( I , NVAL , TYPE , CANAL )
C        ELSEIF (TYPE(1:2).EQ.'CH') THEN
C           CETTE RECOPIE SEMBLE EVITER UN BUG DE ECRIBM
C           CHAINE(1:NVAL) = C(1:NVAL)
C           CALL ECRIBM( CHAINE , NVAL , TYPE , CANAL )
C        ELSE
C           IF(LNG.EQ.1) WRITE(LU,20) TYPE
C           IF(LNG.EQ.2) WRITE(LU,21) TYPE
C           CALL PLANTE(0)
C           STOP
C        ENDIF
C
C-----------------------------------------------------------------------
C
C     ELSEIF(STD(1:3).EQ.'I3E') THEN
C
C ATTENTION : ICI LA DOUBLE PRECISION CRAY N'EST PAS PREVUE
C
C        IF (TYPE(1:2).EQ.'R4') THEN
C           CALL ECRI3E( X , NVAL , 'F' , CANAL , ISTAT )
C        ELSEIF (TYPE(1:2).EQ.'R8') THEN
C           CALL ECRI3E( X , NVAL , 'F' , CANAL , ISTAT )
C        ELSEIF (TYPE(1:1).EQ.'I') THEN
C           CALL ECRI3E( I , NVAL , 'I' , CANAL , ISTAT )
C        ELSEIF (TYPE(1:2).EQ.'CH') THEN
C           CALL ECRI3E( C(1:NVAL) , NVAL , 'C' , CANAL , ISTAT )
C        ELSE
C           IF(LNG.EQ.1) WRITE(LU,20) TYPE
C           IF(LNG.EQ.2) WRITE(LU,21) TYPE
C           CALL PLANTE(0)
C           STOP
C        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF(LNG.EQ.1) THEN
          WRITE(LU,10) STD
10        FORMAT(1X,'ECRI2 : STANDARD D''ECRITURE INCONNU :',A8)
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,11) STD
11        FORMAT(1X,'ECRI2 : UNKNOWN STANDARD:',A8)
        ENDIF
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
