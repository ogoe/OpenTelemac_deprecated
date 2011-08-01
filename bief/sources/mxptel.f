C                       *****************
                        SUBROUTINE MXPTEL
C                       *****************
C
     *(MXPTVS,MXELVS,IKLES,IELM,NPOIN,NELEM,NDP,IPOBO,LISTIN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C Pour le changement prevu, on va faire l'allocation de ITRAV ici en 
C interne, de toute facon il s'agit d'une variable de travail LOCALE.
C On peut passer en plus en argument le type de l'element pour faire
C le traitement specifique 3D.
CCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C***********************************************************************
C BIEF VERSION 5.1           24/08/95    J-M HERVOUET (LNH) 30 71 80 18
C
C***********************************************************************
C
C     CALCULE LE NOMBRE MAXIMUM DE POINTS ET D'ELEMENTS VOISINS D'UN
C     POINT POUR UN MAILLAGE DE TRIANGLES DONNE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   MXPTVS       |<-- | NOMBRE MAXIMUM DE POINTS VOISINS.
C |   MXELVS       |<-- | NOMBRE MAXIMUM D'ELEMENTS VOISINS.
C |   IKLES        | -->| TABLE DE CONNECTIVITE (DU FORMAT SELAFIN)
C |   ITRAV        | -->| TABLEAU DE TRAVAIL ENTIER DE DIMENSION NPOIN
C |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE.
C |   NELEM        | -->| NOMBRE D'ELEMENTS DU MAILLAGE.
C |   NDP          | -->| NOMBRE DE POINTS PAR ELEMENT.
C |   IPOBO        | -->| TABLEEAU QUI VAUT 0 POUR LES POINTS INTERIEURS
C |                |    | ET NON NUL POUR LES POINTS DE BORD.
C |                |    | POUR SE PLACER SUR LES ENREGISTREMENTS DES
C |   LISTIN       | -->| LOGIQUE : IMPRESSION DE MXELVS ET MXPTVS
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(INOUT) :: MXPTVS,MXELVS
      INTEGER, INTENT(IN)    :: IELM,NDP,NPOIN,NELEM
      INTEGER, INTENT(IN)    :: IKLES(NDP,NELEM),IPOBO(NPOIN)
      LOGICAL, INTENT(IN)    :: LISTIN
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C                       ITRAV(NPOIN) : TABLEAU AUTOMATIQUE
      INTEGER I,J,IELEM,ITRAV(NPOIN)
C
C-----------------------------------------------------------------------
C
C 1) INITIALISATION DU NOMBRE D'ELEMENTS VOISINS A 0 :
C
      DO 10 I = 1 , NPOIN
        ITRAV(I) = 0
10    CONTINUE
C
C 2) DECOMPTE DU NOMBRE D'ELEMENTS VOISINS PAR OPERATION D'ASSEMBLAGE :
C
      DO 22 J = 1, NDP
        DO 20 IELEM = 1 , NELEM
          ITRAV(IKLES(j,ielem)) = ITRAV(IKLES(j,ielem)) + 1
20      CONTINUE
22    CONTINUE
C
C 3) RECHERCHE DU MAXIMUM :
C
      MXELVS = ITRAV(1)
      DO 30 I = 2 , NPOIN
        MXELVS = MAX(MXELVS,ITRAV(I))
30    CONTINUE
C
C 4) NOMBRE DE POINTS VOISINS : IL FAUT AJOUTER 1 AU NOMBRE DE POINTS
C                               VOISINS POUR LES POINTS DE BORD.
C    ET RECHERCHE SIMULTANEE DU MAXIMUM
C      
      IF (IELM.EQ.31) THEN
        CALL MXPTEL31(NELEM,NPOIN,MXELVS,IKLES,MXPTVS)
      ELSE
        MXPTVS = MXELVS
        DO 40 I = 1 , NPOIN
          IF(IPOBO(I).NE.0) MXPTVS = MAX(MXPTVS,ITRAV(I)+1)
40      CONTINUE   
      ENDIF  
C
C-----------------------------------------------------------------------
C
      IF(LISTIN) THEN
        IF(LNG.EQ.1) WRITE(LU,97) MXELVS,MXPTVS
        IF(LNG.EQ.2) WRITE(LU,98) MXELVS,MXPTVS
      ENDIF
97    FORMAT(1X,'MXPTEL (BIEF) : NOMBRE MAXIMUM D''ELEMENTS VOISINS D''
     *UN POINT : ',1I3,/,1X,
     *          '                NOMBRE MAXIMUM DE POINTS VOISINS D''UN
     *POINT : ',1I3)
98    FORMAT(1X,'MXPTEL (BIEF) : MAXIMUM NUMBER OF ELEMENTS AROUND A POI
     *NT: ',1I3,/,1X,
     *          '                MAXIMUM NUMBER OF POINTS AROUND A POINT
     *: ',1I3)   
C
C-----------------------------------------------------------------------
C
      RETURN
      END
