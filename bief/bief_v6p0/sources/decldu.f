C                       *****************
                        SUBROUTINE DECLDU
C                       *****************
C
     *(B,A,MESH,COPY,LV)
C
C***********************************************************************
C BIEF VERSION 5.5           26/02/04    J-M HERVOUET (LNH) 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : DECOMPOSITION L D U DES MATRICES ELEMENTAIRES CONTENUES
C            DANS LA MATRICE A.
C
C            ON EXIGE ICI QUE LA DIAGONALE DE A SOIT L'IDENTITE
C
C            CHAQUE MATRICE ELEMENTAIRE EST DECOMPOSEE SOUS LA FORME :
C
C            LE X DE X UE
C
C            LE : TRIANGULAIRE INFERIEURE AVEC DES 1 SUR LA DIAGONALE.
C            DE : DIAGONALE
C            UE : TRIANGULAIRE SUPERIEURE AVEC DES 1 SUR LA DIAGONALE.
C
C                                                T
C            SI LA MATRICE EST SYMETRIQUE : LE =  UE
C
C            LES MATRICES "DE" SONT CONSIDEREES COMME DES DIAGONALES
C            DE TAILLE NPOIN X NPOIN QU'IL FAUT IMAGINER COMPLETEES
C            AVEC DES 1 POUR LES POINTS QUI N'APPARTIENNENT PAS A
C            L'ELEMENT CONSIDERE.
C
C            ON EFFECTUE ENSUITE LE PRODUIT DE TOUTES CES DIAGONALES
C            QUI DONNE LA DIAGONALE DB.
C
C   ||||||   ENFIN : DB EST INVERSEE CAR C'EST SOUS CETTE FORME QU'ELLE
C                    SERA UTILISEE DANS DESREM.
C
C-----------------------------------------------------------------------
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  21 : QUADRILATERE Q1        4                       OUI
C  41 : PRISMES TELEMAC-3D     6                       OUI
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      B         |<-- |  RESULTAT
C |      A         | -->|  MATRICE A
C |      MESH      | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |      COPY      | -->|  SI .TRUE. A EST COPIEE DANS B
C |                |    |  SINON ON CONSIDERE QUE B EST DEJA REMPLIE
C |      LV        | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : DLDU11 , DLDU21 , DLDU41 , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_DECLDU => DECLDU
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ) , INTENT(INOUT) :: B
      TYPE(BIEF_OBJ) , INTENT(IN)    :: A
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      LOGICAL        , INTENT(IN)    :: COPY
      INTEGER        , INTENT(IN)    :: LV
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NELMAX,IELM,NPOIN,NELEM
C
      CHARACTER*1 TYPDA,TYPEA
C
C-----------------------------------------------------------------------
C
      IELM   = A%ELMLIN
      NELEM  = MESH%NELEM
      NELMAX = MESH%NELMAX
C
      TYPDA = A%TYPDIA
      TYPEA = A%TYPEXT
C
      NPOIN = A%D%DIM1
C
C-----------------------------------------------------------------------
C
      IF(A%STO.EQ.1) THEN
C
      IF(IELM.EQ.11) THEN
C
        CALL DLDU11(B%D%R,B%X%R,TYPDA,A%X%R,TYPEA,
     *              MESH%IKLE%I,NELEM,NELMAX,NPOIN,MESH%W%R,COPY,LV)
C
      ELSEIF(IELM.EQ.21.OR.IELM.EQ.12.OR.IELM.EQ.31.OR.IELM.EQ.51) THEN
C
        CALL DLDU21(B%D%R,B%X%R,TYPDA,A%X%R,TYPEA,
     *              MESH%IKLE%I,NELEM,NELMAX,NPOIN,MESH%W%R,COPY,LV)
C
      ELSEIF(IELM.EQ.41) THEN
C
        CALL DLDU41(B%D%R,B%X%R,TYPDA,A%X%R,TYPEA,
     *              MESH%IKLE%I,NELEM,NELMAX,NPOIN,MESH%W%R,COPY,LV)
C
C  IELM NON PREVU : ERREUR
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELM
       IF (LNG.EQ.2) WRITE(LU,101) IELM
100    FORMAT(1X,'DECLDU (BIEF) : IELM = ',1I6,' ELEMENT NON PREVU')
101    FORMAT(1X,'DECLDU (BIEF) : IELM = ',1I6,' ELEMENT NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
C
      ELSEIF(A%STO.EQ.3) THEN
        CALL DLDUSEG(B%D%R,B%X%R,TYPDA,A%X%R,TYPEA,
     *               MESH%GLOSEG%I,MESH%NSEG,NPOIN,COPY)
      ELSE
        WRITE(LU,*) 'UNKNOWN MATRIX STORAGE IN DECLDU'
        CALL PLANTE(1)
        STOP 
      ENDIF
C
C-----------------------------------------------------------------------
C
C  DESCRIPTION DE B
C
      B%TYPDIA='Q'
      B%TYPEXT=TYPEA
C
C-----------------------------------------------------------------------
C
      RETURN
      END
