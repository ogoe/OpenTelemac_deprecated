C                       *****************
                        SUBROUTINE DLDU21
C                       *****************
C
     *(DB,XB,TYPDIA,XA,TYPEXA,
     * IKLE,NELEM,NELMAX,NPOIN,W,COPY,LV)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : POUR LES QUADRILATERES Q1
C
C            DECOMPOSITION L D U DES MATRICES ELEMENTAIRES CONTENUES
C            DANS LA MATRICE A.
C
C            ON EXIGE ICI QUE LA DIAGONALE DE A SOIT L'IDENTITE
C
C            CHAQUE MATRICE ELEMENTAIRE EST DECOMPOSEE SOUS LA FORME :
C
C            LE * DE * UE
C
C            LE : TRIANGULAIRE INFERIEURE AVEC DES 1 SUR LA DIAGONALE.
C            DE : DIAGONALE
C            UE : TRIANGULAIRE SUPERIEURE AVEC DES 1 SUR LA DIAGONALE.
C
C                                                   T
C            SI LA MATRICE EST SYMETRIQUE : LE =  UE
C
C            LES MATRICES "DE" SONT CONSIDEREES COMME DES DIAGONALES
C            DE TAILLE NPOIN * NPOIN QU'IL FAUT IMAGINER COMPLETEES
C            AVEC DES 1 POUR LES POINTS QUI N'APPARTIENNENT PAS A
C            L'ELEMENT CONSIDERE.
C
C            ON EFFECTUE ENSUITE LE PRODUIT DE TOUTES CES DIAGONALES
C            QUI DONNE LA DIAGONALE DB.
C
C  ATTENTION :
C
C  POUR LES MATRICES NON SYMETRIQUES UE <> LE
C  UE (LES BETAS) EST STOCKE DANS XB(.,1  A  6)
C  LE (LES ALFAS) EST STOCKE DANS XB(.,7  A 12)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      DB        |<-- |  DIAGONALE DE LA MATRICE RESULTAT
C |      XB        |<-- |  TERMES EXTRADIAGONAUX DE LA MATRICE RESULTAT
C |      TYPDIA    | -->|  TYPE DE DIAGONALE ( 'Q', 'I' , OU '0' )
C |      XA        | -->|  TERMES EXTRADIAGONAUX DE LA MATRICE A
C |      TYPEXA    | -->|  TYPE DE TERMES EXTRADIAGONAUX ('Q','S',OU'0')
C |      X,Y,Z     | -->|  COORDONNEES DU MAILLAGE.
C |      SURFAC    | -->|  SURFACE DES TRIANGLES.
C |      IKLE      | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      NPOIN     | -->|  DIMENSION DES TABLEAU
C |      W         |<-->|  TABLEAU CONTENANT DB NON ASSEMBLEE
C |      COPY      | -->|  SI .TRUE. A EST COPIEE DANS B
C |                |    |  SINON ON CONSIDERE QUE B EST DEJA REMPLIE
C |      LV        | -->|  LONGUEUR DU VECTEUR POUR LA VECTORISATION
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF, EX_DLDU21 => DLDU21
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: NELEM,NELMAX,LV,NPOIN
      DOUBLE PRECISION, INTENT(OUT) :: DB(NPOIN),XB(NELMAX,*)
      DOUBLE PRECISION, INTENT(IN)  :: XA(NELMAX,*)
      CHARACTER(LEN=1), INTENT(IN)  :: TYPDIA,TYPEXA
      INTEGER, INTENT(IN)           :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(OUT) :: W(NELMAX,4)
      LOGICAL, INTENT(IN)           :: COPY
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
      DOUBLE PRECISION Z(1),C
C
      DOUBLE PRECISION     A12,A13,A14
      DOUBLE PRECISION A21,A22,A23,A24
      DOUBLE PRECISION A31,A32,A33,A34
      DOUBLE PRECISION A41,A42,A43,A44
      DOUBLE PRECISION        BETA12,BETA13,BETA14
      DOUBLE PRECISION ALFA21,BETA22,BETA23,BETA24
      DOUBLE PRECISION ALFA31,ALFA32,BETA33,BETA34
      DOUBLE PRECISION ALFA41,ALFA42,ALFA43,BETA44
C
C-----------------------------------------------------------------------
C
C ON EXIGE UNE MATRICE A A DIAGONALE IDENTITE
C
      IF(TYPDIA(1:1).NE.'I'.AND.NCSIZE.LE.1) THEN
         IF (LNG.EQ.1) WRITE(LU,1000) TYPDIA(1:1)
         IF (LNG.EQ.2) WRITE(LU,1001) TYPDIA(1:1)
1000     FORMAT(1X,'DLDU21 (BIEF) : DIAGONALE DE A NON EGALE A I :',A1)
1001     FORMAT(1X,'DLDU21 (BIEF) : DIAGONAL OF A NOT EQUAL TO I :',A1)
         CALL PLANTE(0)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(TYPEXA(1:1).EQ.'S') THEN
C
        IF(COPY) CALL OV('X=Y     ' , XB , XA , Z , C  , NELMAX*6 )
C
        DO 100 IELEM=1,NELEM
C
C MATRICE A DECOMPOSER ( SYMETRIQUE AVEC DES 1 SUR LA DIAGONALE)
C
C LIGNE 1
C          A11 = 1.D0
           A12 = XA(IELEM,1)
           A13 = XA(IELEM,2)
           A14 = XA(IELEM,3)
C LIGNE 2
           A22 = 1.D0
           A23 = XA(IELEM,4)
           A24 = XA(IELEM,5)
C LIGNE 3
           A33 = 1.D0
           A34 = XA(IELEM,6)
C LIGNE 4
           A44 = 1.D0
C
C DECOMPOSITION L*U DE CROUT
C
           ALFA21 = A12
           ALFA31 = A13
           ALFA41 = A14
C
           BETA12 =  A12
           BETA22 =  A22 - ALFA21*BETA12
           ALFA32 = (A23 - ALFA31*BETA12)/BETA22
           ALFA42 = (A24 - ALFA41*BETA12)/BETA22
C
           BETA13 =  A13
           BETA23 =  A23 - ALFA21*BETA13
           BETA33 =  A33 - ALFA31*BETA13 - ALFA32*BETA23
           ALFA43 = (A34 - ALFA41*BETA13 - ALFA42*BETA23)/BETA33
C
           BETA14 =  A14
           BETA24 =  A24 - ALFA21*BETA14
           BETA34 =  A34 - ALFA31*BETA14 - ALFA32*BETA24
           BETA44 =  A44 - ALFA41*BETA14 - ALFA42*BETA24 - ALFA43*BETA34
C
C ON STOCKE DANS XB ET W2,...,W4
C
           XB(IELEM,1 ) = ALFA21
           XB(IELEM,2 ) = ALFA31
           XB(IELEM,3 ) = ALFA41
           XB(IELEM,4 ) = ALFA32
           XB(IELEM,5 ) = ALFA42
           XB(IELEM,6 ) = ALFA43
C
           W(IELEM,2)    = BETA22
           W(IELEM,3)    = BETA33
           W(IELEM,4)    = BETA44
C
100     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(TYPEXA(1:1).EQ.'Q') THEN
C
        IF(COPY) CALL OV('X=Y     ' , XB , XA , Z , C  , NELMAX*12 )
C
        DO 200 IELEM=1,NELEM
C
C MATRICE A DECOMPOSER ( AVEC DES 1 SUR LA DIAGONALE)
C
C          A11 = 1.D0
           A22 = 1.D0
           A33 = 1.D0
           A44 = 1.D0
C
           A12 = XA(IELEM,1 )
           A13 = XA(IELEM,2 )
           A14 = XA(IELEM,3 )
           A23 = XA(IELEM,4 )
           A24 = XA(IELEM,5 )
           A34 = XA(IELEM,6 )
C
           A21 = XA(IELEM,7 )
           A31 = XA(IELEM,8 )
           A41 = XA(IELEM,9 )
           A32 = XA(IELEM,10)
           A42 = XA(IELEM,11)
           A43 = XA(IELEM,12)
C
C DECOMPOSITION L*U DE CROUT
C
           ALFA21 = A21
           ALFA31 = A31
           ALFA41 = A41
C
           BETA12 =  A12
           BETA22 =  A22 - ALFA21*BETA12
           ALFA32 = (A32 - ALFA31*BETA12)/BETA22
           ALFA42 = (A42 - ALFA41*BETA12)/BETA22
C
           BETA13 =  A13
           BETA23 =  A23 - ALFA21*BETA13
           BETA33 =  A33 - ALFA31*BETA13 - ALFA32*BETA23
           ALFA43 = (A43 - ALFA41*BETA13 - ALFA42*BETA23)/BETA33
C
           BETA14 =  A14
           BETA24 =  A24 - ALFA21*BETA14
           BETA34 =  A34 - ALFA31*BETA14 - ALFA32*BETA24
           BETA44 =  A44 - ALFA41*BETA14 - ALFA42*BETA24 - ALFA43*BETA34
C
C ON STOCKE DANS XB ET W2,...,W4
C ET ON FAIT EN MEME TEMPS LA DECOMPOSITION L D U
C
           XB(IELEM,1 ) = BETA12
           XB(IELEM,2 ) = BETA13
           XB(IELEM,3 ) = BETA14
           XB(IELEM,4 ) = BETA23/BETA22
           XB(IELEM,5 ) = BETA24/BETA22
           XB(IELEM,6 ) = BETA34/BETA33
C
           XB(IELEM,07) = ALFA21
           XB(IELEM,08) = ALFA31
           XB(IELEM,09) = ALFA41
           XB(IELEM,10) = ALFA32
           XB(IELEM,11) = ALFA42
           XB(IELEM,12) = ALFA43
C
           W(IELEM,2)    = BETA22
           W(IELEM,3)    = BETA33
           W(IELEM,4)    = BETA44
C
200     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
         IF (LNG.EQ.1) WRITE(LU,2000) TYPEXA(1:1)
         IF (LNG.EQ.2) WRITE(LU,2001) TYPEXA(1:1)
2000     FORMAT(1X,'DLDU21 (BIEF) : TYPE DE MATRICE NON PREVU :',A1)
2001     FORMAT(1X,'DLDU21 (BIEF) : TYPE OF MATRIX NOT AVAILABLE :',A1)
         CALL PLANTE(0)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  ASSEMBLAGE MULTIPLICATIF DE LA DIAGONALE AVEC INITIALISATION
C  DE DB A 1. ON SAUTE IKLE1 CAR W1 = 1.
C
      CALL ASMVEC(DB,IKLE(1,2),NPOIN,NELEM,NELMAX,3,W(1,2),.TRUE.,LV)
C
C  INVERSION DE DB
C
      CALL OV( 'X=1/Y   ' , DB , DB , Z , C , NPOIN )
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
