C                       *****************
                        SUBROUTINE DLDU11
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
C FONCTION : POUR LES TRIANGLES P1
C
C            DECOMPOSITION L D U DES MATRICES ELEMENTAIRES CONTENUES
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
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      DB        |<-- |  DIAGONALE DE LA MATRICE RESULTAT
C |      XB        |<-- |  TERMES EXTRADIAGONAUX DE LA MATRICE RESULTAT
C |      TYPDIA    |<-- |  TYPE DE DIAGONALE ( 'Q', 'I' , OU '0' )
C |      XA        |<-- |  TERMES EXTRADIAGONAUX DE LA MATRICE A
C |      TYPEXA    |<-- |  TYPE DE TERMES EXTRADIAGONAUX ('Q','S',OU'0')
C |      X,Y,Z     | -->|  COORDONNEES DU MAILLAGE.
C |      SURFAC    | -->|  SURFACE DES TRIANGLES.
C |      IKLE      | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      NPOIN     | -->|  DIMENSION DES TABLEAUX
C |      W         | -->|  TABLEAU DE TRAVAIL DE DIMENSION (NELMAX,3)
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
      USE BIEF, EX_DLDU11 => DLDU11
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
      DOUBLE PRECISION, INTENT(OUT) :: W(NELMAX,3)
      LOGICAL, INTENT(IN)           :: COPY
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
C
      DOUBLE PRECISION Z(1),C
C
C-----------------------------------------------------------------------
C
C ON EXIGE UNE MATRICE A A DIAGONALE IDENTITE (SAUF EN PARALLELISME)
C
      IF(TYPDIA(1:1).NE.'I'.AND.NCSIZE.LE.1) THEN
         IF (LNG.EQ.1) WRITE(LU,100) TYPDIA(1:1)
         IF (LNG.EQ.2) WRITE(LU,101) TYPDIA(1:1)
100      FORMAT(1X,'DLDU11 (BIEF) : DIAGONALE DE A NON EGALE A I :',A1)
101      FORMAT(1X,'DLDU11 (BIEF) : DIAGONAL OF A NOT EQUAL TO I :',A1)
         CALL PLANTE(0)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(TYPEXA(1:1).EQ.'S') THEN
C
        IF(COPY) CALL OV('X=Y     ' , XB , XA , Z , C  , NELMAX*3 )
C
        DO 10 IELEM = 1 , NELEM
         W(IELEM,2) = 1.D0 - XB(IELEM,1)**2
         XB(IELEM,3) = (XB(IELEM,3)-XB(IELEM,1)*XB(IELEM,2))/W(IELEM,2)
         W(IELEM,3) = 1.D0 - XB(IELEM,2)**2 -XB(IELEM,3)**2
10      CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(TYPEXA(1:1).EQ.'Q') THEN
C
        IF(COPY) CALL OV('X=Y     ' , XB , XA , Z , C  , NELMAX*6 )
C
        DO 20 IELEM = 1 , NELEM
C DECOMPOSITION L U
         W(IELEM,2)=1.D0 - XB(IELEM,1)*XB(IELEM,4)
         XB(IELEM,6) = (XB(IELEM,6)-XB(IELEM,1)*XB(IELEM,5))/W(IELEM,2)
         XB(IELEM,3) =  XB(IELEM,3)-XB(IELEM,4)*XB(IELEM,2)
         W(IELEM,3)=1.D0-XB(IELEM,2)*XB(IELEM,5)-XB(IELEM,3)*XB(IELEM,6)
C PASSAGE A LA DECOMPOSITION L D U
         XB(IELEM,3) = XB(IELEM,3) / W(IELEM,2)
20      CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSE
         IF (LNG.EQ.1) WRITE(LU,200) TYPEXA(1:1)
         IF (LNG.EQ.2) WRITE(LU,201) TYPEXA(1:1)
200      FORMAT(1X,'DLDU11 (BIEF) : TYPE DE MATRICE NON PREVU :',A1)
201      FORMAT(1X,'DLDU11 (BIEF) : TYPE OF MATRIX NOT AVAILABLE :',A1)
         CALL PLANTE(0)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  ASSEMBLAGE MULTIPLICATIF DE LA DIAGONALE AVEC INITIALISATION
C  DE DB A 1. ON SAUTE IKLE1 CAR W1 = 1.
C
      CALL ASMVEC(DB,IKLE(1,2),NPOIN,NELEM,NELMAX,2,W(1,2),.TRUE.,LV)
C
C  INVERSION DE DB
C
      CALL OV( 'X=1/Y   ' , DB , DB , Z , C , NPOIN )
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
 
