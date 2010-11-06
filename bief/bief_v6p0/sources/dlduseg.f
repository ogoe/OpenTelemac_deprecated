C                       ******************
                        SUBROUTINE DLDUSEG
C                       ******************
C
     *(DB,XB,TYPDIA,XA,TYPEXA,GLOSEG,NSEG,NPOIN,COPY)
C
C***********************************************************************
C BIEF VERSION 5.5           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : POUR LES SEGMENTS
C
C            DECOMPOSITION L D U DES MATRICES ELEMENTAIRES PAR SEGMENT
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
C
C      (  1   X12 )   (  1   0 ) (  1       0     ) (  1   X12 )
C      (          ) = (        ) (                ) (          )
C      ( X21   1  )   ( X21  1 ) (  0   1-X12*X21 ) (  0     1 )
C
C
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
C |      GLOSEG    | -->|  PASSAGE DE LA NUMEROTATION SEGMENT A GLOBALE
C |      NSEG      | -->|  NOMBRE DE SEGMENTS
C |      NPOIN     | -->|  DIMENSION DES TABLEAUX
C |      COPY      | -->|  SI .TRUE. A EST COPIEE DANS B
C |                |    |  SINON ON CONSIDERE QUE B EST DEJA REMPLIE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES :
C
C**********************************************************************
C
      USE BIEF, EX_DLDUSEG => DLDUSEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: NSEG,NPOIN
      DOUBLE PRECISION, INTENT(OUT) :: DB(NPOIN),XB(NSEG,*)
      DOUBLE PRECISION, INTENT(IN)  :: XA(NSEG,*)
      CHARACTER(LEN=1), INTENT(IN)  :: TYPDIA,TYPEXA
      INTEGER, INTENT(IN)           :: GLOSEG(NSEG,2)
      LOGICAL, INTENT(IN)           :: COPY
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISEG
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
100      FORMAT(1X,'DLDUSEG (BIEF) : DIAGONALE DE A NON IDENTITE :',A1)
101      FORMAT(1X,'DLDUSEG (BIEF) : DIAGONAL OF A NOT IDENTITY :',A1)
         CALL PLANTE(1)
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      IF(TYPEXA(1:1).EQ.'S') THEN
C
        IF(COPY) THEN
          CALL OV('X=Y     ' , XB , XA , Z , C , NSEG )
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(TYPEXA(1:1).EQ.'Q') THEN
C
        IF(COPY) THEN
          CALL OV('X=Y     ' , XB , XA , Z , C , 2*NSEG )
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,200) TYPEXA(1:1)
        IF (LNG.EQ.2) WRITE(LU,201) TYPEXA(1:1)
200     FORMAT(1X,'DLDUSEG (BIEF) : TYPE DE MATRICE NON PREVU :',A1)
201     FORMAT(1X,'DLDUSEG (BIEF) : TYPE OF MATRIX NOT TREATED:',A1)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
C  ASSEMBLAGE MULTIPLICATIF DE LA DIAGONALE AVEC INITIALISATION
C  DE DB A 1.
C
      CALL OV('X=C     ' , DB , DB , DB , 1.D0 , NPOIN )
C
      IF(TYPEXA(1:1).EQ.'S') THEN
C
      DO ISEG=1,NSEG
        DB(GLOSEG(ISEG,2))=DB(GLOSEG(ISEG,2))*(1.D0-XB(ISEG,1)**2)
      ENDDO
C
      ELSE
C
      DO ISEG=1,NSEG
        DB(GLOSEG(ISEG,2))=
     *  DB(GLOSEG(ISEG,2))*(1.D0-XB(ISEG,1)*XB(ISEG,2))
      ENDDO
C
      ENDIF
C
C  INVERSION DE DB (RISQUE DE DIVISION PAR ZERO)
C
      CALL OV( 'X=1/Y   ' , DB , DB , Z , C , NPOIN )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
