C                       *****************
                        SUBROUTINE DESSEG
C                       *****************
C
     *(X, XA,TYPEXA,B,GLOSEG,NSEG,NPOIN,DITR,COPY)
C
C***********************************************************************
C BIEF VERSION 5.5           25/02/04    J-M HERVOUET (LNH) 30 87 80 18
C                                        
C***********************************************************************
C
C FONCTION : RESOLUTION DU SYSTEME L X = B (SEGMENT PAR SEGMENT)
C
C            ICI LA MATRICE L EST LE RESULTAT D'UNE DECOMPOSITION
C            EFFECTUEE PAR LE SOUS-PROGRAMME DECLDU.
C
C            CHAQUE MATRICE ELEMENTAIRE A ETE DECOMPOSEE SOUS LA FORME :
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
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  SOLUTION DU SYSTEME AX = B
C |      XA        |<-- |  TERMES EXTRADIAGONAUX DE LA MATRICE A
C |      B         |<-- |  SECOND MEMBRE DU SYSTEME A RESOUDRE.
C |      GLOSEG    | -->|  PASSAGE DE LA NUMEROTATION SEGMENT A GLOBALE
C |      NSEG      | -->|  NOMBRE DE SEGMENTS
C |      NPOIN     | -->|  DIMENSION DES TABLEAUX
C |      DITR      | -->|  CARACTERE  'D' : ON CALCULE AVEC A
C |                |    |             'T' : ON CALCULE AVEC A TRANSPOSEE
C |      COPY      | -->|  SI .TRUE. B EST RECOPIE DANS X AU PREALABLE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : OV , PLANTE
C
C**********************************************************************
C
      USE BIEF, EX_DESSEG => DESSEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER         , INTENT(IN)    :: NPOIN,NSEG
      INTEGER         , INTENT(IN)    :: GLOSEG(NSEG,2)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XA(NSEG,*),B(NPOIN)
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXA,DITR
      LOGICAL         , INTENT(IN)    :: COPY
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
C 1) INITIALISATION : X = SECOND MEMBRE
C
      IF(COPY) CALL OV( 'X=Y     ' , X , B , B , 0.D0 , NPOIN )
C
C-----------------------------------------------------------------------
C
C 2) INVERSION DES MATRICES TRIANGULAIRES INFERIEURES (DESCENTE)
C
      IF(TYPEXA(1:1).EQ.'S' .OR.
     *  (TYPEXA(1:1).EQ.'Q'.AND.DITR(1:1).EQ.'T')) THEN
C
        DO I=1,NSEG
          X(GLOSEG(I,2))=X(GLOSEG(I,2))-XA(I,1)*X(GLOSEG(I,1))
        ENDDO
C
      ELSEIF(TYPEXA(1:1).EQ.'Q'.AND.DITR(1:1).EQ.'D') THEN
C
        DO I=1,NSEG
          X(GLOSEG(I,2))=X(GLOSEG(I,2))-XA(I,2)*X(GLOSEG(I,1))
        ENDDO
C
      ELSE
        WRITE(LU,*) 'DESSEG, CASE NOT IMPLEMENTED'
        WRITE(LU,*) '        TYPEXA=',TYPEXA,' DITR=',DITR(1:1)
        CALL PLANTE(1)
        STOP 
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
