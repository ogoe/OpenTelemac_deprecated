C                       *****************
                        SUBROUTINE MT01AA
C                       *****************
C
     *( A11 , A12 , A13 ,
     *        A22 , A23 ,
     *              A33 ,
     *  XMUL,SURFAC,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           28/11/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE DE MASSE POUR DES TRIANGLES P1
C
C-----------------------------------------------------------------------
C
C      FONCTION:
C      ========:
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                                 /
C                    A    = XMUL /  (P *P )*J(X,Y) DXDY
C                     I J       /S    I  J
C
C    PAR MAILLE ELEMENTAIRE.
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     L'ELEMENT EST LE TRIANGLE P1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : AUCUN
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT01AA => MT01AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM
      DOUBLE PRECISION SUR12,DET
C
C=======================================================================
C
      SUR12 = XMUL/12.D0
C
C-----------------------------------------------------------------------
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM = 1 , NELEM
C
      DET = SURFAC(IELEM) * SUR12
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = DET
      A13(IELEM) = DET
      A23(IELEM) = DET
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = DET + DET
      A22(IELEM) = DET + DET
      A33(IELEM) = DET + DET
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
