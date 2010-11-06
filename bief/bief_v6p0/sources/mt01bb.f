C                       *****************
                        SUBROUTINE MT01BB
C                       *****************
C
     *( A11 , A12 , A13 , A14 ,
     *        A22 , A23 , A24 ,
     *              A33 , A34 ,
     *                    A44 ,
     *  XMUL,SURFAC,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95     J-M HERVOUET (LNH) 30 87 80 18
C                                           C   MOULIN (LNH) 30 87 83 81
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE DE MASSE POUR LES
C            TRIANGLES QUASI-BULLE
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
      USE BIEF, EX_MT01BB => MT01BB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*),A14(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*),A34(*)
      DOUBLE PRECISION, INTENT(INOUT) ::                      A44(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELEM
      DOUBLE PRECISION XSUR36
C
C=======================================================================
C
      XSUR36 = XMUL / 36.D0
C
      DO 5 IELEM = 1 , NELEM
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) =     SURFAC(IELEM)*XSUR36
         A13(IELEM) =     A12(IELEM)
         A23(IELEM) =     A12(IELEM)
         A14(IELEM) = 2 * A12(IELEM)
         A24(IELEM) =     A14(IELEM)
         A34(IELEM) =     A14(IELEM)
C
C  TERMES DIAGONAUX
C
         A11(IELEM) = 2 * A14(IELEM)
         A22(IELEM) =     A11(IELEM)
         A33(IELEM) =     A11(IELEM)
         A44(IELEM) = 3 * A14(IELEM)
C
5     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
