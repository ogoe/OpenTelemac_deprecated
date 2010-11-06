C                       *****************
                        SUBROUTINE MT01CC
C                       *****************
C
     *( A11 , A12 , A13 , A14 , A15 , A16 ,
     *        A22 , A23 , A24 , A25 , A26 ,
     *              A33 , A34 , A35 , A36 ,
     *                    A44 , A45 , A46 ,
     *                          A55 , A56 ,
     *                                A66 ,
     *  XMUL,SURFAC,NELEM,NELMAX )
C
C***********************************************************************
C BIEF VERSION 5.9   29/02/08   ALGIANE FROEHLY (MATMECA) 01 30 87 80 18
C                                           
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE DE MASSE POUR LES
C            TRIANGLES P2
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
      USE BIEF !, EX_MT01CC => MT01CC
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
      DOUBLE PRECISION, INTENT(INOUT) :: A15(*),A16(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A25(*),A26(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A33(*),A34(*),A35(*),A36(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A44(*),A45(*),A46(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A55(*),A56(*),A66(*)                   
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
      DOUBLE PRECISION XSUR180
C
C=======================================================================
C
      XSUR180 = XMUL / 180.D0
C
      DO IELEM = 1 , NELEM
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = - SURFAC(IELEM) * XSUR180
         A13(IELEM) =          A12(IELEM)
         A14(IELEM) =   0.D0
         A15(IELEM) =   4.D0 * A12(IELEM)
         A16(IELEM) =   0.D0
         A23(IELEM) =          A12(IELEM)
         A24(IELEM) =   0.D0
         A25(IELEM) =   0.D0
         A26(IELEM) =          A15(IELEM)
         A34(IELEM) =          A15(IELEM)
         A35(IELEM) =   0.D0
         A36(IELEM) =   0.D0
         A45(IELEM) = - 4.D0 * A15(IELEM)
         A46(IELEM) =          A45(IELEM)
         A56(IELEM) =          A45(IELEM)
C
C  TERMES DIAGONAUX
C
         A11(IELEM) = - 6.D0 * A12(IELEM)
         A22(IELEM) =          A11(IELEM)
         A33(IELEM) =          A11(IELEM)
         A44(IELEM) = - 8.D0 * A15(IELEM)
         A55(IELEM) =          A44(IELEM)
         A66(IELEM) =          A44(IELEM)
C
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
