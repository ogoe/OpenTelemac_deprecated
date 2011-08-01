C                       *****************
                        SUBROUTINE MT13AA
C                       *****************
C
     *(  A11 , A12 , A13 ,
     *   A21 , A22 , A23 ,
     *   A31 , A32 , A33 ,
     *   XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        C  MOULIN    (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  |||||||||| ATTENTION : SIGNE CHANGE PAR RAPPORT A 3.0
C  |||||||||| ATTENTION : TRANSPOSITION PAR RAPPORT A 3.0
C
C                    /            D
C    A(I,J)=   XMUL /  PSI2(I) *  --( PSI1(J) ) D(OMEGA)
C                  /OMEGA         DX
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C  AVEC ICOORD=3 ON AURAIT UNE DERIVEE SUIVANT Z
C
C  PSI1 : LINEAIRE
C  PSI2 : LINEAIRE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SF,SG,SH   | -->|  STRUCTURES DE F,G ET H.
C |     SU,SV,SW   | -->|  STRUCTURES DE U,V ET W.
C |     F,G,H      | -->|  FONCTIONS INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     U,V,W      | -->|  COMPOSANTES D'UN VECTEUR INTERVENANT DANS LE
C |                |    |  CALCUL DE LA MATRICE.
C |     XEL,YEL,ZEL| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     IKLE1..3   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |     ICOORD     | -->|  1: DERIVEE SUIVANT X, 2:SUIVANT Y
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : ASSVEC , OV
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT13AA => MT13AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION XSUR6
C
C-----------------------------------------------------------------------
C
      XSUR6 = XMUL/6.D0
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
      IF(ICOORD.EQ.1) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM = 1 , NELEM
C
C   TERMES DIAGONAUX
C
      A22(IELEM) =    YEL(IELEM,3) * XSUR6
      A33(IELEM) = -  YEL(IELEM,2) * XSUR6
      A11(IELEM) = - A22(IELEM) - A33(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A21(IELEM) = A11(IELEM)
      A31(IELEM) = A11(IELEM)
      A32(IELEM) = A22(IELEM)
      A12(IELEM) = A22(IELEM)
      A13(IELEM) = A33(IELEM)
      A23(IELEM) = A33(IELEM)
C
1     CONTINUE
C
      ELSEIF(ICOORD.EQ.2) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
      DO 2 IELEM = 1 , NELEM
C
C   TERMES DIAGONAUX
C
      A22(IELEM) = - XEL(IELEM,3) * XSUR6
      A33(IELEM) =   XEL(IELEM,2) * XSUR6
      A11(IELEM) = - A22(IELEM) - A33(IELEM)
C
C   TERMES EXTRADIAGONAUX
C
      A21(IELEM) = A11(IELEM)
      A31(IELEM) = A11(IELEM)
      A32(IELEM) = A22(IELEM)
      A12(IELEM) = A22(IELEM)
      A13(IELEM) = A33(IELEM)
      A23(IELEM) = A33(IELEM)
C
2     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'MT13AA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT13AA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(0)
          STOP
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
