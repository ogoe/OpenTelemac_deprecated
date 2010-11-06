C                       *****************
                        SUBROUTINE SLOP10
C                       *****************
C
     *(COEF,XEL,YEL,Z,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C  BIEF VERSION 5.1          27/01/95    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C FONCTION : CALCUL DU COEFFICIENT 1 / COS(ALFA)
C
C            OU ALFA EST LA PENTE D'UN ELEMENT TRIANGULAIRE.
C
C            CE COEFFICIENT EST UTILISE POUR LES TERMES DE FROTTEMENT
C            SUR LE FOND.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      COEF      |<-- |  RESULTAT
C |      XEL,YEL   | -->|  COORDONNEES DES POINTS DES ELEMENTS
C |      Z         | -->|  COTES DU FOND (DONNEES PAR POINT)
C |      IKLE      | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.              |
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.      |
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)               |
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : PROPAG
C
C  SOUS-PROGRAMME APPELE : NEANT
C
C-----------------------------------------------------------------------
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C 
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*)
C
      DOUBLE PRECISION, INTENT(INOUT) :: COEF(NELEM)
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,*),YEL(NELMAX,*),Z(*)     
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION X2,X3,Y2,Y3,Z2,Z3,A,B,C
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
         DO 1 IELEM = 1 , NELEM
C
         X2 = XEL(IELEM,2)
         X3 = XEL(IELEM,3)
C
         Y2 = YEL(IELEM,2)
         Y3 = YEL(IELEM,3)
C
         Z2 = Z(IKLE(IELEM,2)) -  Z(IKLE(IELEM,1))
         Z3 = Z(IKLE(IELEM,3)) -  Z(IKLE(IELEM,1))
C
         A = (X2*Y3-X3*Y2)**2
         B = (Y2*Z3-Z2*Y3)**2
         C = (X3*Z2-Z3*X2)**2
C
         COEF(IELEM) = SQRT( (A+B+C)/A )
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
