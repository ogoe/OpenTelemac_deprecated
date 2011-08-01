C                       *****************
                        SUBROUTINE MT06AA
C                       *****************
C
     *( A11 , A12 , A13 ,
     *        A22 , A23 ,
     *              A33 ,
     *  XMUL,SF,F,SURFAC,IKLE1,IKLE2,IKLE3,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           29/10/99    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C      FONCTION:
C      ========:
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                              /
C                    A    =   /  F * (P *P )*J(X,Y) DXDY
C                     I J    /S        I  J
C
C     PAR MAILLE ELEMENTAIRE .
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
C |     SF,        | -->|  STRUCTURE DE F
C |     F          | -->|  FONCTION INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     SURFAC     | -->|  SURFACE DES TRIANGLES.
C |     IKLE1..3   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C     -------------
C     | ATTENTION | : LE JACOBIEN DOIT ETRE POSITIF .
C     -------------
C**********************************************************************
C
      USE BIEF, EX_MT06AA => MT06AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*) 
C
C     STRUCTURE DE F 
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELMF,IELEM
C
      DOUBLE PRECISION SUR12,SUR60,DET1,DET2,F123,F1,F2,F3
C
C-----------------------------------------------------------------------
C
      SUR60 = XMUL/60.D0
      SUR12 = XMUL/12.D0
C
C-----------------------------------------------------------------------
C
      IELMF = SF%ELM
C
C-----------------------------------------------------------------------
C
C  CAS OU F EST LINEAIRE
C
      IF(IELMF.EQ.11) THEN
C
      DO 1 IELEM = 1 , NELEM
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F3 = F(IKLE3(IELEM))
C
      F123 = F1 + F2 + F3
C
      DET1 = SURFAC(IELEM) * SUR60
      DET2 = DET1 + DET1
C
C***********************************************************************
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = DET1 * (F123+F123-F3)
      A13(IELEM) = DET1 * (F123+F123-F2)
      A23(IELEM) = DET1 * (F123+F123-F1)
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = DET2 * (F123+F1+F1)
      A22(IELEM) = DET2 * (F123+F2+F2)
      A33(IELEM) = DET2 * (F123+F3+F3)
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
C  CAS OU F EST LINEAIRE
C
      ELSEIF(IELMF.EQ.10) THEN
C
      DO 2 IELEM = 1 , NELEM
C
      F1 = F(IELEM)
C
      DET1 = SURFAC(IELEM) * SUR12
      DET2 = DET1 + DET1
C
C***********************************************************************
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = DET1 * F1
      A13(IELEM) = DET1 * F1
      A23(IELEM) = DET1 * F1
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = DET2 * F1
      A22(IELEM) = DET2 * F1
      A33(IELEM) = DET2 * F1
C
2     CONTINUE
C
C     AUTRES TYPES DE DISCRETISATION DE F
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'MT06AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT06AA (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *        1X,'REAL NAME: ',A6)
       CALL PLANTE(1)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
