C                       *****************
                        SUBROUTINE MT06OO
C                       *****************
C
     *(A11,A12,A22,XMUL,SF,F,LGSEG,IKLE1,IKLE2,NBOR,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.9        20/03/08    J-M HERVOUET (LNHE) 01 30 87 80 18
C                         04/09/92    F  LEPEINTRE (LNH) 30 87 78 54
C
C***********************************************************************
C
C FONCTION :
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                              /
C                    A    =   /  F (P *P )*J(X,Y) DX
C                     I J    /L      I  J
C
C     PAR MAILLE ELEMENTAIRE .
C
C     J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C     L'ELEMENT EST LE SEGMENT P1
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
      USE BIEF, EX_MT06OO => MT06OO
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(*),IKLE2(*),NBOR(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
C
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN)    :: LGSEG(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: A11(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: A12(NELMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: A22(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION SUR12,DET1,F1,F2,F12
C
C-----------------------------------------------------------------------
C
      SUR12  = XMUL/12.D0
C
C-----------------------------------------------------------------------
C
      IELMF = SF%ELM
C
C     F CONSTANTE PAR SEGMENT, DANS UN TABLEAU DE BORD
C
      IF(IELMF.EQ.0) THEN
C
      DO 1 IELEM = 1 , NELEM
C  A OPTIMISER
      F1 = F(IELEM)
      F2 = F(IELEM)
C
      F12 = F1 + F2
C
      DET1 = LGSEG(IELEM) * SUR12
C
      A11(IELEM) = DET1 * (F12+2*F1)
      A12(IELEM) = DET1 * F12
      A22(IELEM) = DET1 * (F12+2*F2)
C
1     CONTINUE
C
C     F LINEAIRE PAR SEGMENT, DANS UN TABLEAU DE BORD
C     NOTE : IKLE EST ICI UN IKLE DE BORD
C
      ELSEIF(IELMF.EQ.1) THEN
C
      DO 2 IELEM = 1 , NELEM
C
      F1 = F(IKLE1(IELEM))
      F2 = F(IKLE2(IELEM))
      F12 = F1 + F2
      DET1 = LGSEG(IELEM) * SUR12
      A11(IELEM) = DET1 * (F12+2*F1)
      A12(IELEM) = DET1 * F12
      A22(IELEM) = DET1 * (F12+2*F2)
C
2     CONTINUE
C
C     F LINEAIRE, DANS UN TABLEAU DEFINI SUR LE DOMAINE
C
      ELSEIF(IELMF.EQ.11.OR.IELMF.EQ.21) THEN
C
      DO 3 IELEM = 1 , NELEM
C
      F1 = F(NBOR(IKLE1(IELEM)))
      F2 = F(NBOR(IKLE2(IELEM)))
      F12 = F1 + F2
      DET1 = LGSEG(IELEM) * SUR12
      A11(IELEM) = DET1 * (F12+2*F1)
      A12(IELEM) = DET1 * F12
      A22(IELEM) = DET1 * (F12+2*F2)
C
3     CONTINUE
C
C     AUTRES TYPES DE DISCRETISATION DE F
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'MT0600 (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT0600 (BIEF) :',/,
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
