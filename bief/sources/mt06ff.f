C                       *****************
                        SUBROUTINE MT06FF
C                       *****************
C
     *( A11 , A12 , A13 , A14 ,
     *        A22 , A23 , A24 ,
     *              A33 , A34 ,
     *                    A44 ,
     *  XMUL,SF,F,X,Y,Z,IKLE1,IKLE2,IKLE3,IKLE4,NBOR,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.3           18/08/94    J-M HERVOUET (LNH) 30 87 80 18
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
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     L'ELEMENT EST LE QUADRILATERE Q1,  DANS UN MAILLAGE DE PRISMES
C
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      USE BIEF, EX_MT06FF => MT06FF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NBOR(*),NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*),A14(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*),A34(*)
      DOUBLE PRECISION, INTENT(INOUT) ::                      A44(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTRINSIC SQRT
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      INTEGER IELMF,I1,I2,I3,I4,IELEM
C
      DOUBLE PRECISION SUR720,AL,S1,S2,S11112,S11122,S11222,S12222
      DOUBLE PRECISION F14,F23,F1114,F2223,F2333,F1444
C
C**********************************************************************
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C
C     F LINEAIRE PAR FACE DE BORD
C
      IF(IELMF.EQ.71) THEN
C
         SUR720  = XMUL/720.D0
C
C   BOUCLE SUR LES FACES DE BORD
C
         DO 1 IELEM = 1,NELEM
C
C  NUMEROTATION GLOBALE DES SOMMETS DE LA FACE
C
            I1 = IKLE1(IELEM)
            I2 = IKLE2(IELEM)
            I3 = IKLE3(IELEM)
            I4 = IKLE4(IELEM)
C
            AL = SQRT((X(NBOR(I2))-X(NBOR(I1)))**2
     *               +(Y(NBOR(I2))-Y(NBOR(I1)))**2) * SUR720
C
            S1 = (Z(NBOR(I4)) - Z(NBOR(I1))) * AL
            S2 = (Z(NBOR(I3)) - Z(NBOR(I2))) * AL
            S11112 = S1 + S1 + S1 + S1 + S2
            S11122 = S1 + S1 + S1 + S2 + S2
            S11222 = S1 + S1 + S2 + S2 + S2
            S12222 = S1 + S2 + S2 + S2 + S2
C
            F14 = F(I1) + F(I4)
            F23 = F(I2) + F(I3)
            F1114 = F(I1) + F(I1) + F14
            F2223 = F(I2) + F(I2) + F23
            F2333 = F23 + F(I3) + F(I3)
            F1444 = F14 + F(I4) + F(I4)
C
C  TERMES DIAGONAUX
C
            A11(IELEM) = 3*F1114*S11112 + F2223*S11122
            A22(IELEM) = 3*F2223*S12222 + F1114*S11222
            A33(IELEM) = 3*F2333*S12222 + F1444*S11222
            A44(IELEM) = 3*F1444*S11112 + F2333*S11122
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
            A12(IELEM) = F1114*S11122 + F2223*S11222
            A13(IELEM) =   F14*S11122 +   F23*S11222
            A14(IELEM) = 3*F14*S11112 +   F23*S11122
            A23(IELEM) = 3*F23*S12222 +   F14*S11222
            A24(IELEM) = A13(IELEM)
            A34(IELEM) = F2333*S11222 + F1444*S11122
C
1        CONTINUE
C
C-----------------------------------------------------------------------
C
C     AUTRES TYPES DE DISCRETISATION DE F
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
         IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100      FORMAT(1X,'MT06FF (BIEF) :',/,
     *          1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *          1X,'NOM REEL : ',A6)
101      FORMAT(1X,'MT06FF (BIEF) :',/,
     *          1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *          1X,'REAL NAME: ',A6)
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
