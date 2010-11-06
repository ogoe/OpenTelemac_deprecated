C                       *****************
                        SUBROUTINE MT11AA
C                       *****************
C
     *(  A11 , A12 , A13 ,
     *   A21 , A22 , A23 ,
     *   A31 , A32 , A33 ,
     *   XMUL,SF,F,XEL,YEL,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        C  MOULIN    (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  EXEMPLE AVEC ICOORD=1
C
C
C                 /           D
C A(I,J)=- XMUL  /  PSI2(J) * -- ( PSI1(I) F ) ) D(OMEGA)
C               /OMEGA        DX
C
C
C  ATTENTION AU SIGNE MOINS ||
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : LINEAIRES
C  PSI2 : LINEAIRES
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
      USE BIEF, EX_MT11AA => MT11AA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX),IKLE3(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF
C
      DOUBLE PRECISION SUR24,X2,X3,Y2,Y3,F1,F2,F3
C
C-----------------------------------------------------------------------
C
      SUR24 = XMUL/24.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C  MEME RESULTAT SI F EST LINEAIRE OU QUASI-BULLE
C
      IF(IELMF.EQ.11.OR.IELMF.EQ.12) THEN
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
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      Y2 = YEL(IELEM,2)
      Y3 = YEL(IELEM,3)
C
      F1  =  F(IKLE1(IELEM)) * SUR24
      F2  =  F(IKLE2(IELEM)) * SUR24
      F3  =  F(IKLE3(IELEM)) * SUR24
C
C   TERMES DIAGONAUX
C
      A11(IELEM) =     Y2  * (F3-F2-4*F1)  +     Y3  * ( F3-F2+4*F1)
      A22(IELEM) = (Y2+Y2) * (F3-F1)       +     Y3  * (-F3-4*F2+F1)
      A33(IELEM) =     Y2  * (4*F3+F2-F1)  + (Y3+Y3) * (-F2+F1)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM)  =-(Y2+Y2) * (F2+F1)        +     Y3  * (F3+F2+F1+F1)
      A13(IELEM)  =      Y2 * (-F3-F2-F1-F1) + (Y3+Y3) * (F3+F1)
      A23(IELEM)  =      Y2 * (F3-F1)        - (Y3+Y3) * (F3+F2)
      A21(IELEM)  =      Y2 * (F3-F1)        +      Y3 * (-F3-F2-F2-F1)
      A31(IELEM)  =      Y2 * (F3+F3+F2+F1)  +      Y3 * (-F2+F1)
      A32(IELEM)  = (Y2+Y2) * (F3+F2)        +      Y3 * (-F2+F1)
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
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      X2  =  XEL(IELEM,2)
      X3  =  XEL(IELEM,3)
C
      F1  =  F(IKLE1(IELEM)) * SUR24
      F2  =  F(IKLE2(IELEM)) * SUR24
      F3  =  F(IKLE3(IELEM)) * SUR24
C
C   TERMES DIAGONAUX
C
      A11(IELEM) =     X2  * (-F3+F2+4*F1) +   X3  * (-F3+F2-4*F1)
      A22(IELEM) = (X2+X2) * (-F3+F1)      +   X3  * (F3+4*F2-F1)
      A33(IELEM) =     X2  * (-4*F3-F2+F1) + (X3+X3) * (F2-F1)
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM)  = (X2+X2) * (F2+F1)        +     X3  * (-F3-F2-F1-F1)
      A13(IELEM)  =      X2 * (F3+F2+F1+F1)  - (X3+X3) * (F3+F1)
      A23(IELEM)  =      X2 * (-F3+F1)       + (X3+X3) * (F3+F2)
      A21(IELEM)  =      X2 * (-F3+F1)       +      X3 * (F3+F2+F2+F1)
      A31(IELEM)  =      X2 * (-F3-F3-F2-F1) +      X3 * (F2-F1)
      A32(IELEM)  =-(X2+X2) * (F3+F2)        +      X3 * (F2-F1)
C
2     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'MT11AA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT11AA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(0)
        ENDIF
C
C     ELSEIF(IELMF.EQ. ) THEN
C     AUTRES TYPES DE FONCTIONS F
C
C-----------------------------------------------------------------------
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT11AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT11AA (BIEF) :',/,
     *        1X,'DISCRETISATION OF F: ',1I6,' NOT AVAILABLE')
       CALL PLANTE(0)
       STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
