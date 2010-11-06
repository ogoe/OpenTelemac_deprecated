C                       *****************
                        SUBROUTINE MT11AC
C                       *****************
C
     *(  A11 , A12 , A13 , A14 , A15, A16,
     *   A21 , A22 , A23 , A24 , A25, A26,
     *   A31 , A32 , A33 , A34 , A35, A36,
     *   XMUL,SF,F,XEL,YEL,IKLE1,IKLE2,IKLE3,
     *   IKLE4,IKLE5,IKLE6,
     *   NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.9        01/07/08    A FROEHLY (MATMECA) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  EXEMPLE AVEC ICOORD = 1
C
C            NPOIN
C              _          /            D
C A(I,J)=-XMUL>_   F  *  /  PSI2(J) *  --( PSI1(K) * PSI1(I) ) D(OMEGA)
C                   K   /OMEGA         DX
C             K=1
C
C  ATTENTION AU SIGNE MOINS ||
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE P1
C  PSI2 : BASES DE TYPE TRIANGLE P2
C
C-----------------------------------------------------------------------
C
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
      USE BIEF!, EX_MT11AC => MT11AC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX,ICOORD
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
      INTEGER, INTENT(IN) :: IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A14(*),A15(*),A16(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A21(*),A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A24(*),A25(*),A26(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A31(*),A32(*),A33(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A34(*),A35(*),A36(*)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C     STRUCTURE DE F
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION X2,X3,Y2,Y3,F1,F2,F3,F4,F5,F6
      DOUBLE PRECISION XSUR30,XSUR60,XSUR120
      DOUBLE PRECISION XSUR90,XSUR360
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
      XSUR30  = XMUL/30.D0
      XSUR60  = XMUL/60.D0
      XSUR120 = XMUL/120.D0
      XSUR90  = XMUL/90.D0
      XSUR360 = XMUL/360.D0
C
      IELMF=SF%ELM
C
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION P1
C-----------------------------------------------------------------------
C
      IF(IELMF.EQ.11) THEN
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
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))        
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM) = ( 2.D0*(F1-F2)*Y2 + (3.D0*F2-2.D0*F1-F3)*Y3) 
     *             * XSUR120 
        A13(IELEM) = ( 2.D0*(F3-F1)*Y3 + (2.D0*F1-3.D0*F3+F2)*Y2) 
     *             * XSUR120 
        A14(IELEM) = (( 4.D0*F1+F3)*Y3 + (F3-4.D0*F1-2.D0*F2)*Y2) 
     *             * XSUR30 
        A15(IELEM) = ((-2.D0*(F1+F2)-F3)*Y2  + 
     *                ( 2.D0*(F1+F3)+F2)*Y3) * XSUR30 
        A16(IELEM) = ((-4.D0*F1-F2)*Y2 + (4.D0*F1-F2+2.D0*F3)*Y3) 
     *             * XSUR30 
        A21(IELEM) = (      (F1-F3)*Y2 + (F3-3.D0*F1+2.D0*F2)*Y3) 
     *             * XSUR120 
        A23(IELEM) = (      2.D0*(F2-F3)*Y3 +        (F1-F3) *Y2) 
     *             * XSUR120 
        A24(IELEM) = ( 2.D0*(F3-F1)*Y2 - (        F3+4.D0*F2)*Y3) 
     *             * XSUR30 
        A25(IELEM) = ( 2.D0*(F3-F1)*Y2 + (F1-2.D0*F3-4.D0*F2)*Y3) 
     *             * XSUR30 
        A26(IELEM) = (      (F3-F1)*Y2 - (F1+2.D0*F2+2.D0*F3)*Y3) 
     *             * XSUR30 
        A31(IELEM) = (      (F2-F1)*Y3 + (3.D0*F1-F2-2.D0*F3)*Y2) 
     *             * XSUR120 
        A32(IELEM) = ( 2.D0*(F2-F3)*Y2 +             (-F1+F2)*Y3) 
     *             * XSUR120 
        A34(IELEM) = (      (F1-F2)*Y3 + (   F1+2.D0*(F2+F3))*Y2) 
     *             * XSUR30 
        A35(IELEM) = ( 2.D0*(F1-F2)*Y3 + (4.D0*F3-F1+2.D0*F2)*Y2) 
     *             * XSUR30 
        A36(IELEM) = (( 4.D0*F3+F2)*Y2 +         2.D0*(F1-F2)*Y3) 
     *             * XSUR30  
C
C   TERMES DIAGONAUX
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
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
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM) = (  2.D0*(F2-F1)*X2 + (2.D0*F1-3.D0*F2+F3)*X3) 
     *             * XSUR120 
        A13(IELEM) = (  2.D0*(F1-F3)*X3 - (2.D0*F1-3.D0*F3+F2)*X2) 
     *             * XSUR120
        A14(IELEM) = ((-4.D0*F1-F3 )*X3 + (4.D0*F1-F3+2.D0*F2)*X2) 
     *             * XSUR30
        A15(IELEM) = (( 2.D0*(F1+F2)+F3)*X2  - 
     *                ( 2.D0*(F1+F3)+F2)*X3) * XSUR30 
        A16(IELEM) = ((4.D0*F1+F2  )*X2 + (F2-4.D0*F1-2.D0*F3)*X3) 
     *             * XSUR30
        A21(IELEM) = (       (F3-F1)*X2 + (3.D0*F1-2.D0*F2-F3)*X3) 
     *             * XSUR120
        A23(IELEM) = (       (F3-F1)*X2 +         2.D0*(F3-F2)*X3) 
     *             * XSUR120 
        A24(IELEM) = (  2.D0*(F1-F3)*X2 + (        F3+4.D0*F2)*X3) 
     *             * XSUR30 
        A25(IELEM) = (  2.D0*(F1-F3)*X2 + (2.D0*F3-F1+4.D0*F2)*X3) 
     *             * XSUR30
        A26(IELEM) = (       (F1-F3)*X2 + (   F1+2.D0*(F2+F3))*X3) 
     *             * XSUR30
        A31(IELEM) = (       (F1-F2)*X3 + (F2+2.D0*F3-3.D0*F1)*X2) 
     *             * XSUR120  
        A32(IELEM) = (  2.D0*(F3-F2)*X2 + (             F1-F2)*X3) 
     *             * XSUR120
        A34(IELEM) = (       (F2-F1)*X3 - (F1+2.D0*F2+2.D0*F3)*X2) 
     *             * XSUR30
        A35(IELEM) = (  2.D0*(F2-F1)*X3 - (4.D0*F3-F1+2.D0*F2)*X2) 
     *             * XSUR30
        A36(IELEM) = ((-4.D0*F3-F2 )*X2 +         2.D0*(F2-F1)*X3) 
     *             * XSUR30
C
C   TERMES DIAGONAUX
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
C
2       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION P2
C-----------------------------------------------------------------------
C
       ELSEIF(IELMF.EQ.13) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT X  =
C================================
C
        IF(ICOORD.EQ.1) THEN
C
C   BOUCLE SUR LES ELEMENTS
C
        DO 3 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
        F4  =  F(IKLE4(IELEM))
        F5  =  F(IKLE5(IELEM))
        F6  =  F(IKLE6(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
       A12(IELEM) = ((3.D0*F2-6.D0*F1-8.D0*(F6-F4)-F3+4.D0*F5)*Y3 
     *            +   6.D0*(F1-F2)*Y2 ) * XSUR360
       A13(IELEM) = ((8.D0*(F4-F6)-4.D0*F5+F2+6.D0*F1-3.D0*F3)*Y2
     *            +   6.D0*(F3-F1)*Y3 ) * XSUR360
       A14(IELEM) = ((-16.D0*F4+4.D0*(F5+F6)-6.D0*F1-F3      )*Y2
     *            +  (6.D0*F1-2.D0*F2+4.D0*F4+8.D0*F6-F3)*Y3) 
     *             * XSUR90
       A15(IELEM) = ((F3-8.D0*F4-4.D0*(F5+F6))*Y2
     *            +  (4.D0*(F4+F5)+8.D0*F6-F2)*Y3) * XSUR90
       A16(IELEM) = ((F2+6.D0*F1+16.D0*F6-4.D0*(F4+F5)       )*Y3
     *            -  (8.D0*F4-F2+6.D0*F1-2.D0*F3+4.D0*F6)*Y2) 
     *             * XSUR90
       A21(IELEM) = ((8.D0*(F4-F5)-3.D0*F1-F3+4.D0*F6        )*Y2
     *            - (3.D0*F1-6.D0*F2+8.D0*(F4-F5)+4.D0*F6-F3 )*Y3) 
     *            * XSUR360
       A23(IELEM) = ((8.D0*(F4-F5)+F1+3.D0*F3-4.D0*F6           )*Y2
     *            +  6.D0*(F2-F3)*Y3 ) * XSUR360
       A24(IELEM) = ((12.D0*(F5-F4)-2.D0*(F1+F3)+4.D0*F6     )*Y2
     *            + (2.D0*F1-4.D0*F4+F3-6.D0*F2-8.D0*F5)*Y3 ) 
     *            * XSUR90
       A25(IELEM) = ((12.D0*(F5-F4)+2.D0*(F1+F3)-4.D0*F6     )*Y2
     *            +  (-F1+4.D0*(F4+F6)-6.D0*F2-16.D0*F5)*Y3 ) 
     *            * XSUR90
       A26(IELEM) = ((F1-8.D0*F5-4.D0*(F6+F4)                )*Y3 
     *            - (4.D0*(F4-F5)+F1-F3)*Y2                 ) 
     *            * XSUR90
       A31(IELEM) = ((4.D0*F4-8.D0*(F5-F6)+3.D0*F1-F2-6.D0*F3)*Y2
     *            +  (3.D0*F1-8.D0*(F6-F5)-4.D0*F4+F2)*Y3 ) 
     *            * XSUR360
       A32(IELEM) = ((4.D0*F4-F1-8.D0*(F6-F5)-3.D0*F2)*Y3
     *            +  6.D0*(F2-F3)*Y2) * XSUR360
       A34(IELEM) = ((4.D0*(F4+F6)-F1+8.D0*F5)*Y2
     *            +  (F1-4.D0*(F5-F6)-F2)*Y3 ) * XSUR90
       A35(IELEM) = ((16.D0*F5-4.D0*(F4+F6)+F1+6.D0*F3)*Y2
     *            +  (-2.D0*(F1+F2)+4.D0*F4+12.D0*(F6-F5))*Y3) 
     *            * XSUR90
       A36(IELEM) = ((8.D0*F5-F2-2.D0*F1+6.D0*F3+4.D0*F6)*Y2
     *            +  (2.D0*(F2+F1)+12.D0*(F6-F5)-4.D0*F4)*Y3) 
     *            * XSUR90
C
C   TERMES DIAGONAUX
C
        A11(IELEM) = ((24.D0*(F6-F1)-5.D0*F3+4.D0*F5+F2)*Y2
     *             +  (24.D0*(F1-F4)+5.D0*F2-4.D0*F5-F3)*Y3) * XSUR360
        A22(IELEM) = ((6.D0*(F1-F3)+24.D0*(F5-F4)      )*Y2
     *             +  (24.D0*(F4-F2)+4.D0*F6+F3-5.D0*F1)*Y3) * XSUR360
        A33(IELEM) = ((24.D0*(F3-F6)+5.D0*F1-4.D0*F4-F2)*Y2
     *             +  (6.D0*(F2-F1)+24.D0*(F6-F5)      )*Y3) * XSUR360
C
3       CONTINUE
C
        ELSEIF(ICOORD.EQ.2) THEN
C
C================================
C  CAS DE LA DERIVEE SUIVANT Y  =
C================================
C
        DO 4 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
        X2  =  XEL(IELEM,2)
        X3  =  XEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM))
        F3  =  F(IKLE3(IELEM))
        F4  =  F(IKLE4(IELEM))
        F5  =  F(IKLE5(IELEM))
        F6  =  F(IKLE6(IELEM))
C
C   TERMES EXTRADIAGONAUX
C
       A12(IELEM) = ((F3+8.D0*(F6-F4)+6.D0*F1-3.D0*F2-4.D0*F5)*X3 
     *            +   6.D0*(F2-F1)*X2 ) * XSUR360
       A13(IELEM) = ((8.D0*(F6-F4)-6.D0*F1+3.D0*F3+4.D0*F5-F2)*X2
     *            -   6.D0*(F3-F1)*X3 ) * XSUR360
       A14(IELEM) = ((6.D0*F1+F3-4.D0*(F6+F5)+16.D0*F4  )*X2
     *            +  (F3-8.D0*F6-6.D0*F1+2.D0*F2-4.D0*F4)*X3) 
     *            *  XSUR90
       A15(IELEM) = ((4.D0*(F6+F5)+8.D0*F4-F3)*X2
     *            +  (F2-8.D0*F6-4.D0*(F4+F5))*X3 ) * XSUR90
       A16(IELEM) = ((6.D0*F1-2.D0*F3+4.D0*F6+8.D0*F4-F2)*X2
     *            +  (-16.D0*F6-6.D0*F1-F2+4.D0*(F5+F4) )*X3) 
     *             * XSUR90 
       A21(IELEM) = ((3.D0*F1+F3-4.D0*F6-8.D0*(F4-F5))*X2
     *            +  (-F3+4.D0*F6+3.D0*F1-6.D0*F2+8.D0*(F4-F5))*X3)
     *            *  XSUR360
       A23(IELEM) = ((4.D0*F6-F1-3.D0*F3-8.D0*(F4-F5))*X2
     *            -   6.D0*(F2-F3)*X3 ) * XSUR360
       A24(IELEM) = ((2.D0*(F1+F3)-4.D0*F6+12.D0*(F4-F5) )*X2
     *            +  (-F3-2.D0*F1+6.D0*F2+4.D0*F4+8.D0*F5)*X3) 
     *            * XSUR90 
       A25(IELEM) = ((4.D0*F6-2.D0*(F1+F3)+12.D0*(F4-F5) )*X2
     *            -  (4.D0*F6-F1-6.D0*F2+4.D0*F4-16.D0*F5)*X3) 
     *            * XSUR90 
       A26(IELEM) = ((F1-F3+4.D0*(F4-F5)     )*X2 
     *            +  (4.D0*(F6+F4)-F1+8.D0*F5)*X3 ) 
     *            * XSUR90 
       A31(IELEM) = ((F2-3.D0*F1+6.D0*F3-8.D0*(F6-F5)-4.D0*F4)*X2 
     *            -  (-8.D0*F6+3.D0*F1+F2-4.D0*F4+8.D0*F5    )*X3)
     *            * XSUR360 
       A32(IELEM) = ((F1-8.D0*(F5-F6)+3.D0*F2-4.D0*F4)*X3
     *            +   6.D0*(F3-F2)*X2 ) * XSUR360
       A34(IELEM) = ((F1-4.D0*(F6+F4)-8.D0*F5)*X2
     *            +  (-4.D0*(F6-F5)-F1+F2    )*X3 ) * XSUR90
       A35(IELEM) = ((4.D0*(F6+F4)-F1-6.D0*F3-16.D0*F5  )*X2
     *            -  (12.D0*(F6-F5)-2.D0*(F1+F2)+4.D0*F4)*X3) 
     *            * XSUR90
       A36(IELEM) = ((2.D0*F1-6.D0*F3-4.D0*F6-8.D0*F5+F2)*X2
     *            +  (12.D0*(F5-F6)-2.D0*(F1+F2)+4.D0*F4)*X3)  
     *            * XSUR90
C
C   TERMES DIAGONAUX
C
        A11(IELEM) = ((24.D0*(F1-F6)-4.D0*F5+5.D0*F3-F2)*X2
     *             +  (24.D0*(F4-F1)-5.D0*F2+4.D0*F5+F3)*X3)*XSUR360
        A22(IELEM) = ((24.D0*(F4-F5)+6.D0*(F3-F1)      )*X2
     *             +  (5.D0*F1+24.D0*(F2-F4)-4.D0*F6-F3)*X3)*XSUR360
        A33(IELEM) = ((4.D0*F4-5.D0*F1+24.D0*(F6-F3)+F2)*X2 
     *             +  (6.D0*(F1-F2)+24.D0*(F5-F6)      )*X3)*XSUR360
C
4       CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT11AC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT11AC (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
C
200       FORMAT(1X,'MT11AC (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT11AC (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
