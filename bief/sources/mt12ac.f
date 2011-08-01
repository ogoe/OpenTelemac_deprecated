C                       *****************
                        SUBROUTINE MT12AC
C                       *****************
C
     *(  A11 , A12 , A13 , A14 , A15, A16,
     *   A21 , A22 , A23 , A24 , A25, A26,
     *   A31 , A32 , A33 , A34 , A35, A36,
     *   XMUL,SF,SU,SV,F,U,V,
     *   XEL,YEL,SURFAC,
     *   IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,
     *   NELEM,NELMAX,ICOORD)
C
C***********************************************************************
C BIEF VERSION 5.9   19/06/08   ALGIANE FROEHLY (MATMECA) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  EXEMPLE AVEC ICOORD = 1
C
C
C                /           D       -> -->
C A(I,J)= XMUL  /  PSI2(J) * --(F) * U.GRAD(PSI1(I))  D(OMEGA)
C              /OMEGA        DX
C
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE P1
C  PSI2 : BASES DE TYPE P2
C  F    : FONCTION DE TYPE TRIANGLE P1
C  U    : VECTEUR DE TYPE P0 OU P2
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
C |     IKLE1..6   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
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
      USE BIEF !, EX_MT12AC => MT12AC
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
      
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: F(*),U(*),V(*)
C
C     STRUCTURES DE F,U,V
      TYPE(BIEF_OBJ), INTENT(IN) :: SF,SU,SV
C
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF,IELMU,IELMV
      DOUBLE PRECISION X2,X3,Y2,Y3,F1,F2,F3
      DOUBLE PRECISION U1,U2,U3,U4,U5,U6,V1,V2,V3,V4,V5,V6,UX,UY
      DOUBLE PRECISION XSUR12 ,XSUR180,XSUR720 
      DOUBLE PRECISION AUX12  ,AUX180 ,AUX720 , UNSURF 
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
      IELMU=SU%ELM
      IELMV=SV%ELM
C
      XSUR12  = XMUL /  12.D0
      XSUR720 = XMUL / 720.D0
      XSUR180 = XMUL / 180.D0
C
C-----------------------------------------------------------------------
C  CAS OU F EST DE DISCRETISATION P1
C-----------------------------------------------------------------------
C
      IF(IELMF.EQ.11) THEN
C
      IF(IELMU.EQ.10.AND.IELMV.EQ.10) THEN
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
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        UX  =  U(IELEM)
        UY  =  V(IELEM)
C
        UNSURF= 1.D0 / SURFAC(IELEM)
        AUX12 = XSUR12 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)= 0.D0
        A13(IELEM)= 0.D0
        A21(IELEM)= 0.D0
        A23(IELEM)= 0.D0
        A24(IELEM)= (Y3*UX-UY*X3) * (Y3*F2-Y2*F3) * AUX12
        A25(IELEM)= A24(IELEM)
        A26(IELEM)= A24(IELEM)
        A31(IELEM)= 0.D0
        A32(IELEM)= 0.D0
        A34(IELEM)= (X2*UY-UX*Y2) * (Y3*F2-Y2*F3) * AUX12
        A35(IELEM)= A34(IELEM)
        A36(IELEM)= A34(IELEM)        
C 
C   LA SOMME DE CHAQUE COLONNE EST NULLE
C
        A14(IELEM)= - A24(IELEM) - A34(IELEM)
        A15(IELEM)= A14(IELEM)
        A16(IELEM)= A14(IELEM)
C 
CC  ELEMENTS DIAGONAUX
C
        A11(IELEM) = 0.D0
        A22(IELEM) = 0.D0
        A33(IELEM) = 0.D0
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
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        UX  =  U(IELEM)
        UY  =  V(IELEM)
C
        UNSURF = 1.D0 / SURFAC(IELEM)
        AUX12 = XSUR12 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)= 0.D0
        A13(IELEM)= 0.D0
        A21(IELEM)= 0.D0
        A23(IELEM)= 0.D0
        A24(IELEM)= (X3*F2-X2*F3) * (UY*X3-UX*Y3) * AUX12
        A25(IELEM)=  A24(IELEM)
        A26(IELEM)=  A24(IELEM)
        A31(IELEM)= 0.D0
        A32(IELEM)= 0.D0
        A34(IELEM)= (X3*F2-X2*F3) * (UX*Y2-UY*X2) * AUX12
        A35(IELEM)=  A34(IELEM)
        A36(IELEM)=  A34(IELEM)
C
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A14(IELEM)= - A24(IELEM) - A34(IELEM)
        A15(IELEM)=   A14(IELEM)
        A16(IELEM)=   A14(IELEM)
C
CC   TERMES DIAGONAUX
C
        A11(IELEM) = 0.D0
        A22(IELEM) = 0.D0
        A33(IELEM) = 0.D0
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
      ELSEIF(IELMU.EQ.13) THEN
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
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        U1  =  U(IKLE1(IELEM))
        U2  =  U(IKLE2(IELEM))
        U3  =  U(IKLE3(IELEM))
        U4  =  U(IKLE4(IELEM))
        U5  =  U(IKLE5(IELEM))
        U6  =  U(IKLE6(IELEM))
        V1  =  V(IKLE1(IELEM))
        V2  =  V(IKLE2(IELEM))
        V3  =  V(IKLE3(IELEM))
        V4  =  V(IKLE4(IELEM))
        V5  =  V(IKLE5(IELEM))
        V6  =  V(IKLE6(IELEM))
C
        UNSURF = 1.D0 / SURFAC(IELEM)
        AUX720 = XSUR720 * UNSURF
        AUX180 = XSUR180 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)= (Y3*F2-Y2*F3) * ((U3+U1-6.D0*U2+4.D0*U6) * (Y3-Y2)
     *            + (V3+V1-6.D0*V2+4.D0*V6) * (X2-X3)) * AUX720
C
        A13(IELEM)= (Y3*F2-Y2*F3) * ((4.D0*V4-6.D0*V3+V2+V1) * (X2-X3)
     *            + (6.D0*U3-4.D0*U4-U2-U1) * (Y2-Y3)) * AUX720
C       
        A21(IELEM)= (Y3*F2-Y2*F3) * ((4.D0*V5+V3+V2-6.D0*V1) * X3
     *            - (U3-6.D0*U1+4.D0*U5+U2) * Y3     ) * AUX720
C        
        A23(IELEM)= (Y3*F2-Y2*F3) * ((V1-6.D0*V3+4.D0*V4+V2) * X3
     *            + (6.D0*U3-4.D0*U4-U2-U1) * Y3     ) * AUX720
C     
        A24(IELEM)= (Y3*F2-Y2*F3) * ((V3-4.D0*V6-8.D0*V4-4.D0*V5) * X3
     *            - (U3-8.D0*U4-4.D0*U6-4.D0*U5) * Y3) * AUX180
C        
        A25(IELEM)= (Y3*F2-Y2*F3) * ((V1-4.D0*V6-8.D0*V5-4.D0*V4) * X3
     *            + (4.D0*U6-U1+4.D0*U4+8.D0*U5) * Y3) * AUX180
C     
        A26(IELEM)= (Y3*F2-Y2*F3) * ((V2-4.D0*V4-4.D0*V5-8.D0*V6) * X3
     *            + (8.D0*U6+4.D0*U4+4.D0*U5-U2) * Y3) * AUX180
C       
        A31(IELEM)= (Y3*F2-Y2*F3) * ((6.D0*V1-4.D0*V5-V3-V2) * X2
     *            + (U3-6.D0*U1+4.D0*U5+U2) * Y2         ) * AUX720
C     
        A32(IELEM)= (Y3*F2-Y2*F3) * ((-4.D0*V6+6.D0*V2-V1-V3) * X2
     *            + (U3+U1+4.D0*U6-6.D0*U2) * Y2        ) * AUX720
C       
        A34(IELEM)= (Y3*F2-Y2*F3) * ((4.D0*V6+8.D0*V4+4.D0*V5-V3) * X2
     *            + (U3-8.D0*U4-4.D0*U6-4.D0*U5) * Y2    ) * AUX180
C      
        A35(IELEM)= (Y2*F3-Y3*F2) * ((V1-8.D0*V5-4.D0*V4-4.D0*V6) * X2
     *            + (4.D0*U6-U1+4.D0*U4+8.D0*U5) * Y2   ) * AUX180
C 
        A36(IELEM)= (Y2*F3-Y3*F2) * ((V2-4.D0*V4-4.D0*V5-8.D0*V6) * X2
     *            + (8.D0*U6+4.D0*U4+4.D0*U5-U2) * Y2   ) * AUX180
C
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A14(IELEM) = - A24(IELEM) - A34(IELEM)
        A15(IELEM) = - A25(IELEM) - A35(IELEM)
        A16(IELEM) = - A26(IELEM) - A36(IELEM)
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
C
3     CONTINUE
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
        Y2  =  YEL(IELEM,2)
        Y3  =  YEL(IELEM,3)
C
        F1  =  F(IKLE1(IELEM))
        F2  =  F(IKLE2(IELEM)) - F1
        F3  =  F(IKLE3(IELEM)) - F1
C
        U1  =  U(IKLE1(IELEM))
        U2  =  U(IKLE2(IELEM))
        U3  =  U(IKLE3(IELEM))
        U4  =  U(IKLE4(IELEM))
        U5  =  U(IKLE5(IELEM))
        U6  =  U(IKLE6(IELEM))
                
        V1  =  V(IKLE1(IELEM))
        V2  =  V(IKLE2(IELEM))
        V3  =  V(IKLE3(IELEM))
        V4  =  V(IKLE4(IELEM))
        V5  =  V(IKLE5(IELEM))
        V6  =  V(IKLE6(IELEM))
C
        UNSURF= 1.D0 / SURFAC(IELEM)
        AUX180= XSUR180 * UNSURF
        AUX720= XSUR720 * UNSURF
C
C   TERMES EXTRADIAGONAUX
C
        A12(IELEM)= (X2*F3-X3*F2) * ((4.D0*V6-6.D0*V2+V3+V1) * (X2-X3)
     *            + (U1+U3+4.D0*U6-6.D0*U2) * (Y3-Y2)) * AUX720
C       
        A13(IELEM)= (X2*F3-X3*F2) * ((V2-6.D0*V3+4.D0*V4+V1) * (X2-X3)
     *            + (6.D0*U3-U1-U2-4.D0*U4) * (Y2-Y3)) * AUX720
C
        A21(IELEM)= (X3*F2-X2*F3) * ((6.D0*V1-V3-V2-4.D0*V5) * X3
     *            + (U3-6.D0*U1+U2+4.D0*U5) * Y3      ) * AUX720
C
        A23(IELEM)= (X2*F3-X3*F2) * ((V2-6.D0*V3+4.D0*V4+V1) * X3
     *            + (6.D0*U3-U1-U2-4.D0*U4) * Y3      ) * AUX720
C
        A24(IELEM)= (X2*F3-X3*F2) * ((V3-4.D0*V6-4.D0*V5-8.D0*V4) * X3
     *            + (4.D0*U5+8.D0*U4-U3+4.D0*U6) * Y3 ) * AUX180
C       
        A25(IELEM)= (X2*F3-X3*F2) * ((V1-8.D0*V5-4.D0*V4-4.D0*V6) * X3
     *            + (4.D0*U4-U1+8.D0*U5+4.D0*U6) * Y3 ) * AUX180
C
        A26(IELEM)= (X2*F3-X3*F2) * ((V2-4.D0*V5-4.D0*V4-8.D0*V6) * X3
     *            + (4.D0*U5-U2+8.D0*U6+4.D0*U4) * Y3 ) * AUX180
C
        A31(IELEM)= (X2*F3-X3*F2) * ((6.D0*V1-V3-V2-4.D0*V5) * X2
     *            + (U3-6.D0*U1+U2+4.D0*U5) * Y2      ) * AUX720
C
        A32(IELEM)= (X3*F2-X2*F3) * ((4.D0*V6-6.D0*V2+V3+V1) * X2
     *            + (6.D0*U2-U1-U3-4.D0*U6) * Y2      ) * AUX720
C
        A34(IELEM)= (X3*F2-X2*F3) * ((V3-4.D0*V6-4.D0*V5-8.D0*V4) * X2
     *            + (4.D0*U5+8.D0*U4-U3+4.D0*U6) * Y2 ) * AUX180
C
        A35(IELEM)= (X3*F2-X2*F3) * ((V1-8.D0*V5-4.D0*V4-4.D0*V6) * X2
     *            + (4.D0*U4-U1+8.D0*U5+4.D0*U6) * Y2 ) * AUX180
C
        A36(IELEM)= (X3*F2-X2*F3) * ((V2-4.D0*V5-4.D0*V4-8.D0*V6) * X2
     *            + (4.D0*U5-U2+8.D0*U6+4.D0*U4) * Y2 ) * AUX180
C
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A14(IELEM) = - A24(IELEM) - A34(IELEM)
        A15(IELEM) = - A25(IELEM) - A35(IELEM)
        A16(IELEM) = - A26(IELEM) - A36(IELEM)
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
C
4       CONTINUE
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,300) IELMU
        IF (LNG.EQ.2) WRITE(LU,301) IELMU
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
       IF (LNG.EQ.1) WRITE(LU,100) IELMF
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT12AC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT12AC (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
C
200       FORMAT(1X,'MT12AC (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT12AC (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
300    FORMAT(1X,'MT12AC (BIEF) :',/,
     *        1X,'DISCRETISATION DE U : ',1I6,' NON PREVUE')
301    FORMAT(1X,'MT12AC (BIEF) :',/,
     *        1X,'DISCRETIZATION OF U : ',1I6,' NOT AVAILABLE')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
