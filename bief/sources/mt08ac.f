C                       *****************
                        SUBROUTINE MT08AC
C                       *****************
C    
     *(  A11 , A12 , A13 , A14 , A15, A16,
     *   A21 , A22 , A23 , A24 , A25, A26,
     *   A31 , A32 , A33 , A34 , A35, A36,
     *   XMUL,SF,F,XEL,YEL,IKLE1,IKLE2,IKLE3,
     *                     IKLE4,IKLE5,IKLE6,
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
C                  /                     D
C A(I,J)= -XMUL   /  PSI2(J) *    F    * --( PSI1(I) ) D(OMEGA)
C                /OMEGA                  DX
C
C  ATTENTION AU SIGNE MOINS !!
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE P1
C  PSI2 : BASES DE TYPE TRIANGLE P2
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
      USE BIEF!, EX_MT08AC => MT08AC
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
C
C-----------------------------------------------------------------------
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
        A12(IELEM) =  (Y3-Y2) * (-F1+2.D0*F2-F3) * (XMUL/120.D0)
        A13(IELEM) = -(Y3-Y2) * (F1+F2-2.D0*F3)  * (XMUL/120.D0)
        A14(IELEM) =  (Y3-Y2) * (2.D0*F1+F3+2.D0*F2) * (XMUL/30.D0)
        A15(IELEM) =  (Y3-Y2) * (F1+2.D0*F3+2.D0*F2) * (XMUL/30.D0)
        A16(IELEM) =  (Y3-Y2) * (F2+2.D0*F1+2.D0*F3) * (XMUL/30.D0)
        A21(IELEM) =  Y3      * (F2-2.D0*F1+F3)  * (XMUL/120.D0)
        A23(IELEM) = -Y3      * (-F1-F2+2.D0*F3) * (XMUL/120.D0)
        A24(IELEM) = -Y3      * (2.D0*F1+F3+2.D0*F2) * (XMUL/30.D0) 
        A25(IELEM) = -Y3      * (F1+2.D0*F3+2.D0*F2) * (XMUL/30.D0)
        A26(IELEM) = -Y3      * (F2+2.D0*F1+2.D0*F3) * (XMUL/30.D0)
        A31(IELEM) = -Y2      * (F2-2.D0*F1+F3)      * (XMUL/120.D0) 
        A32(IELEM) =  Y2      * (-F1+2.D0*F2-F3)     * (XMUL/120.D0)
        A34(IELEM) =  Y2      * (2.D0*F1+F3+2.D0*F2) * (XMUL/30.D0)
        A35(IELEM) =  Y2      * (F1+2.D0*F3+2.D0*F2) * (XMUL/30.D0)
        A36(IELEM) =  Y2      * (F2+2.D0*F1+2.D0*F3) * (XMUL/30.D0)
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
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
        A12(IELEM) =  (X3-X2 ) * (F1-2.D0*F2+F3 ) * (XMUL/120.D0)
        A13(IELEM) =  (-X3+X2) * (-F1-F2+2.D0*F3) * (XMUL/120.D0)
        A14(IELEM) =  (-X3+X2) * (2.D0*F1+F3+2.D0*F2) * (XMUL/30.D0)
        A15(IELEM) =  (-X3+X2) * (F1+2.D0*F3+2.D0*F2) * (XMUL/30.D0)
        A16(IELEM) =  (-X3+X2) * (2.D0*F3+2.D0*F1+F2) * (XMUL/30.D0)
        A21(IELEM) = -X3       * (F3-2.D0*F1+F2 ) * (XMUL/120.D0)
        A23(IELEM) =  X3       * (-F1-F2+2.D0*F3) * (XMUL/120.D0)
        A24(IELEM) =  X3       * (2.D0*F1+F3+2.D0*F2) * (XMUL/30.D0) 
        A25(IELEM) =  X3       * (F1+2.D0*F3+2.D0*F2) * (XMUL/30.D0)
        A26(IELEM) =  X3       * (2.D0*F3+2.D0*F1+F2) * (XMUL/30.D0)
        A31(IELEM) =  X2       * (F3-2.D0*F1+F2 ) * (XMUL/120.D0) 
        A32(IELEM) =  X2       * (F1-2.D0*F2+F3 ) * (XMUL/120.D0)
        A34(IELEM) = -X2       * (2.D0*F1+F3+2.D0*F2) * (XMUL/30.D0)
        A35(IELEM) = -X2       * (F1+2.D0*F3+2.D0*F2) * (XMUL/30.D0)
        A36(IELEM) = -X2       * (2.D0*F3+2.D0*F1+F2) * (XMUL/30.D0)
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
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
          CALL PLANTE(0)
          STOP
        ENDIF
C
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
       A12(IELEM) = (-Y3+Y2) *(F1-6.D0*F2+F3+4.D0*F6) * (XMUL/360.D0) 
       A13(IELEM) = (-Y3+Y2) *(F1+F2-6.D0*F3+4.D0*F4) * (XMUL/360.D0)
       A14(IELEM) = (-Y3+Y2) *(F3-8.D0*F4-4.D0*F6-4.D0*F5) *(XMUL/90.D0)
       A15(IELEM) = (-Y3+Y2) *(F1-4.D0*F4-4.D0*F6-8.D0*F5) *(XMUL/90.D0)
       A16(IELEM) = (F2-4.D0*F4-8.D0*F6-4.D0*F5) *(-Y3+Y2) *(XMUL/90.D0)
       A21(IELEM) =-Y3       *(6.D0*F1-F2-F3-4.D0*F5) * (XMUL/360.D0)
       A23(IELEM) = Y3       *(F1+F2-6.D0*F3+4.D0*F4) * (XMUL/360.D0)
       A24(IELEM) = Y3       *(F3-8.D0*F4-4.D0*F6-4.D0*F5) *(XMUL/90.D0)
       A25(IELEM) = Y3       *(F1-4.D0*F4-4.D0*F6-8.D0*F5) *(XMUL/90.D0)
       A26(IELEM) = Y3       *(F2-4.D0*F4-8.D0*F6-4.D0*F5) *(XMUL/90.D0)
       A31(IELEM) = Y2       *(6.D0*F1-F2-F3-4.D0*F5) * (XMUL/360.D0)
       A32(IELEM) =-Y2       *(F1-6.D0*F2+F3+4.D0*F6) * (XMUL/360.D0)
       A34(IELEM) =-Y2       *(F3-8.D0*F4-4.D0*F6-4.D0*F5) *(XMUL/90.D0)
       A35(IELEM) =-Y2*(F1-4.D0*F4-4.D0*F6-8.D0*F5) *( XMUL/90.D0)
       A36(IELEM) =-Y2*(F2-4.D0*F4-8.D0*F6-4.D0*F5) * (XMUL/90.D0)
C
C   TERMES DIAGONAUX
C   LA SOMME DES LIGNES DE LA MATRICE EST LE VECTEUR NUL
C
        A11(IELEM) = - A21(IELEM) - A31(IELEM)
        A22(IELEM) = - A12(IELEM) - A32(IELEM)
        A33(IELEM) = - A13(IELEM) - A23(IELEM)
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
        A12(IELEM) = (X3-X2) *(F1-6.D0*F2+F3+4.D0*F6) * (XMUL/360.D0)
        A13(IELEM) = (X3-X2) *(F1+F2-6.D0*F3+4.D0*F4) * (XMUL/360.D0)
        A14(IELEM) = (X3-X2) *(F3-8.D0*F4-4.D0*F5-4.D0*F6) *(XMUL/90.D0)
        A15(IELEM) = (X3-X2) *(F1-4.D0*F4-8.D0*F5-4.D0*F6) *(XMUL/90.D0)
        A16(IELEM) = (X3-X2) *(F2-4.D0*F4-4.D0*F5-8.D0*F6) *(XMUL/90.D0)
        A21(IELEM) = X3      *(6.D0*F1-F2-F3-4.D0*F5) * (XMUL/360.D0)
        A23(IELEM) =-X3      *(F1+F2-6.D0*F3+4.D0*F4) * (XMUL/360.D0)
        A24(IELEM) =-X3 *(F3-8.D0*F4-4.D0*F5-4.D0*F6) * (XMUL/90.D0)
        A25(IELEM) =-X3 *(F1-4.D0*F4-8.D0*F5-4.D0*F6) * (XMUL/90.D0)
        A26(IELEM) =-X3 *(F2-4.D0*F4-4.D0*F5-8.D0*F6) * (XMUL/90.D0)
        A31(IELEM) =-X2 *(6.D0*F1-F2-F3-4.D0*F5) * (XMUL/360.D0)
        A32(IELEM) = X2 *(F1-6.D0*F2+F3+4.D0*F6) * (XMUL/360.D0)
        A34(IELEM) = X2 *(F3-8.D0*F4-4.D0*F5-4.D0*F6) * (XMUL/90.D0)
        A35(IELEM) = X2 *(F1-4.D0*F4-8.D0*F5-4.D0*F6) * (XMUL/90.D0)
        A36(IELEM) = X2 *(F2-4.D0*F4-4.D0*F5-8.D0*F6) * (XMUL/90.D0)
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
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF
100    FORMAT(1X,'MT08AC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE',
     *        1X,'NOM REEL DE F : ',A6)
101    FORMAT(1X,'MT08AC (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
C
200       FORMAT(1X,'MT08AC (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT08AC (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
C
C-----------------------------------------------------------------------
C
      RETURN
      END
