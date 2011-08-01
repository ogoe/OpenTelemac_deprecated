C                       *****************
                        SUBROUTINE MT06CC
C                       *****************
C
     *( A11 , A12 , A13 , A14 , A15 , A16 ,
     *        A22 , A23 , A24 , A25 , A26 ,
     *              A33 , A34 , A35 , A36 ,
     *                    A44 , A45 , A46 ,
     *                          A55 , A56 ,
     *                                A66 ,
     *  XMUL,SF,F,SURFAC,
     *  IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.9   19/06/08   ALGIANE FROEHLY (MATMECA) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE SOMME(F*PSII*PSIJ)
C
C            AVEC :  P1 QUADRATIQUE
C                    P2 QUADRATIQUE
C                    F  P1 OU QUADRATIQUE
C
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
      USE BIEF!, EX_MT06CC => MT06CC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE1(NELMAX),IKLE2(NELMAX)
      INTEGER, INTENT(IN) :: IKLE3(NELMAX),IKLE4(NELMAX)
      INTEGER, INTENT(IN) :: IKLE5(NELMAX),IKLE6(NELMAX)
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A14(*),A15(*),A16(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A24(*),A25(*),A26(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A34(*),A35(*),A36(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A44(*),A45(*),A46(*)
      DOUBLE PRECISION, INTENT(INOUT) ::        A55(*),A56(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A66(*)
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
      DOUBLE PRECISION F1,F2,F3,F4,F5,F6,XSUR030,XSUR045,XSUR180
      DOUBLE PRECISION XSUR315,XSUR210,XSUR630,XSU1260
      DOUBLE PRECISION AUX315,AUX210,AUX630,AUX1260 
      INTEGER IELMF,IELEM
C
C=======================================================================
C
C     EXTRACTION DU TYPE D'ELEMENT DE F
C
      IELMF = SF%ELM
C
C  CAS OU F EST P0
C
      IF(IELMF.EQ.10) THEN
C
      XSUR030 = XMUL / 30.D0
      XSUR045 = XMUL / 45.D0
      XSUR180 = XMUL /180.D0
C
      DO 1 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
         F1  =  F(IELEM) * SURFAC(IELEM)
C
C  TERMES DIAGONAUX
C
         A11(IELEM) =   XSUR030 * F1
         A22(IELEM) =   A11(IELEM)
         A33(IELEM) =   A11(IELEM)
         A44(IELEM) =   8.D0 * XSUR045 * F1
         A55(IELEM) =   A44(IELEM) 
         A66(IELEM) =   A44(IELEM)
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = - XSUR180 * F1
         A13(IELEM) =   A12(IELEM)
         A14(IELEM) =   0.D0
         A15(IELEM) = - XSUR045*F1
         A16(IELEM) =   0.D0 
C
         A23(IELEM) =   A12(IELEM)
         A24(IELEM) =   0.D0
         A25(IELEM) =   0.D0
         A26(IELEM) =   A15(IELEM)
C        
         A34(IELEM) =   A15(IELEM)
         A35(IELEM) =   0.D0 
         A36(IELEM) =   0.D0
C
         A45(IELEM) =   4.D0*XSUR045*F1
         A46(IELEM) =   A45(IELEM)
C
         A56(IELEM) =   A45(IELEM)
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
C  CAS OU F EST LINEAIRE
C
      ELSEIF(IELMF.EQ.11) THEN
C
      XSUR210 = XMUL /  210.D0
      XSUR315 = XMUL /  315.D0
      XSU1260 = XMUL / 1260.D0
C
      DO 4 IELEM = 1 , NELEM

      AUX210 = SURFAC(IELEM) * XSUR210
      AUX315 = SURFAC(IELEM) * XSUR315
      AUX1260= SURFAC(IELEM) * XSU1260
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      F1  =  F(IKLE1(IELEM))
      F2  =  F(IKLE2(IELEM))
      F3  =  F(IKLE3(IELEM))
C  
C   INITIALISATION DES INTERMEDIAIRES DE CALCULS
C
C
C  TERMES DIAGONAUX
C
         A11(IELEM) =      (5.D0*F1+     F2+     F3)*AUX210
         A22(IELEM) =      (     F1+5.D0*F2+     F3)*AUX210
         A33(IELEM) =      (     F1+     F2+5.D0*F3)*AUX210
         A44(IELEM) = 8.D0*(3.D0*F1+3.D0*F2+     F3)*AUX315
         A55(IELEM) = 8.D0*(     F1+3.D0*F2+3.D0*F3)*AUX315
         A66(IELEM) = 8.D0*(3.D0*F1+     F2+3.D0*F3)*AUX315
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = -(4.D0*F1+4.D0*F2-     F3)*AUX1260
         A13(IELEM) = -(4.D0*F1-     F2+4.D0*F3)*AUX1260
         A14(IELEM) =  (3.D0*F1-2.D0*F2-     F3)*AUX315
         A15(IELEM) = -(     F1+3.D0*F2+3.D0*F3)*AUX315
         A16(IELEM) =  (3.D0*F1-     F2-2.D0*F3)*AUX315
C        
         A23(IELEM) =  (     F1-4.D0*F2-4.D0*F3)*AUX1260
         A24(IELEM) = -(2.D0*F1-3.D0*F2+     F3)*AUX315
         A25(IELEM) = -(     F1-3.D0*F2+2.D0*F3)*AUX315
         A26(IELEM) = -(3.D0*F1+     F2+3.D0*F3)*AUX315
C        
         A34(IELEM) = -(3.D0*F1+3.D0*F2+     F3)*AUX315
         A35(IELEM) = -(     F1+2.D0*F2-3.D0*F3)*AUX315
         A36(IELEM) = -(2.D0*F1+     F2-3.D0*F3)*AUX315
C   
         A45(IELEM) =  (2.D0*F1+3.D0*F2+2.D0*F3)*AUX315*4.D0
         A46(IELEM) =  (3.D0*F1+2.D0*F2+2.D0*F3)*AUX315*4.D0
C
         A56(IELEM) =  (2.D0*F1+3.D0*F3+2.D0*F2)*AUX315*4.D0
C
4     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.13) THEN
C
C-----------------------------------------------------------------------
C
C   DISCRETISATION QUADRATIQUE DE F :
C
      XSUR315 = XMUL /  315.D0
      XSUR630 = XMUL /  630.D0     
      XSU1260 = XMUL / 1260.D0
C
      DO 5 IELEM = 1 , NELEM

      AUX315 = SURFAC(IELEM) * XSUR315
      AUX630 = SURFAC(IELEM) * XSUR630 
      AUX1260= SURFAC(IELEM) * XSU1260
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
      F1  =  F(IKLE1(IELEM))
      F2  =  F(IKLE2(IELEM))
      F3  =  F(IKLE3(IELEM))
      F4  =  F(IKLE4(IELEM))
      F5  =  F(IKLE5(IELEM))
      F6  =  F(IKLE6(IELEM))
C
C  TERMES DIAGONAUX
C
       A11(IELEM) = (6.D0*(F4+F6)+9.D0*F1+2.D0*F5-F2-F3) * AUX630
       A22(IELEM) = (6.D0*(F4+F5)+2.D0*F6+9.D0*F2-F1-F3) * AUX630
       A33(IELEM) = (6.D0*(F6+F5)+9.D0*F3+2.D0*F4-F1-F2) * AUX630
       A44(IELEM) =  4.D0*(3.D0*(F6+F5)+9.D0*F4-F3) * AUX315
       A55(IELEM) =  4.D0*(3.D0*(F4+F6)+9.D0*F5-F1) * AUX315
       A66(IELEM) =  4.D0*(3.D0*(F4+F5)+9.D0*F6-F2) * AUX315 
C
C  TERMES EXTRADIAGONAUX
C
      A12(IELEM) = -(2.D0*(F1+F2)+4.D0*F4-F3) * AUX1260
      A13(IELEM) = -(2.D0*(F1+F3)+4.D0*F6-F2) * AUX1260
      A14(IELEM) =  (3.D0* F1    -2.D0*F5-F2) * AUX315
      A15(IELEM) = -(2.D0*(F4+F6)+4.D0*F5-F1) * AUX315
      A16(IELEM) =  (3.D0*F1     -2.D0*F5-F3) * AUX315
C
      A23(IELEM) =  (-2.D0*(F2+F3)-4.D0*F5+F1) * AUX1260
      A24(IELEM) =  (-2.D0*F6     +3.D0*F2-F1) * AUX315
      A25(IELEM) =  (-2.D0*F6     +3.D0*F2-F3) * AUX315
      A26(IELEM) =  (-2.D0*(F4+F5)-4.D0*F6+F2) * AUX315
C
      A34(IELEM) =  (-2.D0*(F6+F5)-4.D0*F4+F3) * AUX315
      A35(IELEM) =  (-2.D0*F4     +3.D0*F3-F2) * AUX315
      A36(IELEM) =  (-2.D0*F4     +3.D0*F3-F1) * AUX315
C
      A45(IELEM) =  2.D0*(6.D0*(F4+F5)+4.D0*F6-F1-F3) * AUX315
      A46(IELEM) =  2.D0*(6.D0*(F4+F6)+4.D0*F5-F2-F3) * AUX315
C      
      A56(IELEM) =  2.D0*(6.D0*(F6+F5)+4.D0*F4-F2-F1) * AUX315
C
5     CONTINUE
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'MT06CC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT06CC (BIEF) :',/,
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
