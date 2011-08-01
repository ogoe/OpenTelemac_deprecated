C                       *****************
                        SUBROUTINE MT06BB
C                       *****************
C
     *( A11 , A12 , A13 , A14 ,
     *        A22 , A23 , A24 ,
     *              A33 , A34 ,
     *                    A44 ,
     *  XMUL,SF,F,SURFAC,IKLE1,IKLE2,IKLE3,IKLE4,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95     J-M HERVOUET (LNH) 30 87 80 18
C                                           C   MOULIN (LNH) 30 87 83 81
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE SOMME(F*PSII*PSIJ)
C
C            AVEC :  P1 QUASI-BULLE
C                    P2 QUASI-BULLE
C                    F  P1 OU QUASI-BULLE
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
C |     IKLE1..3   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
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
      USE BIEF, EX_MT06BB => MT06BB
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
      DOUBLE PRECISION, INTENT(IN) :: SURFAC(NELMAX)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      DOUBLE PRECISION F1,F2,F3,F4,XMS090,XMS180,XMS540
      DOUBLE PRECISION XMS018,XMS054,XMS006,XMS009,XMS036
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
      XMS009 = XMUL /  9.D0
      XMS006 = XMUL /  6.D0
      XMS018 = XMUL / 18.D0
      XMS036 = XMUL / 36.D0
C
      DO 1 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
         F1  =  F(IELEM) * SURFAC(IELEM)
C
C  TERMES DIAGONAUX
C
         A11(IELEM) = F1*XMS009
         A22(IELEM) = F1*XMS009
         A33(IELEM) = F1*XMS009
         A44(IELEM) = F1*XMS006
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = F1*XMS036
         A13(IELEM) = F1*XMS036
         A14(IELEM) = F1*XMS018
         A23(IELEM) = F1*XMS036
         A24(IELEM) = F1*XMS018
         A34(IELEM) = F1*XMS018
C
1     CONTINUE
C
C
C-----------------------------------------------------------------------
C
C  CAS OU F EST LINEAIRE
C
      ELSEIF(IELMF.EQ.11) THEN
C
      XMS054 = XMUL /  54.D0
      XMS018 = XMUL /  18.D0
      XMS540 = XMUL / 540.D0
C
      DO 4 IELEM = 1 , NELEM
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
         A11(IELEM) = (SURFAC(IELEM)*(F3+F2+4*F1))*XMS054
         A22(IELEM) = (SURFAC(IELEM)*(F3+4*F2+F1))*XMS054
         A33(IELEM) = (SURFAC(IELEM)*(4*F3+F2+F1))*XMS054
         A44(IELEM) = (SURFAC(IELEM)*(F3+F2+F1))  *XMS018
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = (SURFAC(IELEM)*(   F3+ 7*F2+ 7*F1))*XMS540
         A13(IELEM) = (SURFAC(IELEM)*( 7*F3+   F2+ 7*F1))*XMS540
         A14(IELEM) = (SURFAC(IELEM)*( 7*F3+ 7*F2+16*F1))*XMS540
         A23(IELEM) = (SURFAC(IELEM)*( 7*F3+ 7*F2+   F1))*XMS540
         A24(IELEM) = (SURFAC(IELEM)*( 7*F3+16*F2+ 7*F1))*XMS540
         A34(IELEM) = (SURFAC(IELEM)*(16*F3+ 7*F2+ 7*F1))*XMS540
C
4     CONTINUE
C
C-----------------------------------------------------------------------
C
      ELSEIF(IELMF.EQ.12) THEN
C
C-----------------------------------------------------------------------
C
C   DISCRETISATION QUASI-BULLE DE F :
C
C
      XMS090 = XMUL / 90.D0
      XMS180 = XMUL / 180.D0
C
      DO 5 IELEM = 1 , NELEM
C
C   INITIALISATION DES VARIABLES GEOMETRIQUES
C
         F1  =  F(IKLE1(IELEM))
         F2  =  F(IKLE2(IELEM))
         F3  =  F(IKLE3(IELEM))
         F4  =  F(IKLE4(IELEM))
C
C  TERMES DIAGONAUX
C
         A11(IELEM) = (SURFAC(IELEM)*(  F3+2*F4+  F2+6*F1))*XMS090
         A22(IELEM) = (SURFAC(IELEM)*(  F3+2*F4+6*F2+  F1))*XMS090
         A33(IELEM) = (SURFAC(IELEM)*(6*F3+2*F4+  F2+  F1))*XMS090
         A44(IELEM) = (SURFAC(IELEM)*(2*F3+9*F4+2*F2+2*F1))*XMS090
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = (SURFAC(IELEM)*(F4+2*F2+2*F1))*XMS180
         A13(IELEM) = (SURFAC(IELEM)*(2*F3+F4+2*F1))*XMS180
         A14(IELEM) = (SURFAC(IELEM)*(F3+4*F4+F2+4*F1))*XMS180
         A23(IELEM) = (SURFAC(IELEM)*(2*F3+F4+2*F2))*XMS180
         A24(IELEM) = (SURFAC(IELEM)*(F3+4*F4+4*F2+F1))*XMS180
         A34(IELEM) = (SURFAC(IELEM)*(4*F3+4*F4+F2+F1))*XMS180
C
5     CONTINUE
C
C-----------------------------------------------------------------------
C
C   AUTRE DISCRETISATION
C      ELSEIF(IELMF.EQ.XX) THEN
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'MT06BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT06BB (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *        1X,'REAL NAME: ',A6)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
