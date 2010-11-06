C                       *****************
                        SUBROUTINE MT07BB
C                       *****************
C
     *( A11 , A12 , A13 , A14 ,
     *        A22 , A23 , A24 ,
     *              A33 , A34 ,
     *                    A44 ,
     *  XMUL,SF,F,SURFAC,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           10/01/95     J-M HERVOUET (LNH) 30 87 80 18
C                                           C   MOULIN (LNH) 30 87 83 81
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE : T*M + (1-T)*DM
C            POUR LES TRIANGLES QUASI-BULLE.
C
C     M  EST LA MATRICE DE MASSE QB*QB
C     DM EST M MASSE-LUMPEE
C     T EST UN VECTEUR CONSTANT PAR ELEMENT
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
      USE BIEF, EX_MT07BB => MT07BB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
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
      DOUBLE PRECISION T,XSUR06,XSUR09,XSUR36
      INTEGER IELEM
C
C=======================================================================
C
      XSUR06 = XMUL/6.D0
      XSUR09 = XMUL/9.D0
      XSUR36 = XMUL/36.D0
C
C-----------------------------------------------------------------------
C
      IF(SF%ELM.EQ.10) THEN
C
C-----------------------------------------------------------------------
C
C   DISCRETISATION P0 DE F :
C
      DO 5 IELEM = 1 , NELEM
C
C   INITIALISATION DES INTERMEDIAIRES DE CALCULS
C
         T = F(IELEM)
C
C  TERMES DIAGONAUX
C
         A11(IELEM) = -(T-2)*SURFAC(IELEM)*XSUR09
         A22(IELEM) = A11(IELEM)
         A33(IELEM) = A11(IELEM)
         A44(IELEM) = -(T-2)*SURFAC(IELEM)*XSUR06
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = SURFAC(IELEM)*T*XSUR36
         A13(IELEM) = A12(IELEM)
         A14(IELEM) = 2*A12(IELEM)
         A23(IELEM) = A12(IELEM)
         A24(IELEM) = A14(IELEM)
         A34(IELEM) = A14(IELEM)
C
5     CONTINUE
C
C-----------------------------------------------------------------------
C
C  AUTRE DISCRETISATION DE F
C      ELSEIF(SF%ELM.EQ.XX) THEN
C
C-----------------------------------------------------------------------
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) SF%ELM,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) SF%ELM,SF%NAME
100    FORMAT(1X,'MT07BB (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT07BB (BIEF) :',/,
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
 
 
