C                       *****************
                        SUBROUTINE MT07CC
C                       *****************
C
     *( A11 , A12 , A13 , A14 , A15 , A16 ,
     *        A22 , A23 , A24 , A25 , A26 ,
     *              A33 , A34 , A35 , A36 ,
     *                    A44 , A45 , A46 ,
     *                          A55 , A56 ,
     *                                A66 ,
     *  XMUL,SF,F,SURFAC,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.9       01/08/08     A FROEHLY (MATMECA) 01 30 87 80 18
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE : T*M + (1-T)*DM
C            POUR LES TRIANGLES P2.
C
C     M  EST LA MATRICE DE MASSE P2*P2
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
      USE BIEF!, EX_MT07CC => MT07CC
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
C
      DOUBLE PRECISION, INTENT(INOUT) :: A11(*),A12(*),A13(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A14(*),A15(*),A16(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A22(*),A23(*),A24(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A25(*),A26(*),A33(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A34(*),A35(*),A36(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A44(*),A45(*),A46(*)
      DOUBLE PRECISION, INTENT(INOUT) :: A55(*),A56(*),A66(*)
     
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
      DOUBLE PRECISION T,XSUR30,XSUR45,XSUR180
      INTEGER IELEM
C
C=======================================================================
C
      XSUR30 = XMUL/30.D0
      XSUR45 = XMUL/45.D0
      XSUR180 = XMUL/180.D0
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
C  POUR LES ELEMENTS P2, DM(1) = DM(2) = DM(3) = 0
C                     ET DM(4) = DM(5) = DM(6) = S/3

         A11(IELEM) = T*SURFAC(IELEM)*XSUR30
         A22(IELEM) = A11(IELEM)
         A33(IELEM) = A11(IELEM)
         A44(IELEM) = (15.D0-7.D0*T)*SURFAC(IELEM)*XSUR45
         A55(IELEM) = A44(IELEM)
         A66(IELEM) = A44(IELEM)
C
C  TERMES EXTRADIAGONAUX
C
         A12(IELEM) = -SURFAC(IELEM)*T*XSUR180
         A13(IELEM) = A12(IELEM)
         A14(IELEM) = 0.D0
         A15(IELEM) = -SURFAC(IELEM)*T*XSUR45
         A16(IELEM) = 0.D0
         A23(IELEM) = A12(IELEM)
         A24(IELEM) = 0.D0
         A25(IELEM) = 0.D0
         A26(IELEM) = A15(IELEM)
         A34(IELEM) = A15(IELEM)
         A35(IELEM) = 0.D0
         A36(IELEM) = 0.D0      
         A45(IELEM) = (15.D0-11.D0*T)*SURFAC(IELEM)*XSUR45
         A46(IELEM) = A45(IELEM) 
         A56(IELEM) = A45(IELEM) 
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
100    FORMAT(1X,'MT07CC (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT07CC (BIEF) :',/,
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
