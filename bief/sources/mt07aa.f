C                       *****************
                        SUBROUTINE MT07AA
C                       *****************
C
     *( A11 , A12 , A13 ,
     *        A22 , A23 ,
     *              A33 ,
     *  XMUL,SF,F,SURFAC,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.1           30/06/93    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CONSTRUCTION D'UNE MATRICE DE MASSE
C
C            AVEC MASS-LUMPING LOCAL EN FONCTION D'UN COEFFICIENT
C            LOCAL TETA (FONCTION P0) (ICI LA FONCTION F)
C
C            M = TETA     * MATRICE DE MASSE
C              + (1-TETA) * MATRICE DIAGONALE
C
C            CETTE MATRICE EST ENSUITE MULTIPLIEE PAR XMUL
C
C
C     L'ELEMENT EST LE TRIANGLE P1
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
      USE BIEF, EX_MT07AA => MT07AA
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
      DOUBLE PRECISION, INTENT(INOUT) ::        A22(*),A23(*)
      DOUBLE PRECISION, INTENT(INOUT) ::               A33(*)
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
      INTEGER IELMF,IELEM
C
      DOUBLE PRECISION SUR12,DET,T
C
C-----------------------------------------------------------------------
C
      SUR12 = XMUL/12.D0
C
C-----------------------------------------------------------------------
C
      IELMF = SF%ELM
C
C  CAS OU F EST CONSTANTE PAR ELEMENT
C
      IF(IELMF.EQ.10) THEN
C
      DO 1 IELEM = 1 , NELEM
C
      DET = SURFAC(IELEM) * SUR12
      T   = F(IELEM)
C
C***********************************************************************
C
C  ELEMENTS EXTERIEURS A LA DIAGONALE
C
      A12(IELEM) = T * DET
      A13(IELEM) = T * DET
      A23(IELEM) = T * DET
C
C  TERMES DIAGONAUX
C
      A11(IELEM) = ( DET + DET ) * (2.D0 - T)
      A22(IELEM) = ( DET + DET ) * (2.D0 - T)
      A33(IELEM) = ( DET + DET ) * (2.D0 - T)
C
C   FIN DE LA BOUCLE SUR LES ELEMENTS
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
C     AUTRES TYPES DE DISCRETISATION DE F
C
      ELSE
C
       IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
       IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100    FORMAT(1X,'MT07AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101    FORMAT(1X,'MT07AA (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F NOT AVAILABLE:',1I6,
     *        1X,'REAL NAME: ',A6)
       CALL PLANTE(0)
       STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
