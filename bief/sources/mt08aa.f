C                       *****************
                        SUBROUTINE MT08AA
C                       *****************
C
     *(  A11 , A12 , A13 ,
     *   A21 , A22 , A23 ,
     *   A31 , A32 , A33 ,
     *   XMUL,SF,F,XEL,YEL,IKLE1,IKLE2,IKLE3,NELEM,NELMAX,ICOORD)
C
C ---------------------------------
C CHANGEMENT :
C    LA MATRICE EST CALCULEE AUSSI
C    SI ON EST EN QUASI-BULLE
C    (C'EST LA MEME)
C ---------------------------------
C
C***********************************************************************
C BIEF VERSION 5.1           09/12/94    J-M HERVOUET (LNH) 30 87 80 18
C                                        C  MOULIN    (LNH) 30 87 83 81
C***********************************************************************
C
C  FONCTION : CALCUL DES COEFFICIENTS DE LA MATRICE SUIVANTE
C
C  EXEMPLE AVEC ICOORD = 1
C
C                 /                     D
C A(I,J)=-XMUL   /  PSI2(J) *    F    * --( PSI1(I) ) D(OMEGA)
C               /OMEGA                  DX
C
C  ATTENTION AU SIGNE MOINS ||
C
C  AVEC ICOORD=2 ON AURAIT UNE DERIVEE SUIVANT Y
C
C  PSI1 : BASES DE TYPE TRIANGLE P1
C  PSI2 : BASES DE TYPE IELM2
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
      USE BIEF, EX_MT08AA => MT08AA
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
      DOUBLE PRECISION, INTENT(IN) :: XEL(NELMAX,3),YEL(NELMAX,3)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IELMF
      DOUBLE PRECISION SUR24,X2,X3,Y2,Y3,F1,F2,F3,F123
C
C-----------------------------------------------------------------------
C      
      SUR24 = XMUL/24.D0
C
C-----------------------------------------------------------------------
C
      IELMF=SF%ELM
C
C  CAS OU F EST DE DISCRETISATION P1
C
      IF((IELMF.EQ.11).OR.(IELMF.EQ.12)) THEN
C
C TH
C MEME MATRICE SI F QUASI-BULLE
C TH
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
      F123 = F1 + F2 + F3
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (Y3-Y2) * (  F123 + F2  )
      A13(IELEM) = (Y3-Y2) * (  F123 + F3  )
      A23(IELEM) =  Y3     * ( -F123 - F3  )
      A21(IELEM) =  Y3     * ( -F123 - F1  )
      A31(IELEM) =     Y2  * (  F123 + F1  )
      A32(IELEM) =     Y2  * (  F123 + F2  )
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
      F1  =  F(IKLE1(IELEM)) * SUR24
      F2  =  F(IKLE2(IELEM)) * SUR24
      F3  =  F(IKLE3(IELEM)) * SUR24
      F123 = F1 + F2 + F3
C
C   TERMES EXTRADIAGONAUX
C
      A12(IELEM) = (X2-X3) * (  F123 + F2  )
      A13(IELEM) = (X2-X3) * (  F123 + F3  )
      A23(IELEM) =     X3  * (  F123 + F3  )
      A21(IELEM) =     X3  * (  F123 + F1  )
      A31(IELEM) =  X2     * ( -F123 - F1  )
      A32(IELEM) =  X2     * ( -F123 - F2  )
C
C   TERMES DIAGONAUX
C
      A11(IELEM) = -A21(IELEM) -A31(IELEM)
      A22(IELEM) = -A12(IELEM) -A32(IELEM)
      A33(IELEM) = -A13(IELEM) -A23(IELEM)
C
2     CONTINUE
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,200) ICOORD
          IF (LNG.EQ.2) WRITE(LU,201) ICOORD
200       FORMAT(1X,'MT08AA (BIEF) : COMPOSANTE IMPOSSIBLE ',
     *              1I6,' VERIFIER ICOORD')
201       FORMAT(1X,'MT08AA (BIEF) : IMPOSSIBLE COMPONENT ',
     *              1I6,' CHECK ICOORD')
          CALL PLANTE(1)
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
100    FORMAT(1X,'MT08AA (BIEF) :',/,
     *        1X,'DISCRETISATION DE F : ',1I6,' NON PREVUE')
101    FORMAT(1X,'MT08AA (BIEF) :',/,
     *        1X,'DISCRETIZATION OF F : ',1I6,' NOT AVAILABLE')
       CALL PLANTE(1)
       STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END  
 
 
