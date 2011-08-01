C                       *****************
                        SUBROUTINE KSPG11
C                       *****************
C
     *(KX,KY,XEL,YEL,U,V,IKLE,NELEM,NELMAX,XMUL)
C
C***********************************************************************
C BIEF VERSION 5.1           08/12/94    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : CALCUL D'UN VECTEUR QUI SERA UTILISE
C             PAR LA METHODE :
C
C             STREAMLINE UPWIND PETROV GALERKIN
C
C             AVEC UN DECENTREMENT EGAL A UN
C
C                    DX   U
C             KX = -----------
C                  2 NORME(U)
C
C                    DY   V
C             KY = -----------
C                  2 NORME(U)
C
C-----------------------------------------------------------------------
C
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  12 : TRIANGLE P2            6
C  21 : QUADRILATERE Q1        4                       OUI
C  41 : PRISMES TELEMAC-3D     6
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      KX,KY     | -->|  COORDONNEES DU VECTEUR UNITAIRE.
C |      XEL,YEL,..| -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |      U,V,W     | -->|  COMPOSANTES DE LA VITESSE.
C |      IKLE      | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE.
C |      NELEM     | -->|  NOMBRE D'ELEMENTS DU MAILLAGE.
C |      NELMAX    | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE.
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |      XMUL      | -->|  COEFICIENT MULTIPLICATEUR
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C  PROGRAMMES APPELES :
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NELMAX
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,*)
      DOUBLE PRECISION, INTENT(INOUT) :: KX(NELEM),KY(NELEM)
      DOUBLE PRECISION, INTENT(IN)    :: U(*),V(*),XMUL
      DOUBLE PRECISION, INTENT(IN)    :: XEL(NELMAX,*),YEL(NELMAX,*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,I1,I2,I3
C      
      DOUBLE PRECISION UMOY,VMOY,H,SUNORM,X2,X3,Y2,Y3
      DOUBLE PRECISION SURFAC,GP1X,GP1Y,GP2X,GP2Y,GP3X,GP3Y
      DOUBLE PRECISION A1,A2,A3,H1,H2,H3,C1,C2,C3,UNORM,VNORM
C
      INTRINSIC MAX,SQRT
C
C-----------------------------------------------------------------------
C
      DO 10 IELEM = 1 , NELEM
C
        I1 = IKLE(IELEM,1)
        I2 = IKLE(IELEM,2)
        I3 = IKLE(IELEM,3)
C
        X2 = XEL(IELEM,2)
        X3 = XEL(IELEM,3)
        Y2 = YEL(IELEM,2)
        Y3 = YEL(IELEM,3)
C
        GP1X = Y2-Y3
        GP1Y = X3-X2
        SUNORM = 1.D0 / SQRT(GP1X**2+GP1Y**2)
        GP1X = GP1X * SUNORM
        GP1Y = GP1Y * SUNORM
C
        GP2X = Y3
        GP2Y =   -X3
        SUNORM = 1.D0 / SQRT(GP2X**2+GP2Y**2)
        GP2X = GP2X * SUNORM
        GP2Y = GP2Y * SUNORM
C
        GP3X =   -Y2
        GP3Y = X2
        SUNORM = 1.D0 / SQRT(GP3X**2+GP3Y**2)
        GP3X = GP3X * SUNORM
        GP3Y = GP3Y * SUNORM
C
        C3 = SQRT(  X2**2     +  Y2**2 )
        C1 = SQRT( (X3-X2)**2 + (Y3-Y2)**2 )
        C2 = SQRT(  X3**2     +  Y3**2 )
C
        SURFAC = 0.5D0 * (X2*Y3 - X3*Y2)
C
        H1 = 2*SURFAC/C1
        H2 = 2*SURFAC/C2
        H3 = 2*SURFAC/C3
C
        H = MAX(H1,H2,H3)
C
        UMOY = U(I1) + U(I2) + U(I3)
        VMOY = V(I1) + V(I2) + V(I3)
C
        SUNORM = 1.D0 / MAX ( SQRT(UMOY**2+VMOY**2) , 1.D-10 )
C
        UNORM = UMOY * SUNORM
        VNORM = VMOY * SUNORM
C
        A1 = GP1X * UNORM + GP1Y * VNORM
        A2 = GP2X * UNORM + GP2Y * VNORM
        A3 = GP3X * UNORM + GP3Y * VNORM
C
        IF(A1*H.GT.H1) H = H1
        IF(A2*H.GT.H2) H = H2
        IF(A3*H.GT.H3) H = H3
C
        KX(IELEM) = 0.33333333D0 * XMUL * H * UNORM
        KY(IELEM) = 0.33333333D0 * XMUL * H * VNORM
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
