C                       *****************
                        SUBROUTINE MT04TT
C                       *****************
C
     *( T,XM,XMUL,SU,SV,SW,U,V,W,X,Y,Z,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.6        22/08/05    J-M HERVOUET (LNHE) 01 30 87 80 18
C                                      F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C    FONCTION : CALCUL D'UNE MATRICE POUR LA METHODE SUPG
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                                 / ->             ->
C                    A    = XMUL /  U . GRAD(P ) * U . GRAD(P ) * J(X,Y) DXDY
C                     I J       /S            I              J
C
C    PAR MAILLE ELEMENTAIRE.
C
C
C    ATTENTION !!!!!!!!!!!!!!!!!!!!!
C
C    ICI LA COMPOSANTE VERTICALE EST NEGLIGEE, CONFORMEMENT A LA NOTE DE PRINCIPE
C    ET CONTRAIREMENT AUX PRISMES.
C
C    ATTENTION !!!!!!!!!!!!!!!!!!!!!
C
C
C
C    J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C    L'ELEMENT EST LE TETRAEDRE AVEC INTERPOLATION LINEAIRE
C
C-----------------------------------------------------------------------
C
C     CONVENTION POUR LE STOCKAGE DES TERMES EXTRA-DIAGONAUX :
C
C     XM(IELEM, 1)  ---->  M(1,2) = M(2,1)
C     XM(IELEM, 2)  ---->  M(1,3) = M(3,1)
C     XM(IELEM, 3)  ---->  M(1,4) = M(4,1)
C     XM(IELEM, 4)  ---->  M(2,3) = M(3,2)
C     XM(IELEM, 5)  ---->  M(2,4) = M(4,2)
C     XM(IELEM, 6)  ---->  M(3,4) = M(4,3)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     A11,A12... |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SU,SV,SW   | -->|  STRUCTURES DE U,V ET W.
C |     U,V,W      | -->|  FONCTIONS INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     X,Y,Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     IKLE1..4   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
      USE BIEF, EX_MT04TT => MT04TT
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,4)
C
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,4),XM(NELMAX,6)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*),W(*)
C
C     STRUCTURES DE U,V,W 
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SU,SV,SW 
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES 
C     
      DOUBLE PRECISION X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
      DOUBLE PRECISION U1,U2,U3,U4,V1,V2,V3,V4,W1,W2,W3,W4,XJAC       
      DOUBLE PRECISION COEF,SUR120,AUX1,AUX2,AUX4,AUX5,AUX6,AUX7
      DOUBLE PRECISION A1,A2,A3                   
      INTEGER I1,I2,I3,I4,IELEM   
C
C***********************************************************************
C
      SUR120=XMUL/120.D0
C
      IF((SU%ELM.EQ.31.AND.SV%ELM.EQ.31.AND.SW%ELM.EQ.31).OR.
     *   (SU%ELM.EQ.51.AND.SV%ELM.EQ.51.AND.SW%ELM.EQ.51)     ) THEN
C
C-----------------------------------------------------------------------
C
C   LINEAR DISCRETISATION OF DIFFUSION COEFFICIENTS
C
C   BOUCLE SUR LES TETRAEDRES
C
      DO 20 IELEM=1,NELEM

      I1=IKLE(IELEM,1)
      I2=IKLE(IELEM,2)
      I3=IKLE(IELEM,3)
      I4=IKLE(IELEM,4)
C   
C-----------------------------------------------------------------------
C     
      X2=X(I2)-X(I1)
      Y2=Y(I2)-Y(I1)
      Z2=Z(I2)-Z(I1)
      X3=X(I3)-X(I1)
      Y3=Y(I3)-Y(I1)
      Z3=Z(I3)-Z(I1)
      X4=X(I4)-X(I1)
      Y4=Y(I4)-Y(I1)
      Z4=Z(I4)-Z(I1)
C
C-----------------------------------------------------------------------
C
      U1 = U(IKLE(IELEM,1))
      U2 = U(IKLE(IELEM,2))
      U3 = U(IKLE(IELEM,3))
      U4 = U(IKLE(IELEM,4))
      V1 = V(IKLE(IELEM,1))
      V2 = V(IKLE(IELEM,2))
      V3 = V(IKLE(IELEM,3))
      V4 = V(IKLE(IELEM,4))
      W1 = W(IKLE(IELEM,1))
      W2 = W(IKLE(IELEM,2))
      W3 = W(IKLE(IELEM,3))
      W4 = W(IKLE(IELEM,4))    
C
      XJAC = X2*Y3*Z4-X2*Y4*Z3-Y2*X3*Z4+Y2*X4*Z3+Z2*X3*Y4-Z2*X4*Y3
C
      COEF = SUR120/XJAC
C
      A1 = -Y3*Z4+Y4*Z3+Y2*Z4-Z2*Y4-Y2*Z3+Z2*Y3
      A2 = -X3*Z4+X4*Z3+X2*Z4-Z2*X4-X2*Z3+Z2*X3
      A3 = -X3*Y4+X4*Y3+X2*Y4-Y2*X4-X2*Y3+Y2*X3
C
      AUX1 = U3*U2+U1*U3+U1*U2+U4**2+U1*U4+U4*U2+U3**2+U2**2+U1**2+U3*U4
C
      AUX2 = V1*U4+2*V1*U1+V1*U2+V3*U1+V3*U2+V1*U3+V2*U1+V2*U3
     #      +V4*U1+2*V2*U2+2*V3*U3+V4*U2+V2*U4+V4*U3+V3*U4+2*V4*U4
C
      AUX4 = V1**2+V1*V3+V3**2+V2**2+V3*V2+V1*V2+V4**2+V1*V4+V3*V4+V2*V4
      AUX5 = V3*W2+V3*W1+V1*W2+2*V1*W1+V2*W3+V2*W1+V1*W3+2*V3*W3+2*V2*W2
     #+V4*W1+V4*W3+V4*W2+2*V4*W4+V3*W4+V2*W4+V1*W4
      AUX6 = W2**2+W3**2+W2*W3+W1**2+W1*W3+W4**2+W1*W2+W1*W4+W4*W2+W4*W3
      AUX7 = 2*U2*W2+2*U3*W3+U3*W2+U3*W1+U1*W2+2*U1*W1+U2*W3+U2*W1+
     #U1*W3+U4*W2+U4*W1+2*U4*W4+U1*W4+U2*W4+U4*W3+U3*W4

C
      T(IELEM,1)= (A1**2*AUX1-
     #A1*
     #A2*AUX2+
     #A1*
     #A3*AUX7+
     #A2**2*AUX4-
     #A2*
     #A3*AUX5+
     #A3**2*AUX6)*2*COEF
C
      XM(IELEM,1) = (2*(Y3*Z4-Y4*Z3)*
     #A1*AUX1-
     #(Y3*Z4-Y4*Z3)*A2*AUX2+
     #(Y3*Z4-Y4*Z3)*A3*AUX7-
     #(X3*Z4-X4*Z3)*A1*AUX2+
     #2*(X3*Z4-X4*Z3)*A2*AUX4-
     #(X3*Z4-X4*Z3)*A3*AUX5-
     #(-X3*Y4+X4*Y3)*A1*AUX7+
     #(-X3*Y4+X4*Y3)*A2*AUX5-
     #2*(-X3*Y4+X4*Y3)*A3*AUX6)*COEF
C
      XM(IELEM,2) = (2*(-Y2*Z4+Z2*Y4)*
     #A1*AUX1-
     #(-Y2*Z4+Z2*Y4)*A2*AUX2+
     #(-Y2*Z4+Z2*Y4)*A3*AUX7-
     #(-X2*Z4+Z2*X4)*A1*AUX2+
     #2*(-X2*Z4+Z2*X4)*A2*AUX4-
     #(-X2*Z4+Z2*X4)*A3*AUX5-
     #(X2*Y4-Y2*X4)*A1*AUX7+
     #(X2*Y4-Y2*X4)*A2*AUX5-
     #2*(X2*Y4-Y2*X4)*A3*AUX6)*COEF
C
      XM(IELEM,3) = (-2*(-Y2*Z3+Z2*Y3)*
     #A1*AUX1+
     #(-Y2*Z3+Z2*Y3)*A2*AUX2-
     #(-Y2*Z3+Z2*Y3)*A3*AUX7+
     #(-X2*Z3+Z2*X3)*A1*AUX2-
     #2*(-X2*Z3+Z2*X3)*A2*AUX4+
     #(-X2*Z3+Z2*X3)*A3*AUX5+
     #(X2*Y3-Y2*X3)*A1*AUX7-
     #(X2*Y3-Y2*X3)*A2*AUX5+
     #2*(X2*Y3-Y2*X3)*A3*AUX6)*COEF
C
      T(IELEM,2) =((Y3*Z4-Y4*Z3)**2*AUX1-
     #(Y3*Z4-Y4*Z3)*(X3*Z4-X4*Z3)*AUX2
     #-(Y3*Z4-Y4*Z3)*(-X3*Y4+X4*Y3)*AUX7+
     #(X3*Z4-X4*Z3)**2*AUX4+
     #(X3*Z4-X4*Z3)*(-X3*Y4+X4*Y3)*AUX5+
     #(-X3*Y4+X4*Y3)**2*AUX6)*2*COEF
C 
      XM(IELEM,4) =  (2*(-Y2*Z4+Z2*Y4)*(Y3*Z4-Y4*Z3)*AUX1-
     #(-Y2*Z4+Z2*Y4)*(X3*Z4-X4*Z3)*AUX2-
     #(-Y2*Z4+Z2*Y4)*(-X3*Y4+X4*Y3)*AUX7
     #-(Y3*Z4-Y4*Z3)*(-X2*Z4+Z2*X4)*AUX2+
     #2*(-X2*Z4+Z2*X4)*(X3*Z4-X4*Z3)*AUX4+
     #(-X2*Z4+Z2*X4)*(-X3*Y4+X4*Y3)*AUX5-
     #(Y3*Z4-Y4*Z3)*(X2*Y4-Y2*X4)*AUX7+
     #(X3*Z4-X4*Z3)*(X2*Y4-Y2*X4)*AUX5+
     #2*(X2*Y4-Y2*X4)*(-X3*Y4+X4*Y3)*AUX6)*COEF
C
      XM(IELEM,5) = (-2*(-Y2*Z3+Z2*Y3)*(Y3*Z4-Y4*Z3)*AUX1+
     #(-Y2*Z3+Z2*Y3)*(X3*Z4-X4*Z3)*AUX2+
     #(-Y2*Z3+Z2*Y3)*(-X3*Y4+X4*Y3)*AUX7+
     #(Y3*Z4-Y4*Z3)*(-X2*Z3+Z2*X3)*AUX2-
     #2*(-X2*Z3+Z2*X3)*(X3*Z4-X4*Z3)*AUX4-
     #(-X2*Z3+Z2*X3)*(-X3*Y4+X4*Y3)*AUX5+
     #(Y3*Z4-Y4*Z3)*(X2*Y3-Y2*X3)*AUX7-
     #(X3*Z4-X4*Z3)*(X2*Y3-Y2*X3)*AUX5-
     #2*(X2*Y3-Y2*X3)*(-X3*Y4+X4*Y3)*AUX6)*COEF
C
      T(IELEM,3) = ((-Y2*Z4+Z2*Y4)**2*AUX1-
     #(-Y2*Z4+Z2*Y4)*(-X2*Z4+Z2*X4)*AUX2-
     #(-Y2*Z4+Z2*Y4)*(X2*Y4-Y2*X4)*AUX7+
     #(-X2*Z4+Z2*X4)**2*AUX4+(-X2*Z4+Z2*X4)*(X2*Y4-Y2*X4)*AUX5+
     #(X2*Y4-Y2*X4)**2*AUX6)*2*COEF
C
      XM(IELEM,6) = (-2*(-Y2*Z3+Z2*Y3)*(-Y2*Z4+Z2*Y4)*AUX1+
     #(-Y2*Z3+Z2*Y3)*(-X2*Z4+Z2*X4)*AUX2+
     #(-Y2*Z3+Z2*Y3)*(X2*Y4-Y2*X4)*AUX7+
     #(-Y2*Z4+Z2*Y4)*(-X2*Z3+Z2*X3)*AUX2-
     #2*(-X2*Z3+Z2*X3)*(-X2*Z4+Z2*X4)*AUX4-
     #(-X2*Z3+Z2*X3)*(X2*Y4-Y2*X4)*AUX5+
     #(-Y2*Z4+Z2*Y4)*(X2*Y3-Y2*X3)*AUX7-
     #(-X2*Z4+Z2*X4)*(X2*Y3-Y2*X3)*AUX5-
     #2*(X2*Y3-Y2*X3)*(X2*Y4-Y2*X4)*AUX6)*COEF
C
      T(IELEM,4) = ((-Y2*Z3+Z2*Y3)**2*AUX1-
     #(-Y2*Z3+Z2*Y3)*(-X2*Z3+Z2*X3)*AUX2-
     #(-Y2*Z3+Z2*Y3)*(X2*Y3-Y2*X3)*AUX7
     #+(-X2*Z3+Z2*X3)**2*AUX4+(-X2*Z3+Z2*X3)*(X2*Y3-Y2*X3)*AUX5+
     #(X2*Y3-Y2*X3)**2*AUX6)*2*COEF
C
C 
C-----------------------------------------------------------------------
C        
 20   CONTINUE
C
C---------------------------------------------------------------
C
      ELSE IF((SU%ELM.EQ.30.AND.SV%ELM.EQ.30.AND.SW%ELM.EQ.30).OR.
     *       (SU%ELM.EQ.50.AND.SV%ELM.EQ.50.AND.SW%ELM.EQ.50)     ) THEN

C   P0 DISCRETISATION OF DIFFUSION COEFFICIENTS
C
C   BOUCLE SUR LES TETRAEDRES
C
      DO 21 IELEM=1,NELEM
C
      I1=IKLE(IELEM,1)
      I2=IKLE(IELEM,2)
      I3=IKLE(IELEM,3)
      I4=IKLE(IELEM,4)
C   
C-----------------------------------------------------------------------
C     
      X2=X(I2)-X(I1)
      Y2=Y(I2)-Y(I1)
      Z2=Z(I2)-Z(I1)
      X3=X(I3)-X(I1)
      Y3=Y(I3)-Y(I1)
      Z3=Z(I3)-Z(I1)
      X4=X(I4)-X(I1)
      Y4=Y(I4)-Y(I1)
      Z4=Z(I4)-Z(I1)
C
C-----------------------------------------------------------------------
C
      U1 = U(IELEM)
      U2 = U1
      U3 = U1
      U4 = U1
      V1 = V(IELEM)
      V2 = V1
      V3 = V1
      V4 = V1
      W1 = W(IELEM)
      W2 = W1
      W3 = W1
      W4 = W1
C
C
      XJAC = X2*Y3*Z4-X2*Y4*Z3-Y2*X3*Z4+Y2*X4*Z3+Z2*X3*Y4-Z2*X4*Y3
C
      COEF = SUR120/XJAC
C
      A1 = -Y3*Z4+Y4*Z3+Y2*Z4-Z2*Y4-Y2*Z3+Z2*Y3
      A2 = -X3*Z4+X4*Z3+X2*Z4-Z2*X4-X2*Z3+Z2*X3
      A3 = -X3*Y4+X4*Y3+X2*Y4-Y2*X4-X2*Y3+Y2*X3
C
      AUX1 = U3*U2+U1*U3+U1*U2+U4**2+U1*U4+U4*U2+U3**2+U2**2+U1**2+U3*U4
C
      AUX2 = V1*U4+2*V1*U1+V1*U2+V3*U1+V3*U2+V1*U3+V2*U1+V2*U3
     #      +V4*U1+2*V2*U2+2*V3*U3+V4*U2+V2*U4+V4*U3+V3*U4+2*V4*U4
C
      AUX4 = V1**2+V1*V3+V3**2+V2**2+V3*V2+V1*V2+V4**2+V1*V4+V3*V4+V2*V4
      AUX5 = V3*W2+V3*W1+V1*W2+2*V1*W1+V2*W3+V2*W1+V1*W3+2*V3*W3+2*V2*W2
     #+V4*W1+V4*W3+V4*W2+2*V4*W4+V3*W4+V2*W4+V1*W4
      AUX6 = W2**2+W3**2+W2*W3+W1**2+W1*W3+W4**2+W1*W2+W1*W4+W4*W2+W4*W3
      AUX7 = 2*U2*W2+2*U3*W3+U3*W2+U3*W1+U1*W2+2*U1*W1+U2*W3+U2*W1+
     #U1*W3+U4*W2+U4*W1+2*U4*W4+U1*W4+U2*W4+U4*W3+U3*W4

C
      T(IELEM,1)= (A1**2*AUX1-
     #A1*
     #A2*AUX2+
     #A1*
     #A3*AUX7+
     #A2**2*AUX4-
     #A2*
     #A3*AUX5+
     #A3**2*AUX6)*2*COEF
C
      XM(IELEM,1) = (2*(Y3*Z4-Y4*Z3)*
     #A1*AUX1-
     #(Y3*Z4-Y4*Z3)*A2*AUX2+
     #(Y3*Z4-Y4*Z3)*A3*AUX7-
     #(X3*Z4-X4*Z3)*A1*AUX2+
     #2*(X3*Z4-X4*Z3)*A2*AUX4-
     #(X3*Z4-X4*Z3)*A3*AUX5-
     #(-X3*Y4+X4*Y3)*A1*AUX7+
     #(-X3*Y4+X4*Y3)*A2*AUX5-
     #2*(-X3*Y4+X4*Y3)*A3*AUX6)*COEF
C
      XM(IELEM,2) = (2*(-Y2*Z4+Z2*Y4)*
     #A1*AUX1-
     #(-Y2*Z4+Z2*Y4)*A2*AUX2+
     #(-Y2*Z4+Z2*Y4)*A3*AUX7-
     #(-X2*Z4+Z2*X4)*A1*AUX2+
     #2*(-X2*Z4+Z2*X4)*A2*AUX4-
     #(-X2*Z4+Z2*X4)*A3*AUX5-
     #(X2*Y4-Y2*X4)*A1*AUX7+
     #(X2*Y4-Y2*X4)*A2*AUX5-
     #2*(X2*Y4-Y2*X4)*A3*AUX6)*COEF
C
      XM(IELEM,3) = (-2*(-Y2*Z3+Z2*Y3)*
     #A1*AUX1+
     #(-Y2*Z3+Z2*Y3)*A2*AUX2-
     #(-Y2*Z3+Z2*Y3)*A3*AUX7+
     #(-X2*Z3+Z2*X3)*A1*AUX2-
     #2*(-X2*Z3+Z2*X3)*A2*AUX4+
     #(-X2*Z3+Z2*X3)*A3*AUX5+
     #(X2*Y3-Y2*X3)*A1*AUX7-
     #(X2*Y3-Y2*X3)*A2*AUX5+
     #2*(X2*Y3-Y2*X3)*A3*AUX6)*COEF
C
      T(IELEM,2) =((Y3*Z4-Y4*Z3)**2*AUX1-
     #(Y3*Z4-Y4*Z3)*(X3*Z4-X4*Z3)*AUX2
     #-(Y3*Z4-Y4*Z3)*(-X3*Y4+X4*Y3)*AUX7+
     #(X3*Z4-X4*Z3)**2*AUX4+
     #(X3*Z4-X4*Z3)*(-X3*Y4+X4*Y3)*AUX5+
     #(-X3*Y4+X4*Y3)**2*AUX6)*2*COEF
C 
      XM(IELEM,4) =  (2*(-Y2*Z4+Z2*Y4)*(Y3*Z4-Y4*Z3)*AUX1-
     #(-Y2*Z4+Z2*Y4)*(X3*Z4-X4*Z3)*AUX2-
     #(-Y2*Z4+Z2*Y4)*(-X3*Y4+X4*Y3)*AUX7
     #-(Y3*Z4-Y4*Z3)*(-X2*Z4+Z2*X4)*AUX2+
     #2*(-X2*Z4+Z2*X4)*(X3*Z4-X4*Z3)*AUX4+
     #(-X2*Z4+Z2*X4)*(-X3*Y4+X4*Y3)*AUX5-
     #(Y3*Z4-Y4*Z3)*(X2*Y4-Y2*X4)*AUX7+
     #(X3*Z4-X4*Z3)*(X2*Y4-Y2*X4)*AUX5+
     #2*(X2*Y4-Y2*X4)*(-X3*Y4+X4*Y3)*AUX6)*COEF
C
      XM(IELEM,5) = (-2*(-Y2*Z3+Z2*Y3)*(Y3*Z4-Y4*Z3)*AUX1+
     #(-Y2*Z3+Z2*Y3)*(X3*Z4-X4*Z3)*AUX2+
     #(-Y2*Z3+Z2*Y3)*(-X3*Y4+X4*Y3)*AUX7+
     #(Y3*Z4-Y4*Z3)*(-X2*Z3+Z2*X3)*AUX2-
     #2*(-X2*Z3+Z2*X3)*(X3*Z4-X4*Z3)*AUX4-
     #(-X2*Z3+Z2*X3)*(-X3*Y4+X4*Y3)*AUX5+
     #(Y3*Z4-Y4*Z3)*(X2*Y3-Y2*X3)*AUX7-
     #(X3*Z4-X4*Z3)*(X2*Y3-Y2*X3)*AUX5-
     #2*(X2*Y3-Y2*X3)*(-X3*Y4+X4*Y3)*AUX6)*COEF
C
      T(IELEM,3) = ((-Y2*Z4+Z2*Y4)**2*AUX1-
     #(-Y2*Z4+Z2*Y4)*(-X2*Z4+Z2*X4)*AUX2-
     #(-Y2*Z4+Z2*Y4)*(X2*Y4-Y2*X4)*AUX7+
     #(-X2*Z4+Z2*X4)**2*AUX4+(-X2*Z4+Z2*X4)*(X2*Y4-Y2*X4)*AUX5+
     #(X2*Y4-Y2*X4)**2*AUX6)*2*COEF
C
      XM(IELEM,6) = (-2*(-Y2*Z3+Z2*Y3)*(-Y2*Z4+Z2*Y4)*AUX1+
     #(-Y2*Z3+Z2*Y3)*(-X2*Z4+Z2*X4)*AUX2+
     #(-Y2*Z3+Z2*Y3)*(X2*Y4-Y2*X4)*AUX7+
     #(-Y2*Z4+Z2*Y4)*(-X2*Z3+Z2*X3)*AUX2-
     #2*(-X2*Z3+Z2*X3)*(-X2*Z4+Z2*X4)*AUX4-
     #(-X2*Z3+Z2*X3)*(X2*Y4-Y2*X4)*AUX5+
     #(-Y2*Z4+Z2*Y4)*(X2*Y3-Y2*X3)*AUX7-
     #(-X2*Z4+Z2*X4)*(X2*Y3-Y2*X3)*AUX5-
     #2*(X2*Y3-Y2*X3)*(X2*Y4-Y2*X4)*AUX6)*COEF
C
      T(IELEM,4) = ((-Y2*Z3+Z2*Y3)**2*AUX1-
     #(-Y2*Z3+Z2*Y3)*(-X2*Z3+Z2*X3)*AUX2-
     #(-Y2*Z3+Z2*Y3)*(X2*Y3-Y2*X3)*AUX7
     #+(-X2*Z3+Z2*X3)**2*AUX4+(-X2*Z3+Z2*X3)*(X2*Y3-Y2*X3)*AUX5+
     #(X2*Y3-Y2*X3)**2*AUX6)*2*COEF
C 
C-----------------------------------------------------------------------
C        
 21   CONTINUE
C
C---------------------------------------------------------------
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,1000) SU%ELM,SV%ELM,SW%ELM
        IF (LNG.EQ.2) WRITE(LU,1001) SU%ELM,SV%ELM,SW%ELM
1000    FORMAT(1X,'MT04TT (BIEF) : MAUVAIS TYPE DE U,V OU W : ',
     *  I6,1X,I6,1X,I6)
1001    FORMAT(1X,'MT04TT (BIEF) : WRONG TYPE OF U,V OR W: ',
     *  I6,1X,I6,1X,I6)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
