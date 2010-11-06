C                       *****************
                        SUBROUTINE MT06TT
C                       *****************
C
     *( T,XM,XMUL,SF,F,X,Y,Z,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.6           10/01/06    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C FONCTION : CALCUL DE LA MATRICE FMATMA POUR DES TETRAEDRES
C
C            LE VECTEUR F PEUT ETRE P0 OU P1
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
C |     X,Y,Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     IKLE1..4   | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |     SW         | -->|  Switch 1:carre magique, 2:classic coeff
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
      USE BIEF, EX_MT06TT => MT06TT
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
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
C
      TYPE(BIEF_OBJ), INTENT(IN) :: SF
      DOUBLE PRECISION, INTENT(IN) :: F(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES 
C     
      DOUBLE PRECISION X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,F1,F2,F3,F4       
      INTEGER I1,I2,I3,I4
      INTEGER IELMF,IELEM
C
      DOUBLE PRECISION XSUR720,JACOB        
C
C***********************************************************************
C
C     DISCRETISATION DU VECTEUR F
      IELMF=SF%ELM

      XSUR720=XMUL/720.D0
C
C-----------------------------------------------------------------------
C
C CAS OU F EST LINEAIRE P1
C
      IF (IELMF .EQ. 31) THEN

C     BOUCLE SUR LES TETRAEDRES
C
          DO 20 IELEM=1,NELEM
C
              I1=IKLE(IELEM,1)
              I2=IKLE(IELEM,2)
              I3=IKLE(IELEM,3)
              I4=IKLE(IELEM,4)
C
              F1 = F(I1)
              F2 = F(I2)
              F3 = F(I3)
              F4 = F(I4) 
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
C     JACOBIEN : Z2*(X3*Y4-X4*Y3)+Y2*(X4*Z3-X3*Z4)+X2*(Y3*Z4-Y4*Z3)
C
C     VOLUME DU TETRAEDRE :
C
C     (Z2*(X3*Y4-X4*Y3)+Y2*(X4*Z3-X3*Z4)+X2*(Y3*Z4-Y4*Z3))/6
C
              JACOB = (   Z2*(X3*Y4-X4*Y3)
     *                  + Y2*(X4*Z3-X3*Z4)
     *                  + X2*(Y3*Z4-Y4*Z3) ) * XSUR720
C
              T(IELEM,1)  = JACOB*(2*F3+6*F1+2*F4+2*F2)
              T(IELEM,2)  = JACOB*(2*F3+2*F1+2*F4+6*F2)
              T(IELEM,3)  = JACOB*2*(F1+3*F3+F2+F4)
              T(IELEM,4)  = JACOB*(2*F1+2*F3+6*F4+2*F2)
C
              XM(IELEM,1) = JACOB*(F3+2*F1+F4+2*F2)
              XM(IELEM,2) = JACOB*(2*F3+F2+2*F1+F4)
              XM(IELEM,3) = JACOB*(F2+F3+2*F4+2*F1)
              XM(IELEM,4) = JACOB*(F1+2*F2+F4+2*F3)
              XM(IELEM,5) = JACOB*(2*F4+F3+F1+2*F2)
              XM(IELEM,6) = JACOB*(F1+2*F3+F2+2*F4)
C
C-----------------------------------------------------------------------
C        
20        CONTINUE 
C
C-----------------------------------------------------------------------
C
C CAS OU F EST CONSTANT PAR ELEMENT P0
C
      ELSEIF (IELMF .EQ. 30) THEN

C     BOUCLE SUR LES TETRAEDRES
C
          DO 30 IELEM=1,NELEM
C
              I1=IKLE(IELEM,1)
              I2=IKLE(IELEM,2)
              I3=IKLE(IELEM,3)
              I4=IKLE(IELEM,4)
C
C             MEME VALEUR POUR LES QUATRE Fi
C             NOTE THAT THE CALCULATIONS FOR T AND XM
C             COULD BE SIMPLIFIED ...
C
              F1 = F(IELEM)
              F2 = F1
              F3 = F1
              F4 = F1
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
C     JACOBIEN : Z2*(X3*Y4-X4*Y3)+Y2*(X4*Z3-X3*Z4)+X2*(Y3*Z4-Y4*Z3)
C
C     VOLUME DU TETRAEDRE :
C
C     (Z2*(X3*Y4-X4*Y3)+Y2*(X4*Z3-X3*Z4)+X2*(Y3*Z4-Y4*Z3))/6
C
              JACOB = (   Z2*(X3*Y4-X4*Y3)
     *                  + Y2*(X4*Z3-X3*Z4)
     *                  + X2*(Y3*Z4-Y4*Z3) ) * XSUR720
C
              T(IELEM,1)  = JACOB*(2*F3+6*F1+2*F4+2*F2)
              T(IELEM,2)  = JACOB*(2*F3+2*F1+2*F4+6*F2)
              T(IELEM,3)  = JACOB*2*(F1+3*F3+F2+F4)
              T(IELEM,4)  = JACOB*(2*F1+2*F3+6*F4+2*F2)
C
              XM(IELEM,1) = JACOB*(F3+2*F1+F4+2*F2)
              XM(IELEM,2) = JACOB*(2*F3+F2+2*F1+F4)
              XM(IELEM,3) = JACOB*(F2+F3+2*F4+2*F1)
              XM(IELEM,4) = JACOB*(F1+2*F2+F4+2*F3)
              XM(IELEM,5) = JACOB*(2*F4+F3+F1+2*F2)
              XM(IELEM,6) = JACOB*(F1+2*F3+F2+2*F4)
C
C-----------------------------------------------------------------------
C        
30    CONTINUE 
C
C-----------------------------------------------------------------------
      ELSE
C
          IF (LNG.EQ.1) WRITE(LU,100) IELMF,SF%NAME
          IF (LNG.EQ.2) WRITE(LU,101) IELMF,SF%NAME
100       FORMAT(1X,'MT06TT (BIEF) :',/,
     *        1X,'DISCRETISATION DE F NON PREVUE : ',1I6,
     *        1X,'NOM REEL : ',A6)
101       FORMAT(1X,'MT06TT (BIEF) :',/,
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
      END SUBROUTINE MT06TT
