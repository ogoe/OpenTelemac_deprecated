C                       *****************
                        SUBROUTINE MT01TT
C                       *****************
C
     *( T,XM,XMUL,X,Y,Z,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.3         04/01/02    J-M HERVOUET (LNH) 01 30 87 80 18
C                                       
C***********************************************************************
C
C FONCTION : CALCUL D'UNE MATRICE DE MASSE EN TETRAEDRES
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
C |     T,XM       |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     X,Y,Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |     IKLE       | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |     NELEM      | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |     NELMAX     | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS D'UN MAILLAGE ADAPTATIF)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
      USE BIEF, EX_MT01TT => MT01TT
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELEM,NELMAX
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,4)
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,4),XM(NELMAX,6)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL
      DOUBLE PRECISION, INTENT(IN)    :: X(*),Y(*),Z(*)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES 
C     
      DOUBLE PRECISION X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,VOLSUR20      
      INTEGER I1,I2,I3,I4,IELEM
C
      DOUBLE PRECISION XSUR120        
C
C***********************************************************************
C
      XSUR120=XMUL/120.D0
C
C-----------------------------------------------------------------------
C
C     BOUCLE SUR LES TETRAEDRES
C
      DO 20 IELEM=1,NELEM
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
C     TERMES EXTRA-DIAGONAUX
C
C     VOLUME DU TETRAEDRE :
C
C     (Z2*(X3*Y4-X4*Y3)+Y2*(X4*Z3-X3*Z4)+X2*(Y3*Z4-Y4*Z3))/6
C
C     XM(IELEM,1) = VOLUME / 20
C
C     SOMME DES TERMES (AVEC LES SYMETRIQUES) = VOLUME DU TETRAEDRE
C
      VOLSUR20 = 
     *(Z2*(X3*Y4-X4*Y3)+Y2*(X4*Z3-X3*Z4)+X2*(Y3*Z4-Y4*Z3))*XSUR120
      XM(IELEM,1) = MAX(VOLSUR20,1.D-4)
      XM(IELEM,2) = XM(IELEM,1)
      XM(IELEM,3) = XM(IELEM,1)
      XM(IELEM,4) = XM(IELEM,1)
      XM(IELEM,5) = XM(IELEM,1)
      XM(IELEM,6) = XM(IELEM,1)
C
C     TERMES DIAGONAUX
C    
      T(IELEM,1) = 2 * XM(IELEM,1)
      T(IELEM,2) = T(IELEM,1)
      T(IELEM,3) = T(IELEM,1)
      T(IELEM,4) = T(IELEM,1)
C
C-----------------------------------------------------------------------
C        
20    CONTINUE 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
