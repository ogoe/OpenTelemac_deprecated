C                       *****************
                        SUBROUTINE MT08TT
C                       *****************
C
     *( T,XM,XMUL,X,Y,Z,SF,F,IKLE,NELEM,NELMAX)
C
C***********************************************************************
C BIEF VERSION 5.3           21/03/02  J-M HERVOUET (LNH) 01 30 87 80 18
C                                        
C***********************************************************************
C
C FONCTION : 
C
C-----------------------------------------------------------------------
C
C      FONCTION:
C      ========:
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C                                 /                  
C                    A    = XMUL / P   F . GRAD(P ) * J(X,Y) DXDY
C                     I J       /S  J            I               
C
C    PAR MAILLE ELEMENTAIRE.
C
C
C    !!!!!!!!!!!!!!!!!!!!!!   ATTENTION
C
C    SEULE LA COMPOSANTE Z EST TRAITEE ICI
C
C    !!!!!!!!!!!!!!!!!!!!!!   ATTENTION
C
C    J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C    L'ELEMENT EST LE TRIANGLE P1
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |     T,XM       |<-- |  ELEMENTS DE LA MATRICE
C |     XMUL       | -->|  FACTEUR MULTIPLICATIF
C |     SF,SG,SH   | -->|  STRUCTURES DE F,G ET H.
C |     SU,SV,SW   | -->|  STRUCTURES DE U,V ET W.
C |     F,G,H      | -->|  FONCTIONS INTERVENANT DANS LE CALCUL DE LA
C |                |    |  MATRICE.
C |     U,V,W      | -->|  COMPOSANTES D'UN VECTEUR INTERVENANT DANS LE
C |                |    |  CALCUL DE LA MATRICE.
C |     X,Y,Z      | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
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
      USE BIEF, EX_MT08TT => MT08TT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER,INTENT(IN) :: NELEM,NELMAX
      INTEGER,INTENT(IN) :: IKLE(NELMAX,4)
C
      DOUBLE PRECISION,INTENT(INOUT) :: T(NELMAX,4),XM(NELMAX,12)
C
      DOUBLE PRECISION,INTENT(IN) :: XMUL
      DOUBLE PRECISION,INTENT(IN) :: F(*),X(*),Y(*),Z(*)
C
C     STRUCTURE DE F
C
      TYPE(BIEF_OBJ),INTENT(IN) :: SF
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      DOUBLE PRECISION X2,X3,X4,Y2,Y3,Y4,F1,F2,F3,F4,XSUR120
C
      INTEGER I1,I2,I3,I4,IELEM
C
C**********************************************************************
C
      XSUR120 = XMUL/120.D0
C
      IF(SF%ELM.NE.31.AND.SF%ELM.NE.51) THEN
        IF (LNG.EQ.1) WRITE(LU,1000) SF%ELM
        IF (LNG.EQ.2) WRITE(LU,1001) SF%ELM
1000    FORMAT(1X,'MT08TT (BIEF) : TYPE DE F NON PREVU : ',I6)
1001    FORMAT(1X,'MT08TT (BIEF): TYPE OF F NOT IMPLEMENTED: ',I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C   BOUCLE SUR LES ELEMENTS
C
      DO 1 IELEM=1,NELEM
C
         I1 = IKLE(IELEM,1)
         I2 = IKLE(IELEM,2)
         I3 = IKLE(IELEM,3)
         I4 = IKLE(IELEM,4)
C
         X2 = X(I2)-X(I1)
         X3 = X(I3)-X(I1)
         X4 = X(I4)-X(I1)
C
         Y2 = Y(I2)-Y(I1)
         Y3 = Y(I3)-Y(I1)
         Y4 = Y(I4)-Y(I1)
C
         F1 = F(I1)
         F2 = F(I2)
         F3 = F(I3)
         F4 = F(I4)
C
         T(IELEM,1)=(
     # X2*Y4*F3-X3*Y4*F3+X4*Y3*F2-Y2*X4*F2+X2*Y4*F2+2*Y2*X3*F1
     #-2*X3*Y4*F1+Y2*X3*F2+Y2*X3*F3+2*X4*Y3*F1+2*X2*Y4*F1-2*Y2*X4*F1
     #-2*X2*Y3*F1-Y2*X4*F3-X2*Y3*F2-X2*Y3*F3+X4*Y3*F3-X3*Y4*F4
     #+X2*Y4*F4-Y2*X4*F4-X2*Y3*F4+Y2*X3*F4+X4*Y3*F4-X3*Y4*F2 )*XSUR120
C
         T(IELEM,2)=(
     # X3*Y4*F3-2*X4*Y3*F2+X3*Y4*F1-X4*Y3*F1-X4*Y3*F3+X3*Y4*F4
     #-X4*Y3*F4+2*X3*Y4*F2                        )*XSUR120
C
         T(IELEM,3)=(-X2*Y4+Y2*X4)*(F1+2*F3+F2+F4) *XSUR120
C
         T(IELEM,4)=(
     # -Y2*X3*F1-Y2*X3*F2-Y2*X3*F3+X2*Y3*F1+X2*Y3*F2+X2*Y3*F3+
     #2*X2*Y3*F4-2*Y2*X3*F4                       )*XSUR120
C
         XM(IELEM,01)=(
     # X2*Y4*F3-X3*Y4*F3+2*X4*Y3*F2-2*Y2*X4*F2+2*X2*Y4*F2+Y2*X3*F1
     #-X3*Y4*F1+2*Y2*X3*F2+Y2*X3*F3+X4*Y3*F1+X2*Y4*F1-Y2*X4*F1-X2*Y3
     #*F1-Y2*X4*F3-2*X2*Y3*F2-X2*Y3*F3+X4*Y3*F3
     #-X3*Y4*F4+X2*Y4*F4-Y2*X4*F4
     #-X2*Y3*F4+Y2*X3*F4+X4*Y3*F4-2*X3*Y4*F2    )*XSUR120
C
         XM(IELEM,02)=(
     #-X3*Y4+X4*Y3+X2*Y4-Y2*X4-X2*Y3+Y2*X3)*(F1+2*F3+F2+F4)*XSUR120
C
         XM(IELEM,03)=(
     # X2*Y4*F3-X3*Y4*F3+X4*Y3*F2-Y2*X4*F2+X2*Y4*F2+Y2*X3*F1-X3*Y4*F1
     #+Y2*X3*F2+Y2*X3*F3+X4*Y3*F1+X2*Y4*F1-Y2*X4*F1-X2*Y3*F1-Y2*X4*F3
     #-X2*Y3*F2-X2*Y3*F3+X4*Y3*F3-2*X3*Y4*F4+2*X2*Y4*F4-2*Y2*X4*F4-2
     #*X2*Y3*F4+2*Y2*X3*F4+2*X4*Y3*F4-X3*Y4*F2)*XSUR120
C
         XM(IELEM,04)= -(-X3*Y4+X4*Y3)*(F1+2*F3+F2+F4)*XSUR120
C
         XM(IELEM,05)=( X3*Y4*F3-X4*Y3*F2+X3*Y4*F1-X4*Y3*F1-X4*Y3*F3
     #                  +2*X3*Y4*F4-2*X4*Y3*F4+X3*Y4*F2)*XSUR120
C
         XM(IELEM,06)=( -X2*Y4*F3+Y2*X4*F2-X2*Y4*F2-X2*Y4*F1+Y2*X4*F1
     #                  +Y2*X4*F3-2*X2*Y4*F4+2*Y2*X4*F4)*XSUR120
C
         XM(IELEM,07)=( X3*Y4*F3-X4*Y3*F2+2*X3*Y4*F1-2*X4*Y3*F1
     #                 -X4*Y3*F3+X3*Y4*F4-X4*Y3*F4+X3*Y4*F2)*XSUR120
C
         XM(IELEM,08)=( -X2*Y4*F3+Y2*X4*F2-X2*Y4*F2-2*X2*Y4*F1
     #                  +2*Y2*X4*F1+Y2*X4*F3-X2*Y4*F4+Y2*X4*F4)*XSUR120
C
         XM(IELEM,09)=( -X2*Y4*F3+2*Y2*X4*F2-2*X2*Y4*F2-X2*Y4*F1
     #                  +Y2*X4*F1+Y2*X4*F3-X2*Y4*F4+Y2*X4*F4 )*XSUR120
C
         XM(IELEM,10)=( -2*Y2*X3*F1-Y2*X3*F2-Y2*X3*F3+2*X2*Y3*F1
     #                  +X2*Y3*F2+X2*Y3*F3+X2*Y3*F4-Y2*X3*F4)*XSUR120
C
         XM(IELEM,11)=(-Y2*X3*F1-2*Y2*X3*F2-Y2*X3*F3+X2*Y3*F1
     #                 +2*X2*Y3*F2+X2*Y3*F3+X2*Y3*F4-Y2*X3*F4)*XSUR120
C
         XM(IELEM,12)= -(-X2*Y3+Y2*X3)*(F1+2*F3+F2+F4)*XSUR120
C
1     CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END
