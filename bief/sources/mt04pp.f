C                       *****************
                        SUBROUTINE MT04PP
C                       *****************
C
     *( T,XM,XMUL,SU,SV,SW,U,V,W,X,Y,Z,SURFAC,IKLE,NELEM,NELMAX,FORMUL)
C
C***********************************************************************
C BIEF VERSION 5.6         07/07/05   J-M HERVOUET (LNHE) 01 30 87 80 18
C                                       F  LEPEINTRE (LNH) 30 87 78 54
C
C 16/09/2005: VERTICAL UPWIND REMOVED
C
C***********************************************************************
C
C FONCTION : CONSTRUCTION DE LA MATRICE U GRAG(PSII) U GRAD(PSIJ)
C
C-----------------------------------------------------------------------
C
C      FONCTION:
C      ========:
C
C    CE SOUS-PROGRAMME CALCULE LES COEFFICIENTS DE LA MATRICE SUIVANTE:
C
C    MAUGUG2 :
C
C                            / ->                 ->
C               A    = XMUL /  U   . GRAD  (P ) * U . GRAD  (P ) * J(X,Y) DXDY
C                I J       /S    2D      2D  I     2D     2D  J
C
C
C    MAUGUG1 : SEE COMPONENTS F0 AND G0
C
C                            / ->                 ->
C               A    = XMUL /  F   . GRAD  (P ) * U . GRAD  (P ) * J(X,Y) DXDY
C                I J       /S    2D      2D  I     2D     2D  J
C
C
C    PAR MAILLE ELEMENTAIRE.
C
C    J(X,Y) : JACOBIEN DE LA TRANSFORMATION ISOPARAMETRIQUE
C
C    L'ELEMENT EST LE PRISME P1
C
C    SIMPLIFICATIONS: AN AVERAGE HEIGHT OF THE PRISM IS TAKEN
C                     VELOCITIES ARE CONSIDERED CONSTANT ON THE ELEMENT
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
      USE BIEF, EX_MT04PP => MT04PP
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NELMAX
      INTEGER, INTENT(IN) :: IKLE(NELMAX,6)
C
      DOUBLE PRECISION, INTENT(INOUT) :: T(NELMAX,6),XM(NELMAX,30)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
C
      DOUBLE PRECISION, INTENT(IN) :: XMUL
      DOUBLE PRECISION, INTENT(IN) :: U(*),V(*),W(*)
C
C     STRUCTURES DE      U, V, W
      TYPE(BIEF_OBJ), INTENT(IN) :: SU,SV,SW
C
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*),Z(*)
C
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     DECLARATIONS SPECIFIQUES A CE SOUS-PROGRAMME
C
      DOUBLE PRECISION X2,X3,Y2,Y3
      DOUBLE PRECISION U0,V0,F0,G0,HH,C,DX,SURNORMU
C
      INTEGER I1,I2,I3,I4,I5,I6,IELEM
C
      INTRINSIC SQRT
C
C**********************************************************************
C
      IF(SU%ELM.NE.41) THEN
        IF (LNG.EQ.1) WRITE(LU,1000) SU%ELM
        IF (LNG.EQ.2) WRITE(LU,1001) SU%ELM
1000    FORMAT(1X,'MT04PP (BIEF) : TYPE DE U NON PREVU : ',I6)
1001    FORMAT(1X,'MT04PP (BIEF) : TYPE OF U NOT IMPLEMENTED: ',I6)
        CALL PLANTE(1)
        STOP
      ENDIF
      IF(SV%ELM.NE.41) THEN
        IF (LNG.EQ.1) WRITE(LU,2000) SV%ELM
        IF (LNG.EQ.2) WRITE(LU,2001) SV%ELM
2000    FORMAT(1X,'MT04PP (BIEF) : TYPE DE V NON PREVU : ',I6)
2001    FORMAT(1X,'MT04PP (BIEF) : TYPE OF V NOT IMPLEMENTED: ',I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C     HERE WSCONV WHICH IS IN FACT DEFINED PER LAYER
      IF(SW%ELM.NE.41) THEN
        IF (LNG.EQ.1) WRITE(LU,3000) SW%ELM
        IF (LNG.EQ.2) WRITE(LU,3001) SW%ELM
3000    FORMAT(1X,'MT04PP (BIEF) : TYPE DE W NON PREVU : ',I6)
3001    FORMAT(1X,'MT04PP (BIEF) : TYPE OF W NOT IMPLEMENTED: ',I6)
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IF(FORMUL(1:7).EQ.'MAUGUG2') THEN
C
C     BOUCLE SUR LES ELEMENTS 3D
C
      DO IELEM = 1 , NELEM
C
      I1 = IKLE(IELEM,1)
      I2 = IKLE(IELEM,2)
      I3 = IKLE(IELEM,3)
      I4 = IKLE(IELEM,4)
      I5 = IKLE(IELEM,5)
      I6 = IKLE(IELEM,6)
C
      X2  =  X(I2) - X(I1)
      X3  =  X(I3) - X(I1)
      Y2  =  Y(I2) - Y(I1)
      Y3  =  Y(I3) - Y(I1)
C
      U0 = (U(I1)+U(I2)+U(I3)+U(I4)+U(I5)+U(I6))/6.D0
      V0 = (V(I1)+V(I2)+V(I3)+V(I4)+V(I5)+V(I6))/6.D0
C     AVERAGE HEIGHT OF PRISM
      HH = (Z(I4)-Z(I1)+Z(I5)-Z(I2)+Z(I6)-Z(I3))/3.D0
      C  = XMUL*HH/24.D0/SURFAC(IELEM)
C
      T(IELEM,1)=2*C*(U0*Y3-U0*Y2-V0*X3+V0*X2)**2
      T(IELEM,2)=2*C*(-U0*Y3+V0*X3)**2
      T(IELEM,3)=2*C*(-U0*Y2+V0*X2)**2
      T(IELEM,4)=T(IELEM,1)
      T(IELEM,5)=T(IELEM,2)
      T(IELEM,6)=T(IELEM,3)
C
      XM(IELEM,01)= 2*(-U0*Y3+V0*X3)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,02)=-2*(-U0*Y2+V0*X2)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,03)=(U0*Y3-U0*Y2-V0*X3+V0*X2)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,04)=(-U0*Y3+V0*X3)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,05)=-(-U0*Y2+V0*X2)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,06)=-2*(-U0*Y2+V0*X2)*(-U0*Y3+V0*X3)*C
      XM(IELEM,07)=(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-U0*Y3+V0*X3)*C
      XM(IELEM,08)=(-U0*Y3+V0*X3)*(-U0*Y3+V0*X3)*C
      XM(IELEM,09)=-(-U0*Y2+V0*X2)*(-U0*Y3+V0*X3)*C
      XM(IELEM,10)=-(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-U0*Y2+V0*X2)*C
      XM(IELEM,11)=-(-U0*Y3+V0*X3)*(-U0*Y2+V0*X2)*C
      XM(IELEM,12)=(-U0*Y2+V0*X2)*(-U0*Y2+V0*X2)*C
      XM(IELEM,13)= 2*(-U0*Y3+V0*X3)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,14)=-2*(-U0*Y2+V0*X2)*(U0*Y3-U0*Y2-V0*X3+V0*X2)*C
      XM(IELEM,15)=-2*(-U0*Y2+V0*X2)*(-U0*Y3+V0*X3)*C      
C
      ENDDO
C
      ELSEIF(FORMUL(1:7).EQ.'MAUGUG1') THEN
C
      DO IELEM = 1 , NELEM
C
      I1 = IKLE(IELEM,1)
      I2 = IKLE(IELEM,2)
      I3 = IKLE(IELEM,3)
      I4 = IKLE(IELEM,4)
      I5 = IKLE(IELEM,5)
      I6 = IKLE(IELEM,6)
C
      X2  =  X(I2) - X(I1)
      X3  =  X(I3) - X(I1)
      Y2  =  Y(I2) - Y(I1)
      Y3  =  Y(I3) - Y(I1)
C     VELOCITIES CONSIDERED CONSTANT
      U0 = (U(I1)+U(I2)+U(I3)+U(I4)+U(I5)+U(I6))/6.D0
      V0 = (V(I1)+V(I2)+V(I3)+V(I4)+V(I5)+V(I6))/6.D0
!
      SURNORMU=1.D0/MAX(SQRT(U0**2+V0**2),1.D-8)
      DX=SQRT(2.D0*SURFAC(IELEM))
!     ICI F0 ET G0 DOIVENT ETRE 
!     DX*U/(2*NORME(U)) ET DY*V/(2*NORME(U)) MAIS ON PREND DY=DX
!     APPROCHE PAR RACINE(2*SURFACE DU TRIANGLE)
      F0 = 0.5D0*DX*U0*SURNORMU
      G0 = 0.5D0*DX*V0*SURNORMU
C     AVERAGE HEIGHT OF PRISM
      HH = (Z(I4)-Z(I1)+Z(I5)-Z(I2)+Z(I6)-Z(I3))/3.D0
      C  = XMUL*HH/24.D0/SURFAC(IELEM)
C
      T(IELEM,1)=2*(U0*Y3-U0*Y2-V0*X3+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      T(IELEM,2)=2*(-U0*Y3+V0*X3)*(-F0*Y3+G0*X3)*C
      T(IELEM,3)=2*(-U0*Y2+V0*X2)*(-F0*Y2+G0*X2)*C
      T(IELEM,4)=T(IELEM,1)
      T(IELEM,5)=T(IELEM,2)
      T(IELEM,6)=T(IELEM,3)
C
      XM(IELEM,01)=2*(-U0*Y3+V0*X3)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,02)=-2*(-U0*Y2+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,03)=(U0*Y3-U0*Y2-V0*X3+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,04)=(-U0*Y3+V0*X3)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,05)=-(-U0*Y2+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,06)=-2*(-U0*Y2+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,07)=(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,08)=(-U0*Y3+V0*X3)*(-F0*Y3+G0*X3)*C
      XM(IELEM,09)=-(-U0*Y2+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,10)=-(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y2+G0*X2)*C
      XM(IELEM,11)=-(-U0*Y3+V0*X3)*(-F0*Y2+G0*X2)*C
      XM(IELEM,12)=(-U0*Y2+V0*X2)*(-F0*Y2+G0*X2)*C
      XM(IELEM,13)=2*(-U0*Y3+V0*X3)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,14)=-2*(-U0*Y2+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,15)=-2*(-U0*Y2+V0*X2)*(-F0*Y3+G0*X3)*C 
C
      XM(IELEM,16)= 2*(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,17)= -2*(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y2+G0*X2)*C
      XM(IELEM,21)= -2*(-U0*Y3+V0*X3)*(-F0*Y2+G0*X2)*C
      XM(IELEM,18)=(U0*Y3-U0*Y2-V0*X3+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,22)= (-U0*Y3+V0*X3)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,25)= -(-U0*Y2+V0*X2)*(F0*Y3-F0*Y2-G0*X3+G0*X2)*C
      XM(IELEM,19)= (U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,23)= (-U0*Y3+V0*X3)*(-F0*Y3+G0*X3)*C
      XM(IELEM,26)= -(-U0*Y2+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,28)= 2*(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y3+G0*X3)*C
      XM(IELEM,20)= -(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y2+G0*X2)*C
      XM(IELEM,24)= -(-U0*Y3+V0*X3)*(-F0*Y2+G0*X2)*C
      XM(IELEM,27)= (-U0*Y2+V0*X2)*(-F0*Y2+G0*X2)*C
      XM(IELEM,29)= -2*(U0*Y3-U0*Y2-V0*X3+V0*X2)*(-F0*Y2+G0*X2)*C
      XM(IELEM,30)= -2*(-U0*Y3+V0*X3)*(-F0*Y2+G0*X2)*C  
C
      ENDDO
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,4000) FORMUL
        IF (LNG.EQ.2) WRITE(LU,4001) FORMUL
4000    FORMAT(1X,'MT04PP (BIEF) : FORMULE NON PREVUE : ',A16)
4001    FORMAT(1X,'MT04PP (BIEF) : UNEXPECTED FORMULA: ',A16)
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
