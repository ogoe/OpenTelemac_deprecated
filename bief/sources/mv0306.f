C                       *****************
                        SUBROUTINE MV0306
C                       *****************
C
     *(OP, X , DA,TYPDIA,
     * XA12,XA13,XA14,XA15,XA16,XA21,XA23,XA24,
     * XA25,XA26,XA31,XA32,XA34,XA35,XA36,
     * TYPEXT, Y,C,
     * IKLE1,IKLE2,IKLE3,IKLE4,IKLE5,IKLE6,
     * NPOIN,NPT2,NELEM,W1,W2,W3,W4,W5,W6)
C
C***********************************************************************
C BIEF VERSION 5.9           31/03/08    ALGIANE FROEHLY (MATMECA)
C                                        J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : OPERATIONS MATRICE VECTEUR POUR TRIANGLES P1*P2
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET LA MATRICE M. LE RESULTAT
C   EST LE VECTEUR X.
C
C   CES OPERATIONS SONT DIFFERENTES SUIVANT LE TYPE DE DIAGONALE
C   ET LE TYPE DES TERMES EXTRADIAGONAUX.
C
C  OPERATIONS PROGRAMMEES :
C
C      OP = 'X=AY    '  : X = AY
C      OP = 'X=-AY   '  : X = - AY
C      OP = 'X=X+AY  '  : X = X + AY
C      OP = 'X=X-AY  '  : X = X - AY
C      OP = 'X=X+CAY '  : X = X + C AY
C      OP = 'X=TAY   '  : X = TA Y (TRANSPOSEE DE A)
C      OP = 'X=-TAY  '  : X = - TA Y (- TRANSPOSEE DE A)
C      OP = 'X=X+TAY '  : X = X + TA Y
C      OP = 'X=X-TAY '  : X = X - TA Y
C      OP = 'X=X+CTAY'  : X = X + C TA Y
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      OP        | -->| OPERATION A EFFECTUER
C |      X         |<-- | VECTEUR IMAGE
C |      DA        | -->| DIAGONALE DE LA MATRICE
C |      TYPDIA    | -->| TYPE DE LA DIAGONALE (CHAINE DE CARACTERES)
C |                |    | TYPDIA = 'Q' : DIAGONALE QUELCONQUE
C |                |    | TYPDIA = 'I' : DIAGONALE IDENTITE.
C |                |    | TYPDIA = '0' : DIAGONALE NULLE.
C |      XA12,.... | -->| TERMES EXTRADIAGONAUX ELEMENTAIRES
C |      TYPEXT    | -->| TYPEXT = 'Q' : QUELCONQUES.
C |                |    | TYPEXT = 'S' : SYMETRIQUES.
C |                |    | TYPEXT = '0' : NULS.
C |      Y         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |      IKLE1,..  | -->| CORRESPONDANCE NUMEROTATIONS LOCALE ET GLOBALE
C |      NPOIN     | -->| NOMBRE DE POINTS.
C |      NELEM     | -->| NOMBRE D'ELEMENTS.
C |      W1,..     |<-- | TABLEAUX DE TRAVAIL DE DIMENSION NELEM
C |                |    | QUI CONTIENDRONT UNE PARTIE DU RESULTAT SOUS
C |                |    | FORME NON ASSEMBLEE.
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : ASSVEC , OV
C
C***********************************************************************
C
      USE BIEF!, EX_MV0306 => MV0306
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NELEM,NPOIN,NPT2
C
      INTEGER, INTENT(IN) :: IKLE1(*),IKLE2(*),IKLE3(*)
      INTEGER, INTENT(IN) :: IKLE4(*),IKLE5(*),IKLE6(*)
C
      DOUBLE PRECISION, INTENT(INOUT) :: W1(*),W2(*),W3(*)
      DOUBLE PRECISION, INTENT(INOUT) :: W4(*),W5(*),W6(*)
      DOUBLE PRECISION, INTENT(IN) :: Y(*),DA(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      DOUBLE PRECISION, INTENT(IN) :: XA12(*),XA13(*),XA14(*),XA15(*)
      DOUBLE PRECISION, INTENT(IN) :: XA16(*),XA21(*),XA23(*),XA24(*)
      DOUBLE PRECISION, INTENT(IN) :: XA25(*),XA26(*),XA31(*),XA32(*)
      DOUBLE PRECISION, INTENT(IN) :: XA34(*),XA35(*),XA36(*)
      DOUBLE PRECISION, INTENT(IN) :: C
C
      CHARACTER(LEN=8), INTENT(IN) :: OP
      CHARACTER(LEN=1), INTENT(IN) :: TYPDIA,TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM
      DOUBLE PRECISION Z(1)
C
C-----------------------------------------------------------------------
C
      IF(OP(1:8).EQ.'X=AY    ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 20 IELEM = 1 , NELEM
             W1(IELEM) =     XA12(IELEM) * Y(IKLE2(IELEM))
     *                     + XA13(IELEM) * Y(IKLE3(IELEM))
     *                     + XA14(IELEM) * Y(IKLE4(IELEM))
     *                     + XA15(IELEM) * Y(IKLE5(IELEM))
     *                     + XA16(IELEM) * Y(IKLE6(IELEM)) 
             W2(IELEM) =     XA21(IELEM) * Y(IKLE1(IELEM))
     *                     + XA23(IELEM) * Y(IKLE3(IELEM))
     *                     + XA24(IELEM) * Y(IKLE4(IELEM))
     *                     + XA25(IELEM) * Y(IKLE5(IELEM))
     *                     + XA26(IELEM) * Y(IKLE6(IELEM))
             W3(IELEM) =     XA31(IELEM) * Y(IKLE1(IELEM))
     *                     + XA32(IELEM) * Y(IKLE2(IELEM))
     *                     + XA34(IELEM) * Y(IKLE4(IELEM))
     *                     + XA35(IELEM) * Y(IKLE5(IELEM))
     *                     + XA36(IELEM) * Y(IKLE6(IELEM))


20     CONTINUE
C
         ELSEIF(TYPEXT(1:1).EQ.'0') THEN
C
           CALL OV ('X=C     ', W1 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W2 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W3 , Y , Z , 0.D0 , NELEM )
C
         ELSE
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE :
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=YZ    ', X , Y , DA , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=Y     ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'0') THEN
           CALL OV ('X=C     ', X , Y , Z  , 0.D0 , NPOIN )
         ELSE
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=-AY   ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 21 IELEM = 1 , NELEM
             W1(IELEM) =   - XA12(IELEM) * Y(IKLE2(IELEM))
     *                     - XA13(IELEM) * Y(IKLE3(IELEM))
     *                     - XA14(IELEM) * Y(IKLE4(IELEM))
     *                     - XA15(IELEM) * Y(IKLE5(IELEM))
     *                     - XA16(IELEM) * Y(IKLE6(IELEM)) 
             W2(IELEM) =   - XA21(IELEM) * Y(IKLE1(IELEM))
     *                     - XA23(IELEM) * Y(IKLE3(IELEM))
     *                     - XA24(IELEM) * Y(IKLE4(IELEM))
     *                     - XA25(IELEM) * Y(IKLE5(IELEM))
     *                     - XA26(IELEM) * Y(IKLE6(IELEM))   
             W3(IELEM) =   - XA31(IELEM) * Y(IKLE1(IELEM))
     *                     - XA32(IELEM) * Y(IKLE2(IELEM))
     *                     - XA34(IELEM) * Y(IKLE4(IELEM))
     *                     - XA35(IELEM) * Y(IKLE5(IELEM))
     *                     - XA36(IELEM) * Y(IKLE6(IELEM)) 
21         CONTINUE
C
         ELSEIF(TYPEXT(1:1).EQ.'0') THEN
C
           CALL OV ('X=C     ', W1 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W2 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W3 , Y , Z , 0.D0 , NELEM )
C
         ELSE
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE :
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=-YZ   ', X , Y , DA , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=-Y    ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'0') THEN
           CALL OV ('X=C     ', X , Y , Z  , 0.D0 , NPOIN )
         ELSE
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+AY  ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 40  IELEM = 1 , NELEM
             W1(IELEM) = W1(IELEM) + XA12(IELEM) * Y(IKLE2(IELEM))
     *                             + XA13(IELEM) * Y(IKLE3(IELEM))
     *                             + XA14(IELEM) * Y(IKLE4(IELEM))
     *                             + XA15(IELEM) * Y(IKLE5(IELEM))
     *                             + XA16(IELEM) * Y(IKLE6(IELEM)) 
             W2(IELEM) = W2(IELEM) + XA21(IELEM) * Y(IKLE1(IELEM))
     *                             + XA23(IELEM) * Y(IKLE3(IELEM))
     *                             + XA24(IELEM) * Y(IKLE4(IELEM))
     *                             + XA25(IELEM) * Y(IKLE5(IELEM))
     *                             + XA26(IELEM) * Y(IKLE6(IELEM)) 
             W3(IELEM) = W3(IELEM) + XA31(IELEM) * Y(IKLE1(IELEM))
     *                             + XA32(IELEM) * Y(IKLE2(IELEM))
     *                             + XA34(IELEM) * Y(IKLE4(IELEM))
     *                             + XA35(IELEM) * Y(IKLE5(IELEM))
     *                             + XA36(IELEM) * Y(IKLE6(IELEM)) 
40         CONTINUE
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE :
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=X+YZ  ', X , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=X+Y   ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X-AY  ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 60 IELEM = 1 , NELEM
             W1(IELEM) = W1(IELEM) - XA12(IELEM) * Y(IKLE2(IELEM))
     *                             - XA13(IELEM) * Y(IKLE3(IELEM))
     *                             - XA14(IELEM) * Y(IKLE4(IELEM))
     *                             - XA15(IELEM) * Y(IKLE5(IELEM))
     *                             - XA16(IELEM) * Y(IKLE6(IELEM)) 
             W2(IELEM) = W2(IELEM) - XA21(IELEM) * Y(IKLE1(IELEM))
     *                             - XA23(IELEM) * Y(IKLE3(IELEM))
     *                             - XA24(IELEM) * Y(IKLE4(IELEM))
     *                             - XA25(IELEM) * Y(IKLE5(IELEM))
     *                             - XA26(IELEM) * Y(IKLE6(IELEM)) 
             W3(IELEM) = W3(IELEM) - XA31(IELEM) * Y(IKLE1(IELEM))
     *                             - XA32(IELEM) * Y(IKLE2(IELEM))
     *                             - XA34(IELEM) * Y(IKLE4(IELEM))
     *                             - XA35(IELEM) * Y(IKLE5(IELEM))
     *                             - XA36(IELEM) * Y(IKLE6(IELEM)) 
60         CONTINUE
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE :
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=X-YZ  ', X , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=X-Y   ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(0)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+CAY ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 80 IELEM = 1 , NELEM
             W1(IELEM) = W1(IELEM)
     *               + C * (      XA12(IELEM) * Y(IKLE2(IELEM))
     *                          + XA13(IELEM) * Y(IKLE3(IELEM))
     *                          + XA14(IELEM) * Y(IKLE4(IELEM)) 
     *                          + XA15(IELEM) * Y(IKLE5(IELEM))
     *                          + XA16(IELEM) * Y(IKLE6(IELEM)) )
             W2(IELEM) = W2(IELEM)
     *               + C * (      XA21(IELEM) * Y(IKLE1(IELEM))
     *                          + XA23(IELEM) * Y(IKLE3(IELEM))
     *                          + XA24(IELEM) * Y(IKLE4(IELEM)) 
     *                          + XA25(IELEM) * Y(IKLE5(IELEM))
     *                          + XA26(IELEM) * Y(IKLE6(IELEM)) )
             W3(IELEM) = W3(IELEM)
     *               + C * (      XA31(IELEM) * Y(IKLE1(IELEM))
     *                          + XA32(IELEM) * Y(IKLE2(IELEM))
     *                          + XA34(IELEM) * Y(IKLE4(IELEM)) 
     *                          + XA35(IELEM) * Y(IKLE5(IELEM))
     *                          + XA36(IELEM) * Y(IKLE6(IELEM)) )
80         CONTINUE
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE :
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=X+CYZ  ', X , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=X+CY   ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=TAY   ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO IELEM = 1 , NELEM
             W1(IELEM) =   + XA21(IELEM) * Y(IKLE2(IELEM))
     *                     + XA31(IELEM) * Y(IKLE3(IELEM))
             W2(IELEM) =   + XA12(IELEM) * Y(IKLE1(IELEM))
     *                     + XA32(IELEM) * Y(IKLE3(IELEM))
             W3(IELEM) =   + XA13(IELEM) * Y(IKLE1(IELEM))
     *                     + XA23(IELEM) * Y(IKLE2(IELEM))
             W4(IELEM) =   + XA14(IELEM) * Y(IKLE1(IELEM))
     *                     + XA24(IELEM) * Y(IKLE2(IELEM))
     *                     + XA34(IELEM) * Y(IKLE3(IELEM))
             W5(IELEM) =   + XA15(IELEM) * Y(IKLE1(IELEM))
     *                     + XA25(IELEM) * Y(IKLE2(IELEM))
     *                     + XA35(IELEM) * Y(IKLE3(IELEM))
             W6(IELEM) =   + XA16(IELEM) * Y(IKLE1(IELEM))
     *                     + XA26(IELEM) * Y(IKLE2(IELEM))
     *                     + XA36(IELEM) * Y(IKLE3(IELEM))
           ENDDO
C
         ELSEIF(TYPEXT(1:1).EQ.'0') THEN
C
           CALL OV ('X=C     ', W1 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W2 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W3 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W4 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W5 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W6 , Y , Z , 0.D0 , NELEM )
C
         ELSE
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=YZ    ', X , Y , DA , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=Y     ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'0') THEN
           CALL OV ('X=C     ', X , Y , DA , 0.D0 , NPOIN )
         ELSE
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C        THE DIAGONAL REACHES ONLY LINEAR POINTS, OTHERS NOT INITIALIZED
C        THEY ARE SET TO ZERO HERE
         CALL OV ('X=C     ',X(NPOIN+1:NPT2),Y,Z,0.D0,NPT2-NPOIN)
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=-TAY  ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 101 IELEM = 1 , NELEM
             W1(IELEM) =   - XA21(IELEM) * Y(IKLE2(IELEM))
     *                     - XA31(IELEM) * Y(IKLE3(IELEM))
             W2(IELEM) =   - XA12(IELEM) * Y(IKLE1(IELEM))
     *                     - XA32(IELEM) * Y(IKLE3(IELEM))
             W3(IELEM) =   - XA13(IELEM) * Y(IKLE1(IELEM))
     *                     - XA23(IELEM) * Y(IKLE2(IELEM))
             W4(IELEM) =   - XA14(IELEM) * Y(IKLE1(IELEM))
     *                     - XA24(IELEM) * Y(IKLE2(IELEM))
     *                     - XA34(IELEM) * Y(IKLE3(IELEM))
             W5(IELEM) =   - XA15(IELEM) * Y(IKLE1(IELEM))
     *                     - XA25(IELEM) * Y(IKLE2(IELEM))
     *                     - XA35(IELEM) * Y(IKLE3(IELEM))
             W6(IELEM) =   - XA16(IELEM) * Y(IKLE1(IELEM))
     *                     - XA26(IELEM) * Y(IKLE2(IELEM))
     *                     - XA36(IELEM) * Y(IKLE3(IELEM))
101        CONTINUE
C
         ELSEIF(TYPEXT(1:1).EQ.'0') THEN
C
           CALL OV ('X=C     ', W1 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W2 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W3 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W4 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W5 , Y , Z , 0.D0 , NELEM )
           CALL OV ('X=C     ', W6 , Y , Z , 0.D0 , NELEM )
C
         ELSE
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(0)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=-YZ   ', X , Y , DA , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=-Y    ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'0') THEN
           CALL OV ('X=C     ', X , Y , DA , 0.D0 , NPOIN )
         ELSE
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C        THE DIAGONAL REACHES ONLY LINEAR POINTS, OTHERS NOT INITIALIZED
C        THEY ARE SET TO ZERO HERE
         CALL OV ('X=C     ',X(NPOIN+1:NPT2),Y,Z,0.D0,NPT2-NPOIN)
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+TAY ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 120 IELEM = 1 , NELEM
             W1(IELEM) = W1(IELEM) + XA21(IELEM) * Y(IKLE2(IELEM))
     *                             + XA31(IELEM) * Y(IKLE3(IELEM))
             W2(IELEM) = W2(IELEM) + XA12(IELEM) * Y(IKLE1(IELEM))
     *                             + XA32(IELEM) * Y(IKLE3(IELEM))
             W3(IELEM) = W3(IELEM) + XA13(IELEM) * Y(IKLE1(IELEM))
     *                             + XA23(IELEM) * Y(IKLE2(IELEM))
             W4(IELEM) = W4(IELEM) + XA14(IELEM) * Y(IKLE1(IELEM))
     *                             + XA24(IELEM) * Y(IKLE2(IELEM))
     *                             + XA34(IELEM) * Y(IKLE3(IELEM))
             W5(IELEM) = W5(IELEM) + XA15(IELEM) * Y(IKLE1(IELEM))
     *                             + XA25(IELEM) * Y(IKLE2(IELEM))
     *                             + XA35(IELEM) * Y(IKLE3(IELEM))
             W6(IELEM) = W6(IELEM) + XA16(IELEM) * Y(IKLE1(IELEM))
     *                             + XA26(IELEM) * Y(IKLE2(IELEM))
     *                             + XA36(IELEM) * Y(IKLE3(IELEM))
120        CONTINUE
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=X+YZ  ', X , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=X+Y   ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X-TAY ') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 140 IELEM = 1 , NELEM
             W1(IELEM) = W1(IELEM) - XA21(IELEM) * Y(IKLE2(IELEM))
     *                             - XA31(IELEM) * Y(IKLE3(IELEM))
             W2(IELEM) = W2(IELEM) - XA12(IELEM) * Y(IKLE1(IELEM))
     *                             - XA32(IELEM) * Y(IKLE3(IELEM))
             W3(IELEM) = W3(IELEM) - XA13(IELEM) * Y(IKLE1(IELEM))
     *                             - XA23(IELEM) * Y(IKLE2(IELEM))
             W4(IELEM) = W4(IELEM) - XA14(IELEM) * Y(IKLE1(IELEM))
     *                             - XA24(IELEM) * Y(IKLE2(IELEM))
     *                             - XA34(IELEM) * Y(IKLE3(IELEM))
             W5(IELEM) = W5(IELEM) - XA15(IELEM) * Y(IKLE1(IELEM))
     *                             - XA25(IELEM) * Y(IKLE2(IELEM))
     *                             - XA35(IELEM) * Y(IKLE3(IELEM))
             W6(IELEM) = W6(IELEM) - XA16(IELEM) * Y(IKLE1(IELEM))
     *                             - XA26(IELEM) * Y(IKLE2(IELEM))
     *                             - XA36(IELEM) * Y(IKLE3(IELEM))
140        CONTINUE
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=X-YZ  ', X , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=X-Y   ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+CTAY') THEN
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q') THEN
C
           DO 160 IELEM = 1 , NELEM
             W1(IELEM) = W1(IELEM)
     *                 + C * (    + XA21(IELEM) * Y(IKLE2(IELEM))
     *                            + XA31(IELEM) * Y(IKLE3(IELEM)) )
             W2(IELEM) = W2(IELEM)
     *                 + C * (    + XA12(IELEM) * Y(IKLE1(IELEM))
     *                            + XA32(IELEM) * Y(IKLE3(IELEM)) )
             W3(IELEM) = W3(IELEM)
     *                 + C * (    + XA13(IELEM) * Y(IKLE1(IELEM))
     *                            + XA23(IELEM) * Y(IKLE2(IELEM)) )
             W4(IELEM) = W4(IELEM)
     *                 + C * (    + XA14(IELEM) * Y(IKLE1(IELEM))
     *                            + XA24(IELEM) * Y(IKLE2(IELEM))
     *                            + XA34(IELEM) * Y(IKLE3(IELEM)) )
             W5(IELEM) = W5(IELEM)
     *                 + C * (    + XA15(IELEM) * Y(IKLE1(IELEM))
     *                            + XA25(IELEM) * Y(IKLE2(IELEM))
     *                            + XA35(IELEM) * Y(IKLE3(IELEM)) )
             W6(IELEM) = W6(IELEM)
     *                 + C * (    + XA16(IELEM) * Y(IKLE1(IELEM))
     *                            + XA26(IELEM) * Y(IKLE2(IELEM))
     *                            + XA36(IELEM) * Y(IKLE3(IELEM)) )
160        CONTINUE
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF (LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF (LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C   CONTRIBUTION DE LA DIAGONALE
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=X+CYZ ', X , Y , DA , C , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=X+CY  ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).NE.'0') THEN
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,3000) OP
        IF (LNG.EQ.2) WRITE(LU,3001) OP
        CALL PLANTE(1)
        STOP
C
C-----------------------------------------------------------------------
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
C
1000  FORMAT(1X,'MV0306 (BIEF) : TERMES EXTRADIAG. TYPE INCONNU: ',A1)
1001  FORMAT(1X,'MV0306 (BIEF) : EXTRADIAG. TERMS  UNKNOWN TYPE : ',A1)
2000  FORMAT(1X,'MV0306 (BIEF) : DIAGONALE : TYPE INCONNU: ',A1)
2001  FORMAT(1X,'MV0306 (BIEF) : DIAGONAL : UNKNOWN TYPE : ',A1)
3000  FORMAT(1X,'MV0306 (BIEF) : OPERATION INCONNUE : ',A8)
3001  FORMAT(1X,'MV0306 (BIEF) : UNKNOWN OPERATION : ',A8)
C
C-----------------------------------------------------------------------
C
      END