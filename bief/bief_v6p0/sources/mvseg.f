C                       ****************
                        SUBROUTINE MVSEG
C                       ****************
C
     *(OP, X , DA,TYPDIA,XA1,XA2,
     * TYPEXT, Y,C,NPOIN,NELEM,NSEG1,NSEG2,GLOSEG1,GLOSEG2,IELM1,IELM2)
C
C***********************************************************************
C BIEF VERSION 5.1           05/02/91    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C  FONCTION : MATRIX VECTOR PRODUCT FOR EDGE-BASED STORAGE
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
C      OP = 'X=CAY   '  : X = CAY
C      OP = 'X=-AY   '  : X = -AY
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
C PROGRAMMES APPELES   : ASSVEC , OV , PLANTE
C
C***********************************************************************
C
      USE BIEF,EX_MVSEG => MVSEG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN
      INTEGER, INTENT(IN) :: GLOSEG1(*),GLOSEG2(*)
      INTEGER, INTENT(IN) :: NSEG1,NSEG2,NELEM,IELM1,IELM2
C
      DOUBLE PRECISION, INTENT(IN)    :: Y(*),DA(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      DOUBLE PRECISION, INTENT(IN)    :: XA1(*),XA2(*)
      DOUBLE PRECISION, INTENT(IN)    :: C
C
      CHARACTER(LEN=8),INTENT(IN) :: OP
      CHARACTER(LEN=1),INTENT(IN) :: TYPDIA,TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISEG,I,MINSEG,MAXSEG
      DOUBLE PRECISION Z(1)
C
      INTRINSIC MIN,MAX
C
C-----------------------------------------------------------------------
C
      MINSEG = MIN(NSEG1,NSEG2)
      MAXSEG = MAX(NSEG1,NSEG2)
C
      IF(OP(1:8).EQ.'X=AY    ') THEN
C
C   CONTRIBUTION DE LA DIAGONALE :
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
C          SQUARE PART
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+XA1(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+XA2(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO 
C
C          THE REST OF THE RECTANGULAR MATRIX
C      
           IF(NSEG1.GT.NSEG2) THEN
C            UNE PARTIE DE X N'A PAS ETE INITIALISEE
C            PAR LA CONTRIBUTION DE LA DIAGONALE
             IF(IELM1.EQ.12.OR.IELM2.EQ.12) THEN
C              OPTIMISATION POUR ELEMENT QUASI-BULLE
               DO I = NPOIN+1,NPOIN+NELEM
                 X(I)=0.D0
               ENDDO
             ELSE
               DO ISEG = MINSEG+1,MAXSEG
                 X(GLOSEG2(ISEG))=0.D0
               ENDDO
             ENDIF
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ENDIF
C
         ELSEIF(TYPEXT(1:1).NE.'0') THEN
C
           IF(LNG.EQ.1) WRITE(LU,1000) TYPEXT
           IF(LNG.EQ.2) WRITE(LU,1001) TYPEXT
           CALL PLANTE(1)
           STOP
C
         ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=CAY   ') THEN
C
C   CONTRIBUTION DE LA DIAGONALE :
C
         IF(TYPDIA(1:1).EQ.'Q') THEN
           CALL OV ('X=CYZ   ', X , Y , DA , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'I') THEN
           CALL OV ('X=CY    ', X , Y , Z  , C  , NPOIN )
         ELSEIF(TYPDIA(1:1).EQ.'0') THEN
           CALL OV ('X=C     ', X , Y , DA , 0.D0 , NPOIN )
         ELSE
           IF (LNG.EQ.1) WRITE(LU,2000) TYPDIA
           IF (LNG.EQ.2) WRITE(LU,2001) TYPDIA
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+C*XA1(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+C*XA2(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             IF(IELM1.EQ.12.OR.IELM2.EQ.12) THEN
C              OPTIMISATION POUR ELEMENT QUASI-BULLE
               DO I = NPOIN+1,NPOIN+NELEM
                 X(I)=0.D0
               ENDDO
             ELSE
               DO ISEG = MINSEG+1,MAXSEG
                 X(GLOSEG2(ISEG))=0.D0
               ENDDO
             ENDIF
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+C*XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+C*XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=-AY   ') THEN
C
C   CONTRIBUTION DE LA DIAGONALE :
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))-XA1(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))-XA2(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             IF(IELM1.EQ.12.OR.IELM2.EQ.12) THEN
C              OPTIMISATION POUR ELEMENT QUASI-BULLE
               DO I = NPOIN+1,NPOIN+NELEM
                 X(I)=0.D0
               ENDDO
             ELSE
               DO ISEG = MINSEG+1,MAXSEG
                 X(GLOSEG2(ISEG))=0.D0
               ENDDO
             ENDIF
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))-XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))-XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+AY  ') THEN
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+XA1(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+XA2(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X-AY  ') THEN
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
           CALL PLANTE(1)
           STOP
         ENDIF
C
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))-XA1(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))-XA2(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))-XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))-XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+CAY ') THEN
C
C   CONTRIBUTION DE LA DIAGONALE :
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+C*XA1(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+C*XA2(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+C*XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+C*XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=TAY   ') THEN
C
C   CONTRIBUTION DE LA DIAGONALE :
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+XA2(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+XA1(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             IF(IELM1.EQ.12.OR.IELM2.EQ.12) THEN
C              OPTIMISATION POUR ELEMENT QUASI-BULLE
               DO I = NPOIN+1,NPOIN+NELEM
                 X(I)=0.D0
               ENDDO
             ELSE
               DO ISEG = MINSEG+1,MAXSEG
                 X(GLOSEG2(ISEG))=0.D0
               ENDDO
             ENDIF
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=-TAY   ') THEN
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))-XA2(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))-XA1(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))-XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             IF(IELM1.EQ.12.OR.IELM2.EQ.12) THEN
C              OPTIMISATION POUR ELEMENT QUASI-BULLE
               DO I = NPOIN+1,NPOIN+NELEM
                 X(I)=0.D0
               ENDDO
             ELSE
               DO ISEG = MINSEG+1,MAXSEG
                 X(GLOSEG2(ISEG))=0.D0
               ENDDO
             ENDIF
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))-XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+TAY ') THEN
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+XA2(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+XA1(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X-TAY ') THEN
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))-XA2(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))-XA1(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))-XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))-XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ENDIF
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
C-----------------------------------------------------------------------
C
      ELSEIF(OP(1:8).EQ.'X=X+CTAY') THEN
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
C   CONTRIBUTION DES TERMES EXTRADIAGONAUX :
C
         IF(TYPEXT(1:1).EQ.'Q'.OR.TYPEXT(1:1).EQ.'S') THEN
C
           DO ISEG = 1 , MINSEG
             X(GLOSEG1(ISEG))=
     *       X(GLOSEG1(ISEG))+C*XA2(ISEG)*Y(GLOSEG2(ISEG))
             X(GLOSEG2(ISEG))=
     *       X(GLOSEG2(ISEG))+C*XA1(ISEG)*Y(GLOSEG1(ISEG))
           ENDDO
           IF(NSEG1.GT.NSEG2) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG1(ISEG))=
     *         X(GLOSEG1(ISEG))+C*XA2(ISEG)*Y(GLOSEG2(ISEG))
             ENDDO
           ELSEIF(NSEG2.GT.NSEG1) THEN
             DO ISEG = MINSEG+1,MAXSEG
               X(GLOSEG2(ISEG))=
     *         X(GLOSEG2(ISEG))+C*XA2(ISEG)*Y(GLOSEG1(ISEG))
             ENDDO
           ENDIF
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
1000  FORMAT(1X,'MVSEG (BIEF) : TERMES EXTRADIAG. TYPE INCONNU: ',A1)
1001  FORMAT(1X,'MVSEG (BIEF) : EXTRADIAG. TERMS  UNKNOWN TYPE : ',A1)
2000  FORMAT(1X,'MVSEG (BIEF) : DIAGONALE : TYPE INCONNU: ',A1)
2001  FORMAT(1X,'MVSEG (BIEF) : DIAGONAL : UNKNOWN TYPE : ',A1)
3000  FORMAT(1X,'MVSEG (BIEF) : OPERATION INCONNUE : ',A8)
3001  FORMAT(1X,'MVSEG (BIEF) : UNKNOWN OPERATION : ',A8)
C
      END
