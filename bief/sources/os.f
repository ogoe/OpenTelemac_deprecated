C                       *************
                        SUBROUTINE OS
C                       *************
C
     * ( OP , X , Y , Z , C , IOPT , INFINI , ZERO )
C
C***********************************************************************
C BIEF VERSION 5.6        18/08/05    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : OPERATIONS SUR LES STRUCTURES
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET Z ET LA CONSTANTE C. LE RESULTAT
C   EST LE VECTEUR X.
C
C   SUR LES TABLEAUX QUELCONQUES OU LES VECTEURS :
C
C   OP = 'X=C     '     :  X MIS A LA VALEUR C
C   OP = 'X=0     '     :  X MIS A 0
C   OP = 'X=Y     '     :  Y COPIE DANS X
C   OP = 'X=+Y    '     :  IDEM
C   OP = 'X=-Y    '     : -Y COPIE DANS X
C   OP = 'X=1/Y   '     :  INVERSE DE Y MIS DANS X
C   OP = 'X=Y+Z   '     :  SOMME DE Y ET Z MISE DANS X
C   OP = 'X=Y-Z   '     :  DIFFERENCE DE Y ET Z MISE DANS X
C   OP = 'X=YZ    '     :  PRODUIT Y PAR  Z MIS DANS X
C   OP = 'X=-YZ   '     :  PRODUIT Y PAR  Z MIS DANS X
C   OP = 'X=XY    '     :  PRODUIT Y PAR  X MIS DANS X
C   OP = 'X=X+YZ  '     :  PRODUIT Y PAR  Z AJOUTE A X
C   OP = 'X=X-YZ  '     :  PRODUIT Y PAR  Z RETRANCHE A X
C   OP = 'X=CXY   '     :  PRODUIT DE C X ET Y MIS DANS X
C   OP = 'X=CYZ   '     :  PRODUIT DE C, Y ET Z MIS DANS X
C   OP = 'X=CXYZ  '     :  PRODUIT DE C, X, Y ET Z MIS DANS X
C   OP = 'X=X+CYZ '     :  PRODUIT DE C, Y ET Z AJOUTE A X
C   OP = 'X=Y/Z   '     :  DIVISION DE Y PAR Z MIS DANS X
C   OP = 'X=CY/Z  '     :  PRODUIT DE C ET Y DIVISE PAR Z ET MIS DANS X
C   OP = 'X=CXY/Z '     :  PRODUIT DE C, X ET Y DIVISE PAR Z ET MIS DS X
C   OP = 'X=X+CY/Z'     :  PRODUIT DE C ET Y DIVISE PAR Z AJOUTE A X
C   OP = 'X=X+Y   '     :  Y AJOUTE A X
C   OP = 'X=X-Y   '     :  Y RETRANCHE A X
C   OP = 'X=CX    '     :  X MULTIPLIE PAR C
C   OP = 'X=CY    '     :  CY MIS DANS X
C   OP = 'X=Y+CZ  '     :  CZ AJOUTE A Y ET MIS DANS X
C   OP = 'X=X+CY  '     :  CY AJOUTE A X
C   OP = 'X=SQR(Y)'     :  RACINE DE Y MIS DANS X
C   OP = 'X=ABS(Y)'     :  VALEUR ABSOLUE DE Y MIS DANS X
C   OP = 'X=N(Y,Z)'     :  X MODULE DU VECTEUR DE COMPOSANTES Y ET Z
C   OP = 'X=Y+C   '     :  C AJOUTE A Y MIS DANS X
C   OP = 'X=X+C   '     :  C AJOUTE A X
C   OP = 'X=Y**C  '     :  Y A LA PUISSANCE C MIS DANS X
C   OP = 'X=COS(Y)'     :  COSINUS DE Y MIS DANS X
C   OP = 'X=SIN(Y)'     :  SINUS DE Y MIS DANS X
C   OP = 'X=ATN(Y)'     :  ARC TANGENTE DE Y MIS DANS X
C   OP = 'X=A(Y,Z)'     :  ATAN2(Y,Z) MIS DANS X
C   OP = 'X=+(Y,C)'     :  X = MAX DE Y ET DE C
C   OP = 'X=-(Y,C)'     :  X = MIN DE Y ET DE C
C   OP = 'X=+(Y,Z)'     :  X = MAX DE Y ET DE Z
C   OP = 'X=-(Y,Z)'     :  X = MIN DE Y ET DE Z
C   OP = 'X=YIFZ<C'     :  X = Y SI Z < C (POUR CHAQUE POINT CONSIDERE)
C   OP = 'X=C(Y-Z)'     :  X = C*(Y-Z)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |      OP        | -->| CHAINE DE CARACTERES INDIQUANT L'OPERATION
C |                |   >| A EFFECTUER.
C |      X         |<-- | STRUCTURE RESULTAT
C |      Y         | -->| STRUCTURE OPERANDE
C |      Z         | -->| STRUCTURE OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELANTS : BEAUCOUP
C PROGRAMMES APPELES   : NEANT
C
C PRECAUTIONS D'EMPLOI : LES OPERATIONS 1/Y ET Y/Z
C                        SUPPRIMENT LES DIVISIONS PAR 0.
C                        UN PASSAGE CORRECT N'EST DONC PAS
C                        UNE PREUVE QUE Y OU Z N'EST JAMAIS NUL.
C
C**********************************************************************
C
      USE BIEF, EX_OS => OS
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     OPTIONAL ARGUMENTS
C
      INTEGER,          INTENT(IN), OPTIONAL :: IOPT
      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: INFINI
      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: ZERO
C
C     ARGUMENTS
C
      TYPE(BIEF_OBJ),   INTENT(INOUT), OPTIONAL, TARGET :: X
      TYPE(BIEF_OBJ),   INTENT(IN)   , OPTIONAL, TARGET :: Y,Z
      DOUBLE PRECISION, INTENT(IN)   , OPTIONAL :: C
      CHARACTER(LEN=8), INTENT(IN)              :: OP
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C     LOCAL VARIABLES
C
      INTEGER IBL,TYPX,IDIM,N,NMAX
      LOGICAL YAY,YAZ,YAC
      TYPE(BIEF_OBJ), POINTER :: YY,ZZ
      DOUBLE PRECISION CC
C
C-----------------------------------------------------------------------
C
      TYPX = X%TYPE
C
      YAY=.FALSE.
      YAZ=.FALSE.
      YAC=.FALSE.
      IF(OP(3:3).EQ.'Y'.OR.OP(4:4).EQ.'Y'.OR.OP(5:5).EQ.'Y'.OR.
     *   OP(6:6).EQ.'Y'.OR.OP(7:7).EQ.'Y'.OR.OP(8:8).EQ.'Y') YAY=.TRUE.
      IF(OP(3:3).EQ.'Z'.OR.OP(4:4).EQ.'Z'.OR.OP(5:5).EQ.'Z'.OR.
     *   OP(6:6).EQ.'Z'.OR.OP(7:7).EQ.'Z'.OR.OP(8:8).EQ.'Z') YAZ=.TRUE.
C     
C     CHECKING THAT CONSTANT C IS IN THE REQUIRED OPERATION
C     I.E. IF THERE IS C IN OP, EXCEPT WHEN IT IS X=COS(Y)
C
      IF((OP(3:3).EQ.'C'.AND.OP(4:4).NE.'O').OR.
     *    OP(4:4).EQ.'C'.OR.OP(5:5).EQ.'C'.OR.
     *    OP(6:6).EQ.'C'.OR.OP(7:7).EQ.'C'.OR.OP(8:8).EQ.'C') YAC=.TRUE.
C
      IF(PRESENT(C)) THEN
        CC=C
      ELSE
        IF(YAC) THEN
          IF (LNG.EQ.1) WRITE(LU,1) OP
          IF (LNG.EQ.2) WRITE(LU,2) OP
1         FORMAT(1X,'OS (BIEF) : C ABSENT ET OPERATION ',A8,' DEMANDEE')
2         FORMAT(1X,'OS (BIEF) : C MISSING AND OPERATION ',A8,' ASKED')
          CALL PLANTE(1)
          STOP
        ENDIF
      ENDIF
C
      IF(YAY) THEN
        IF(PRESENT(Y)) THEN
          YY=>Y
        ELSE
          IF (LNG.EQ.1) WRITE(LU,10) OP
          IF (LNG.EQ.2) WRITE(LU,11) OP
10        FORMAT(1X,'OS (BIEF) : Y ABSENT ET OPERATION ',A8,' DEMANDEE')
11        FORMAT(1X,'OS (BIEF) : Y MISSING AND OPERATION ',A8,' ASKED')
          CALL PLANTE(1)
          STOP
        ENDIF
      ELSE
        YY=>X
      ENDIF
C
C     OPERATION WITH Y AND Z (IF THERE IS Z THERE SHOULD BE Y)
C
      IF(YAZ) THEN
C
        IF(PRESENT(Z)) THEN
C
        ZZ=>Z
C
C       COMPARING TYPES OF Y AND Z
C
        IF(.NOT.CMPOBJ(Y,Z)) THEN
          IF (LNG.EQ.1) WRITE(LU,40) Y%NAME,Y%ELM,Z%NAME,Z%ELM
          IF (LNG.EQ.2) WRITE(LU,41) Y%NAME,Y%ELM,Z%NAME,Z%ELM
40        FORMAT(1X,'OS (BIEF) : TYPES DIFFERENTS POUR ',A6,' (',1I2,
     *              ') ET ',A6,' (',1I2,')')
41        FORMAT(1X,'OS (BIEF) : DIFFERENT TYPES FOR ',A6,' (',1I2,
     *              ') AND ',A6,' (',1I2,')')
          CALL PLANTE(1)
          STOP
        ENDIF
C
        ELSE
C
          IF (LNG.EQ.1) WRITE(LU,20) OP
          IF (LNG.EQ.2) WRITE(LU,21) OP
20        FORMAT(1X,'OS (BIEF) : Z ABSENT ET OPERATION ',A8,' DEMANDEE')
21        FORMAT(1X,'OS (BIEF) : Z MISSING AND OPERATION ',A8,' ASKED')
          CALL PLANTE(1)
          STOP
C
        ENDIF
C
      ELSE
        ZZ=>X      
      ENDIF
C
C-----------------------------------------------------------------------
C     VECTORS
C-----------------------------------------------------------------------
C
      IF(TYPX.EQ.2) THEN
C
C     OPERATION WITH Y : Y IS CHECKED
C
        IF(YAY) THEN
C         TYPES DIFFERENTS : X PREND ALORS LA STRUCTURE DE Y.
          IF(.NOT.CMPOBJ(X,Y)) CALL CPSTVC(Y,X)
        ENDIF
C
C       CHECKING MEMORY
C
        IF(X%DIM1.GT.X%MAXDIM1) THEN
          IF (LNG.EQ.1) WRITE(LU,100) X%NAME
          IF (LNG.EQ.2) WRITE(LU,101) X%NAME
100       FORMAT(1X,'OS (BIEF) : DEPASSEMENT DE MEMOIRE SUR : ',A6)
101       FORMAT(1X,'OS (BIEF) : BEYOND ALLOWED MEMORY IN: ',A6)
          CALL PLANTE(1)
          STOP
        ENDIF
C
        IF(.NOT.PRESENT(IOPT)) THEN
C
        IF(X%DIM2.GT.1) THEN
C
          DO IDIM = 1 , X%DIM2
            CALL OV_2(OP,X%R,IDIM,YY%R,IDIM,
     *                            ZZ%R,IDIM,CC,X%MAXDIM1,X%DIM1)
          END DO
C
        ELSE
C
          CALL OV(OP,X%R,YY%R,ZZ%R,CC,X%DIM1)
C
        ENDIF
C
        ELSE
C
        IF(X%DIM2.GT.1) THEN
C
          DO IDIM = 1 , X%DIM2
            CALL OVD_2(OP,X%R,IDIM,YY%R,IDIM,ZZ%R,IDIM,CC,
     *                 X%MAXDIM1,X%DIM1,IOPT,INFINI,ZERO)
          END DO
C
        ELSE
C
          CALL OVD(OP,X%R,YY%R,ZZ%R,CC,X%DIM1,IOPT,INFINI,ZERO)
C
        ENDIF
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(TYPX.EQ.4) THEN
C
C-----------------------------------------------------------------------
C     BLOCKS
C-----------------------------------------------------------------------
C
       DO 60 IBL = 1 , X%N
         IF(YAY) THEN
           IF(.NOT.CMPOBJ(X%ADR(IBL)%P,Y%ADR(IBL)%P)) THEN
             CALL CPSTVC(Y%ADR(IBL)%P,X%ADR(IBL)%P)
           ENDIF
         ENDIF
C
C        CHECKING MEMORY
C
         N = X%ADR(IBL)%P%DIM1
         NMAX = X%ADR(IBL)%P%MAXDIM1
         IF(N.GT.NMAX) THEN
           IF (LNG.EQ.1) WRITE(LU,100) X%ADR(IBL)%P%NAME
           IF (LNG.EQ.2) WRITE(LU,101) X%ADR(IBL)%P%NAME
           IF (LNG.EQ.1) WRITE(LU,200) X%NAME
           IF (LNG.EQ.2) WRITE(LU,201) X%NAME
200        FORMAT(1X,'            CE VECTEUR EST DANS LE BLOC : ',A6)
201        FORMAT(1X,'            THIS VECTOR IS IN BLOCK: ',A6)
           CALL PLANTE(1)
           STOP
         ENDIF
C
         IF(.NOT.PRESENT(IOPT)) THEN
C
         IF(X%ADR(IBL)%P%DIM2.GT.1) THEN
C
         DO IDIM = 1 , X%ADR(IBL)%P%DIM2
           CALL OV_2(OP,X%ADR(IBL)%P%R,IDIM,
     *                 YY%ADR(IBL)%P%R,IDIM,
     *                 ZZ%ADR(IBL)%P%R,IDIM, CC , NMAX , N )
         END DO
C
         ELSE
C
           CALL OV(OP,X%ADR(IBL)%P%R,
     *               YY%ADR(IBL)%P%R,
     *               ZZ%ADR(IBL)%P%R, CC , N )
C
         ENDIF
C
         ELSE
C
         IF(X%ADR(IBL)%P%DIM2.GT.1) THEN
C
         DO IDIM = 1 , X%ADR(IBL)%P%DIM2
           CALL OVD_2(OP,X%ADR(IBL)%P%R,IDIM,
     *                  YY%ADR(IBL)%P%R,IDIM,
     *                  ZZ%ADR(IBL)%P%R,IDIM, CC , NMAX , N ,
     *                  IOPT,INFINI,ZERO)
         END DO
C
         ELSE
C
           CALL OVD(OP,X%ADR(IBL)%P%R,
     *                YY%ADR(IBL)%P%R,
     *                ZZ%ADR(IBL)%P%R, CC , N ,IOPT,INFINI,ZERO)
C
         ENDIF
C
         ENDIF
C
C
60     CONTINUE
C
C-----------------------------------------------------------------------
C
C     ERREUR OU OBJET NON TRAITE
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,1000) X%TYPE
        IF (LNG.EQ.2) WRITE(LU,1001) X%TYPE
1000    FORMAT(1X,'OS (BIEF) : TYPE D''OBJET NON TRAITE: ',1I3)
1001    FORMAT(1X,'OS (BIEF) : TYPE OF OBJECT NOT IMPLEMENTED: ',1I3)
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END       
