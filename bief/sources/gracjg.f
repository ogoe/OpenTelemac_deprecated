C                       *****************
                        SUBROUTINE GRACJG
C                       *****************
C
     *(X, A,B , MESH, D,AD,G,R, CFG,INFOGR,AUX)
C
C***********************************************************************
C BIEF VERSION 5.9         06/10/08    J-M HERVOUET (LNH) 01 30 87 80 18
C***********************************************************************
C
C  FONCTION : RESOLUTION D'UN SYSTEME LINEAIRE A X = B
C             PAR LA METHODE DU GRADIENT CONJUGUE.
C
C-----------------------------------------------------------------------
C                        PRECONDITIONNEMENT
C-----------------------------------------------------------------------
C VALEUR DE CFG%PRECON                      SIGNIFICATION
C-----------------------------------------------------------------------
C                     I
C        0 OU 1       I  RIEN
C                     I
C        2            I  PRECONDITIONNEMENT DIAGONAL AVEC LA DIAGONALE
C                     I  DE LA MATRICE.
C        3            I  PRECONDITIONNEMENT BLOC-DIAGONAL
C                     I
C        5            I  PRECONDITIONNEMENT DIAGONAL AVEC LA VALEUR
C                     I  ABSOLUE DE LA DIAGONALE DE LA MATRICE.
C                     I
C        7            I  CROUT EBE
C                     I
C       11            I  GAUSS-SEIDEL EBE
C                     I
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  VALEUR INITIALE, PUIS SOLUTION
C |      A         | -->|  MATRICE DU SYSTEME
C |      B         | -->|  SECOND MEMBRE DU SYSTEME.
C |      MESH      | -->|  BLOC DES ENTIERS DU MAILLAGE.
C |      D         |<-->|  DIRECTION DE DESCENTE.
C |      AD        |<-->|  MATRICE A MULTIPLIEE PAR D.
C |      G         |<-->|  GRADIENT DE DESCENTE.
C |      R         |<-->|  RESIDU (CONFONDU AVEC LE GRADIENT SI IL N'Y A
C |                |    |  PAS DE PRECONDITIONNEMENT DANS GRACJG)
C |      INFOGR    | -->|  SI OUI, IMPRESSION D'UN COMPTE-RENDU.
C |      AUX       | -->|  MATRICE POUR LE PRECONDITIONNEMENT.
C |________________|____|______________________________________________-
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : MATRBL , OS , DOTS
C
C**********************************************************************
C
      USE BIEF, EX_GRACJG => GRACJG
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(SLVCFG), INTENT(INOUT)    :: CFG
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: B
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: D,AD,G,R,X
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(IN)     :: A
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: AUX
      LOGICAL, INTENT(IN)            :: INFOGR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER M
C
      DOUBLE PRECISION XL,RMRM,RMDM,RMGM,TESTL
      DOUBLE PRECISION BETA,RO,DAD,RM1GM1,STO1,C
C
      LOGICAL RELAT,PREC,CROUT,GSEB,PRE3D,PREBE
C
      DOUBLE PRECISION RMIN
      DATA RMIN/1.D-15/
C
C-----------------------------------------------------------------------
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
C   INITIALISATIONS
C
      STO1 = 0.D0
      CROUT =.FALSE.
      IF(07*(CFG%PRECON/07).EQ.CFG%PRECON) CROUT=.TRUE.
      GSEB=.FALSE.
      IF(11*(CFG%PRECON/11).EQ.CFG%PRECON) GSEB=.TRUE.
      PREBE=.FALSE.
      IF(13*(CFG%PRECON/13).EQ.CFG%PRECON) PREBE=.TRUE.
      PRE3D=.FALSE.
      IF(17*(CFG%PRECON/17).EQ.CFG%PRECON) PRE3D=.TRUE.
      PREC=.FALSE.
      IF(CROUT.OR.GSEB.OR.PREBE.OR.PRE3D) PREC=.TRUE.
C
C-----------------------------------------------------------------------
C   INITIALISATIONS
C-----------------------------------------------------------------------
C
      M   = 0
C
C  NORME DU SECOND MEMBRE POUR CALCULER LA PRECISION RELATIVE :
C
      XL = P_DOTS(B,B,MESH)
C
C     TEST OLIVIER BOITEAU (SINETICS) : SECOND MEMBRE TROP PETIT
C                                       ==> INCONNUE = 0
C
      TESTL=SQRT(XL)
      IF(TESTL.LT.RMIN) THEN        
        CALL OS('X=0     ',X=X)
        IF(INFOGR) THEN
          IF(LNG.EQ.1) WRITE(LU,105) TESTL
          IF(LNG.EQ.2) WRITE(LU,106) TESTL
        ENDIF
        GOTO 1000
      ENDIF
C
      IF(XL.LT.1.D0) THEN
         XL = 1.D0
         RELAT = .FALSE.
      ELSE
         RELAT = .TRUE.
      ENDIF
C
C CALCUL DU RESIDU INITIAL ET SORTIE EVENTUELLE :
C
      CALL MATRBL( 'X=AY    ',R,A,X,  C,MESH)
C
      CALL OS( 'X=X-Y   ' , R , B , B , C )
      RMRM   = P_DOTS(R,R,MESH)
      RMGM   = RMRM
C
      IF (RMRM.LT.CFG%EPS**2*XL) GO TO 900
C
C-----------------------------------------------------------------------
C PRECONDITIONNEMENT :
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
C       CALCUL DE C G0 = R0
        IF(CROUT.OR.GSEB.OR.PREBE) THEN
          CALL DOWNUP( G, AUX , R , 'D' , MESH )
          IF(NCSIZE.GT.1) CALL PARMOY(G,MESH)
        ELSEIF(PRE3D) THEN
          CALL CPSTVC(R%ADR(1)%P,G%ADR(1)%P)         
          CALL TRID3D(AUX%X%R,G%ADR(1)%P%R,R%ADR(1)%P%R,
     *                MESH%NPOIN,NBPTS(11))
        ENDIF        
C       CALCUL DE RMGM ET STOCKAGE
        RMGM = P_DOTS(R,G,MESH)
        STO1 = RMGM
      ENDIF
C
C-----------------------------------------------------------------------
C CALCUL DE LA DIRECTION DE DESCENTE INITIALE :
C-----------------------------------------------------------------------
C
      CALL OS( 'X=Y     ' , D , G , G , C )
C
C-----------------------------------------------------------------------
C CALCUL DU PRODUIT A D INITIAL :
C-----------------------------------------------------------------------
C
      CALL MATRBL( 'X=AY    ',AD,A,D,C,MESH)
C
C-----------------------------------------------------------------------
C CALCUL DE RO INITIAL :
C-----------------------------------------------------------------------
C
      DAD = P_DOTS(D,AD,MESH)
      RO = RMGM / DAD
C     RO = RMDM / DAD  (ON MET RMGM CAR ICI D0=G0)
C
C-----------------------------------------------------------------------
C
C CALCUL DE X1 = X0 - RO  * D
C
      CALL OS('X=X+CY  ',X=X,Y=D,C=-RO)
C
C-----------------------------------------------------------------------
C  BOUCLE DES ITERATIONS :
C-----------------------------------------------------------------------
C
2     M  = M  + 1
C
C-----------------------------------------------------------------------
C CALCUL DU RESIDU : R(M) = R(M-1) - RO(M-1) A D(M-1)
C-----------------------------------------------------------------------
C
      CALL OS('X=X+CY  ',X=R,Y=AD,C=-RO)
C
C  CERTAINES VALEURS SERONT CHANGEES EN CAS DE PRECONDITIONNEMENT
C
      RM1GM1 = RMGM
      RMRM   = P_DOTS(R,R,MESH)
      RMDM   = RMRM
      RMGM   = RMRM
C
C TEST DE FIN :
C
      IF(RMRM.LE.XL*CFG%EPS**2) GO TO 900
C
C-----------------------------------------------------------------------
C PRECONDITIONNEMENT : RESOLUTION DE C G = R
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
C
C       RESOLUTION DE C G = R
        IF(CROUT.OR.GSEB.OR.PREBE) THEN
          CALL DOWNUP( G, AUX , R , 'D' , MESH )
          IF(NCSIZE.GT.1) CALL PARMOY(G,MESH)
        ELSEIF(PRE3D) THEN
          CALL CPSTVC(R%ADR(1)%P,G%ADR(1)%P)
          CALL TRID3D(AUX%X%R,G%ADR(1)%P%R,R%ADR(1)%P%R,
     *                MESH%NPOIN,NBPTS(11))
        ENDIF        
C       CALCUL DE RMGM ET RM1GM1
        RM1GM1 = STO1
        RMGM = P_DOTS(R,G,MESH)
        STO1=RMGM
C
      ENDIF
C
C-----------------------------------------------------------------------
C CALCUL DE D PAR RECURRENCE :
C-----------------------------------------------------------------------
C
      BETA = RMGM / RM1GM1
!     CALL OS( 'X=CX    ' , D , D , D , BETA )
!     CALL OS( 'X=X+Y   ' , D , G , G , C    )
!     EQUIVALENT OPTIMISE DES DEUX CALL OS :
      CALL OS( 'X=Y+CZ  ' , X=D , Y=G , Z=D , C=BETA )
C
C-----------------------------------------------------------------------
C PRECONDITIONNEMENT :
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
C       CALCUL DE RMDM
        RMDM=P_DOTS(R,D,MESH)
      ENDIF
C
C-----------------------------------------------------------------------
C CALCUL DE A D :
C-----------------------------------------------------------------------
C
      CALL MATRBL( 'X=AY    ',AD,A,D,C,MESH)
C
C-----------------------------------------------------------------------
C CALCUL DE RO
C-----------------------------------------------------------------------
C
      DAD = P_DOTS(D,AD,MESH)
      RO = RMDM/DAD
C
C CALCUL DE X(M) = X(M-1) - RO * D
C
      CALL OS('X=X+CY  ',X=X,Y=D,C=-RO)
C
      IF(M.LT.CFG%NITMAX) GO TO 2
C
C-----------------------------------------------------------------------
C
C     IF(INFOGR) THEN
        TESTL = SQRT( RMRM / XL )
        IF (RELAT) THEN
           IF (LNG.EQ.1) WRITE(LU,103) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,104) M,TESTL
        ELSE
           IF (LNG.EQ.1) WRITE(LU,203) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,204) M,TESTL
        ENDIF
C     ENDIF
      GO TO 1000
C
C-----------------------------------------------------------------------
C
900   CONTINUE
C
      IF(INFOGR) THEN
        TESTL = SQRT( RMRM / XL )
        IF (RELAT) THEN
           IF (LNG.EQ.1) WRITE(LU,101) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,102) M,TESTL
        ELSE
           IF (LNG.EQ.1) WRITE(LU,201) M,TESTL
           IF (LNG.EQ.2) WRITE(LU,202) M,TESTL
        ENDIF
      ENDIF
C
1000  CONTINUE
      RETURN
C
C-----------------------------------------------------------------------
C
C   FORMATS
C
101   FORMAT(1X,'GRACJG (BIEF) : ',
     *                     1I8,' ITERATIONS, PRECISION RELATIVE:',G16.7)
102   FORMAT(1X,'GRACJG (BIEF) : ',
     *                     1I8,' ITERATIONS, RELATIVE PRECISION:',G16.7)
201   FORMAT(1X,'GRACJG (BIEF) : ',
     *                     1I8,' ITERATIONS, PRECISION ABSOLUE :',G16.7)
202   FORMAT(1X,'GRACJG (BIEF) : ',
     *                     1I8,' ITERATIONS, ABSOLUTE PRECISION:',G16.7)
103   FORMAT(1X,'GRACJG (BIEF) : MAX D''ITERATIONS ATTEINT:',
     *                     1I8,' PRECISION RELATIVE:',G16.7)
104   FORMAT(1X,'GRACJG (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     *                     1I8,' RELATIVE PRECISION:',G16.7)
105   FORMAT(1X,'GRACJG (BIEF) : ',
     *         ' SOLUTION X=0 CAR NORME L2 DE B QUASI-NULLE:',G16.7)
106   FORMAT(1X,'GRACJG (BIEF) : ',
     *        ' SOLUTION X=0 BECAUSE L2-NORM OF B VERY SMALL:',G16.7)
203   FORMAT(1X,'GRACJG (BIEF) : MAX D''ITERATIONS ATTEINT:',
     *                     1I8,' PRECISION ABSOLUE :',G16.7)
204   FORMAT(1X,'GRACJG (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     *                     1I8,' ABSOLUTE PRECISION:',G16.7)
C
C-----------------------------------------------------------------------
C
      END
