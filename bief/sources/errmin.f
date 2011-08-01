C                       *****************
                        SUBROUTINE ERRMIN
C                       *****************
C
     *(X, A,B , MESH, D,AD,G,R, CFG,INFOGR,AUX)
C
C***********************************************************************
C BIEF VERSION 5.6           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : RESOLUTION D'UN SYSTEME LINEAIRE A X = B
C             PAR DES METHODES DE TYPE GRADIENT CONJUGUE :
C
C-----------------------------------------------------------------------
C                        PRECONDITIONNEMENT
C-----------------------------------------------------------------------
C  VALEUR DE PRECON   I                  SIGNIFICATION
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
C |                |    |  PAS DE PRECONDITIONNEMENT DANS SOLGRA)
C |      INFOGR    | -->|  SI OUI, IMPRESSION D'UN COMPTE-RENDU.
C |      AUX       | -->|  MATRICE POUR LE PRECONDITIONNEMENT.
C |________________|____|______________________________________________-
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : MATRBL , OS
C
C**********************************************************************
C
      USE BIEF, EX_ERRMIN => ERRMIN
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
      DOUBLE PRECISION XL,RMRM,TESTL,DD
      DOUBLE PRECISION BETA,RO,GMGM,GM1GM1
      DOUBLE PRECISION STO2,TGMTGM,C
C
      LOGICAL RELAT,PREC,CROUT,GSEB
C
C-----------------------------------------------------------------------
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
C   INITIALISATIONS
C     STO2 POUR EVITER UN "CAUTION DU COMPILATEUR CRAY"
      STO2  =0.D0
      CROUT =.FALSE.
      IF(7*(CFG%PRECON/7).EQ.CFG%PRECON) CROUT=.TRUE.
      GSEB=.FALSE.
      IF(11*(CFG%PRECON/11).EQ.CFG%PRECON) GSEB=.TRUE.
      PREC=.FALSE.
      IF(CROUT.OR.GSEB.OR.13*(CFG%PRECON/13).EQ.CFG%PRECON) PREC=.TRUE.
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
      GMGM   = RMRM
C
      IF (RMRM.LT.CFG%EPS**2*XL) GO TO 900
C
C-----------------------------------------------------------------------
C PRECONDITIONNEMENT :
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
C
C       CALCUL DE C G0 = R
        CALL DOWNUP(G, AUX , R , 'D' , MESH)
C
C       ICI C EST CONSIDEREE COMME SYMETRIQUE, SINON IL FAUDRAIT
C       RESOUDRE TC GPRIM = G
C
C                       T -1
C       DANS B  ONT MET  C   G
        CALL DOWNUP(B , AUX , G , 'T' , MESH)
        GMGM   = P_DOTS(G,G,MESH)
        STO2   = GMGM
C
      ENDIF
C
C-----------------------------------------------------------------------
C CALCUL DE LA DIRECTION DE DESCENTE INITIALE :
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
        CALL MATRBL( 'X=TAY   ',D,A,B,  C,MESH)
      ELSE
        CALL MATRBL( 'X=TAY   ',D,A,G,  C,MESH)
      ENDIF
C
      TGMTGM = P_DOTS(D,D,MESH)
C
C-----------------------------------------------------------------------
C CALCUL DU PRODUIT A D INITIAL :
C-----------------------------------------------------------------------
C
      CALL MATRBL( 'X=AY    ',AD,A,D,  C,MESH)
C
C-----------------------------------------------------------------------
C CALCUL DE RO INITIAL :
C-----------------------------------------------------------------------
C
      RO = GMGM/TGMTGM
C
C-----------------------------------------------------------------------
C
C CALCUL DE X1 = X0 - RO  * D
C
      CALL OS( 'X=X+CY  ' , X , D , D , -RO )
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
      CALL OS( 'X=X+CY  ' , R , AD , AD , -RO )
C
C  CERTAINES VALEURS SERONT CHANGEES EN CAS DE PRECONDITIONNEMENT
C
      GM1GM1 = GMGM
      RMRM   = P_DOTS(R,R,MESH)
      GMGM   = RMRM
C
C TEST DE FIN :
C
      IF (RMRM.LE.XL*CFG%EPS**2) GO TO 900
C
C-----------------------------------------------------------------------
C PRECONDITIONNEMENT : RESOLUTION DE C G = R
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
C
C       RESOLUTION DE C G = R
        CALL DOWNUP(G, AUX , R , 'D' , MESH)
C
        CALL DOWNUP(B , AUX , G , 'T' , MESH)
        GM1GM1 = STO2
        GMGM = P_DOTS(G,G,MESH)
        STO2 = GMGM
C
      ENDIF
C
C-----------------------------------------------------------------------
C CALCUL DE D PAR RECURRENCE :
C-----------------------------------------------------------------------
C
      IF(PREC) THEN
C                                          T  T -1          T -1
C                               AD EST ICI  A  C  G   B EST  C   G
        CALL MATRBL( 'X=TAY   ',AD,A,B,  C,MESH)
      ELSE
C                               AD EST ICI TAG
        CALL MATRBL( 'X=TAY   ',AD,A,G,  C,MESH)
      ENDIF
C
      BETA = GMGM/GM1GM1
C
      CALL OS( 'X=CX    ' , D , D , D , BETA )
C
C                               AD EST ICI TAG
      CALL OS( 'X=X+Y   ' , D , AD , AD , C   )
C
C-----------------------------------------------------------------------
C CALCUL DE A D :
C-----------------------------------------------------------------------
C
      CALL MATRBL( 'X=AY    ',AD,A,D,  C,MESH)
C
C-----------------------------------------------------------------------
C CALCUL DE RO
C-----------------------------------------------------------------------
C
      DD = P_DOTS(D,D,MESH)
      RO = GMGM/DD
C
C CALCUL DE X(M) = X(M-1) - RO * D
C
      CALL OS( 'X=X+CY  ' , X , D , D , -RO )
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
1000  RETURN
C
C-----------------------------------------------------------------------
C
C   FORMATS
C
101   FORMAT(1X,'ERRMIN (BIEF) : ',
     *                     1I8,' ITERATIONS, PRECISION RELATIVE:',G16.7)
102   FORMAT(1X,'ERRMIN (BIEF) : ',
     *                     1I8,' ITERATIONS, RELATIVE PRECISION:',G16.7)
201   FORMAT(1X,'ERRMIN (BIEF) : ',
     *                     1I8,' ITERATIONS, PRECISION ABSOLUE :',G16.7)
202   FORMAT(1X,'ERRMIN (BIEF) : ',
     *                     1I8,' ITERATIONS, ABSOLUTE PRECISION:',G16.7)
103   FORMAT(1X,'ERRMIN (BIEF) : MAX D'' ITERATIONS ATTEINT:',
     *                     1I8,' PRECISION RELATIVE:',G16.7)
104   FORMAT(1X,'ERRMIN (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     *                     1I8,' RELATIVE PRECISION:',G16.7)
203   FORMAT(1X,'ERRMIN (BIEF) : MAX D'' ITERATIONS ATTEINT:',
     *                     1I8,' PRECISION ABSOLUE :',G16.7)
204   FORMAT(1X,'ERRMIN (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     *                     1I8,' ABSOLUTE PRECISON:',G16.7)
C
C-----------------------------------------------------------------------
C
      END 
 
 
