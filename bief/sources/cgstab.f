C                       *****************
                        SUBROUTINE CGSTAB
C                       *****************
C
     *(X, A,B , MESH, P,Q,R,S,T,V, CFG,INFOGR,AUX)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    R RATKE      (HANNOVER)
C                                        A MALCHEREK  (HANNOVER)
C                                        J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C  FONCTION : RESOLUTION D'UN SYSTEME LINEAIRE A X = B
C             PAR UNE METHODE DE GRADIENT CONJUGUE CARRE
C-----------------------------------------------------------------------
C                        PRECONDITIONNEMENT
C                        (VOIR AUSSI SOLV01)
C-----------------------------------------------------------------------
C  VALEUR DE PRECON   I                  SIGNIFICATION
C-----------------------------------------------------------------------
C        0            I  RIEN
C        2            I  PRECONDITIONNEMENT DIAGONAL AVEC LA DIAGONALE
C                     I  DE LA MATRICE.
C        3            I  PRECONDITIONNEMENT DIAGONAL AVEC LA MATRICE
C                     I  CONDENSEE.
C        5            I  DIAGONAL MAIS EN POUVANT AVOIR DES VALEURS
C                     I  NULLES OU NEGATIVES SUR LA DIAGONALE
C        7            I  CROUT
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  VALEUR INITIALE, PUIS SOLUTION
C |      A         | -->|  MATRICE DU SYSTEME
C |      B         | -->|  SECOND MEMBRE DU SYSTEME.
C |   P,Q,R,S,T,V  |<-->|  TABLEAUX DE TRAVAIL
C |      INFOGR    | -->|  SI OUI, ON IMPRIME UN COMPTE-RENDU
C |      AUX       | -->|  MATRICE DE PRECONDITIONNEMENT.
C |________________|____|______________________________________________-
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : MATRBL , OS , P_DOTS
C
C**********************************************************************
C
      USE BIEF, EX_CGSTAB => CGSTAB
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: X,P,Q,R,S,T,V
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: AUX,A,B
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH
      TYPE(SLVCFG)    , INTENT(INOUT) :: CFG
      LOGICAL         , INTENT(IN)    :: INFOGR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION ALFA,ALFA1,BETA,BETA1,OMEG,OMEG1,OMEG2
      DOUBLE PRECISION XL,TESTL,RMRM,C
C
      INTEGER M
C
      LOGICAL RELAT,CROUT
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
C   INITIALISATIONS
      CROUT=.FALSE.
      IF(7*(CFG%PRECON/7).EQ.CFG%PRECON) CROUT=.TRUE.
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
      IF (XL.LT.1.D0) THEN
         XL = 1.D0
         RELAT = .FALSE.
      ELSE
         RELAT = .TRUE.
      ENDIF
C
C SI LE SECOND MEMBRE EST NUL, X=0 ET C'EST FINI
C
      IF(XL.LT.CFG%ZERO**2) THEN
        RMRM = 0.D0
        CALL OS( 'X=C     ' , X , X , X , 0.D0 )
        GOTO 900
      ENDIF
C
C CALCUL DU RESIDU INITIAL ET SORTIE EVENTUELLE :
C
      CALL MATRBL( 'X=AY    ',V,A,X,C,  MESH)
C
      CALL OS( 'X=Y-Z   ' , R , B , V , C )
      RMRM   = P_DOTS(R,R,MESH)
C
      IF (RMRM.LT.CFG%EPS**2*XL) GO TO 900
C
C-----------------------------------------------------------------------
C PRECONDITIONNEMENT :
C-----------------------------------------------------------------------
C
      IF(CROUT) THEN
C       CALCUL DE C R  = B
        CALL DOWNUP(R, AUX , B , 'D' , MESH)
        IF(NCSIZE.GT.1) CALL PARMOY(R,MESH)
      ELSE
        CALL OS( 'X=Y     ' , R , B , B , C )
      ENDIF
C
C-----------------------------------------------------------------------
C SUITE DES INITIALISATIONS
C-----------------------------------------------------------------------
C
      IF(CROUT) THEN
        CALL DOWNUP(V, AUX , V , 'D' , MESH)
        IF(NCSIZE.GT.1) CALL PARMOY(V,MESH)
      ENDIF
C
      CALL OS( 'X=X-Y   ' , R , V , V , C    )
      CALL OS( 'X=Y     ' , P , R , R , C    )
      CALL OS( 'X=C     ' , V , V , V , 0.D0 )
      CALL OS( 'X=C     ' , Q , Q , Q , 0.D0 )
C
      ALFA  = 1.D0
      BETA  = 1.D0
      OMEG1 = 1.D0
C
C-----------------------------------------------------------------------
C  BOUCLE DES ITERATIONS :
C-----------------------------------------------------------------------
C
2     M  = M  + 1
C
      BETA1 = P_DOTS(R,P,MESH)
      OMEG2 = OMEG1*BETA1/BETA
      OMEG  = OMEG2/ALFA
      BETA  = BETA1
C
      CALL OS( 'X=Y+CZ  ' , Q , R    , Q ,  OMEG )
      CALL OS( 'X=X+CY  ' , Q , V    , V , -OMEG2)
C
      CALL MATRBL( 'X=AY    ',V,A,Q,C,  MESH)
C
      IF(CROUT) THEN
        CALL DOWNUP(V, AUX , V , 'D' , MESH)
        IF(NCSIZE.GT.1) CALL PARMOY(V,MESH)
      ENDIF
C
      OMEG1 = P_DOTS(P,V,MESH)
      OMEG1 = BETA1/OMEG1
C
      CALL OS( 'X=Y+CZ  ' , S , R    , V , -OMEG1)
C
      CALL MATRBL( 'X=AY    ',T,A,S,C,  MESH)
C
      IF(CROUT) THEN
        CALL DOWNUP(T, AUX , T , 'D' , MESH)
        IF(NCSIZE.GT.1) CALL PARMOY(T,MESH)
      ENDIF
C
      ALFA  = P_DOTS(T,S,MESH)
      ALFA1 = P_DOTS(T,T,MESH)
      ALFA  = ALFA/ALFA1
C
      CALL OS( 'X=X+CY  ' , X , Q , Q ,  OMEG1)
      CALL OS( 'X=X+CY  ' , X , S , S ,  ALFA )
C
      CALL OS( 'X=Y+CZ  ' , R , S , T , -ALFA )
C
      RMRM   = P_DOTS(R,R,MESH)
C
C TEST DE FIN :
C
      IF (RMRM.LE.XL*CFG%EPS**2) GO TO 900
C
      IF(M.LT.CFG%NITMAX) GO TO 2
C
C-----------------------------------------------------------------------
C
C     IF(INFOGR) THEN
        TESTL = SQRT( RMRM / XL )
        IF(RELAT) THEN
          IF (LNG.EQ.1) WRITE(LU,102) M,TESTL
          IF (LNG.EQ.2) WRITE(LU,104) M,TESTL
        ELSE
          IF (LNG.EQ.1) WRITE(LU,202) M,TESTL
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
        IF(RELAT) THEN
          IF (LNG.EQ.1) WRITE(LU,101) M,TESTL
          IF (LNG.EQ.2) WRITE(LU,103) M,TESTL
        ELSE
          IF (LNG.EQ.1) WRITE(LU,201) M,TESTL
          IF (LNG.EQ.2) WRITE(LU,203) M,TESTL
        ENDIF
      ENDIF
C
1000  RETURN
C
C-----------------------------------------------------------------------
C
C   FORMATS
C
101   FORMAT(1X,'CGSTAB (BIEF) : ',1I8,' ITERATIONS',
     *          ' PRECISION RELATIVE: ',G16.7)
103   FORMAT(1X,'CGSTAB (BIEF) : ',1I8,' ITERATIONS',
     *          ' RELATIVE PRECISION: ',G16.7)
201   FORMAT(1X,'CGSTAB (BIEF) : ',1I8,' ITERATIONS',
     *          ' PRECISION ABSOLUE: ',G16.7)
203   FORMAT(1X,'CGSTAB (BIEF) : ',1I8,' ITERATIONS',
     *          ' ABSOLUTE PRECISION: ',G16.7)
102   FORMAT(1X,'CGSTAB (BIEF) : MAX D''ITERATIONS ATTEINT ',1I8,
     *          ' PRECISION RELATIVE: ',G16.7)
104   FORMAT(1X,'CGSTAB (BIEF) : EXCEEDING MAXIMUM ITERATIONS ',1I8,
     *          ' RELATIVE PRECISION: ',G16.7)
202   FORMAT(1X,'CGSTAB (BIEF) : MAX D''ITERATIONS ATTEINT ',1I8,
     *          ' PRECISION ABSOLUE: ',G16.7)
204   FORMAT(1X,'CGSTAB (BIEF) : EXCEEDING MAXIMUM ITERATIONS',1I8,
     *          ' ABSOLUTE PRECISION:',G16.7)
C
C-----------------------------------------------------------------------
C
      END 
 
 
