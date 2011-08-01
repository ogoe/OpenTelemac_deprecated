C                       *****************
                        SUBROUTINE CGSQUA
C                       *****************
C
     *(X,A,B,MESH, G,G0,P,K,H,AHPK,CFG,INFOGR)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C                                        F  LEPEINTRE (LNH) 30 87 78 54
C***********************************************************************
C
C  FONCTION : RESOLUTION D'UN SYSTEME LINEAIRE A X = B
C             PAR LA METHODE DU GRADIENT CONJUGUE CARRE
C
C
C     ALGORITHME:
C
C        |
C        |   INITIALISATION
C        |   ---------------
C        |
C        |    0                          N
C        |   X  VECTEUR QUELCONQUE DANS R , APPROCHANT LA SOLUTION
C        |
C        |      0      0
C        |     G  = A X  - B
C        |
C        |     K0 = P0 = G0
C        |
C        |
C        |   ITERATIONS
C        |   ----------
C        |
C        |              M     0
C        |       M   ( K  ,  G  )
C        |     RO =  ------------
C        |               M    0
C        |           (A P ,  G  )
C        |
C        |      M     M    M     M
C        |     H   = K - RO * A P
C        |
C        |      M+1   M    M       M    M
C        |     G   = G - RO * A ( H  + K  )
C        |
C        |      M+1   M      M     M    M
C        |     X   = X   - RO * ( H  + K  )
C        |
C        |                 0     M+1
C        |              ( G  , G   )
C        |     BETA =   ------------
C        |                 0   M
C        |              ( G , G    )
C        |
C        |      M+1   M+1        M    M     M       M
C        |     P  =  G   + 2*BETA  * H + BETA**2 * P
C        |
C        |      M+1   M+1        M    M
C        |     K  =  G   +   BETA  * H
C        |
C        |
C
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
C        5            I  AUTRE
C-----------------------------------------------------------------------
C
C  SIGNIFICATION DE IELM :
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS          PROGRAMME ICI
C
C  11 : TRIANGLE P1            3                       OUI
C  12 : TRIANGLE P2            6                       OUI
C  13 : TRIANGLE P1-ISO P1     6                       OUI
C  14 : TRIANGLE P2            7
C  21 : QUADRILATERE Q1        4                       OUI
C  22 : QUADRILATERE Q2        8
C  24 : QUADRILATERE Q2        9
C  31 : TETRAEDRE P1           4                       OUI
C  32 : TETRAEDRE P2          10
C  41 : PRISMES MITHRIDATE     6                       OUI
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  VALEUR INITIALE, PUIS SOLUTION
C |      A         | -->|  MATRICE DU SYSTEME
C |      B         | -->|  SECOND MEMBRE DU SYSTEME.
C |      D         |<-->|  DIRECTION DE DESCENTE.
C |      AD        |<-->|  MATRICE A MULTIPLIEE PAR D.
C |      G         |<-->|  GRADIENT DE DESCENTE.
C |      AG        |<-->|  A X (GRADIENT DE DESCENTE).
C |      INFOGR    | -->|  SI OUI, ON IMPRIME UN COMPTE-RENDU
C |________________|____|______________________________________________-
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : MATRBL , OS , P_DOTS
C
C**********************************************************************
C
      USE BIEF, EX_CGSQUA => CGSQUA
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: X,G,G0,P,K,H,AHPK
      TYPE(BIEF_OBJ)  , INTENT(IN   ) :: B
      TYPE(BIEF_OBJ)  , INTENT(IN)    :: A
      TYPE(SLVCFG)    , INTENT(INOUT) :: CFG
      LOGICAL         , INTENT(IN)    :: INFOGR
      TYPE(BIEF_MESH) , INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION XL,RO,TESTL,RL,GMP1G0,BETA,GMG0,C
C
      INTEGER M
C
      LOGICAL RELAT
C
      INTRINSIC SQRT
C
C-----------------------------------------------------------------------
C
C CALCUL DE LA NORME DU SECOND MEMBRE
C
      XL = P_DOTS(B,B,MESH)
C
      IF( XL.LT.1.D0) THEN
            XL = 1.D0
            RELAT = .FALSE.
        ELSE
            RELAT = .TRUE.
      ENDIF
C
      M = 0
C
C INITIALISATION DE G  : A X0 - B
C
      CALL MATRBL( 'X=AY    ',G,A,X,C,  MESH)
C
      CALL OS( 'X=X-Y   ' , G , B , B , C )
C
C ON CONTROLE QUE LA PRECISION N'EST PAS DEJA ATTEINTE
C
      RL   = P_DOTS(G,G,MESH)
C
      IF(RL.LT.CFG%EPS**2*XL) THEN
        TESTL = SQRT(RL/XL)
        IF (INFOGR) THEN
          IF(RELAT) THEN
            IF (LNG.EQ.1) WRITE(LU,100) M,TESTL
            IF (LNG.EQ.2) WRITE(LU,101) M,TESTL
          ELSE
            IF (LNG.EQ.1) WRITE(LU,200) M,TESTL
            IF (LNG.EQ.2) WRITE(LU,201) M,TESTL
          ENDIF
        ENDIF
        GOTO 1000
      ENDIF
C
C INITIALISATION DE G0 , P , ET K
C
      CALL OS('X=Y     ' , G0 , G , G , C )
      CALL OS('X=Y     ' , P  , G , G , C )
      CALL OS('X=Y     ' , K  , G , G , C )
C
      M = 1
20    CONTINUE
C
C CALCUL DE AP (MIS DANS H QUI EST RECALCULE APRES)
C
      CALL MATRBL( 'X=AY    ',H,A,P,C,  MESH)
C
C CALCUL DE RO
C
      RO = P_DOTS(K,G0,MESH) / P_DOTS(H,G0,MESH)
C
C CALCUL DE H+K (MIS DANS H, AP ETANT DEJA DANS H)
C
      CALL OS( 'X=CX    ' , H , H , H , -RO   )
      CALL OS( 'X=X+CY  ' , H , K , K , 2.D0  )
C
C             M+1   M       M       M
C CALCUL DE  X   = X  -   RO * (H+K)     (H+K DANS H)
C
      CALL OS( 'X=X+CY  ' , X , H , H , -RO )
C
C CALCUL DE A(H+K)  (ICI H+K EST DANS H ET ON RANGE A(H+K) DANS AHPK)
C
      CALL MATRBL( 'X=AY    ',AHPK,A,H,C,  MESH)
C
C              M    0
C CALCUL DE ( G   , G  )
C
      GMG0 = P_DOTS(G,G0,MESH)
C
C CALCUL DE GM
C
      CALL OS( 'X=X+CY  ' , G , AHPK , AHPK , -RO )
C
      RL   = P_DOTS(G,G,MESH)
      IF (RL.GT.CFG%EPS**2*XL) THEN
         IF (M.GE.CFG%NITMAX) THEN
C          IF(INFOGR) THEN
                TESTL=SQRT(RL/XL)
                IF(RELAT) THEN
                  IF (LNG.EQ.1) WRITE(LU,102) M,TESTL
                  IF (LNG.EQ.2) WRITE(LU,103) M,TESTL
                ELSE
                  IF (LNG.EQ.1) WRITE(LU,202) M,TESTL
                  IF (LNG.EQ.2) WRITE(LU,203) M,TESTL
                ENDIF
C          ENDIF
           GOTO 1000
         ELSE
           M = M + 1
         ENDIF
      ELSE
         IF(INFOGR) THEN
           TESTL=SQRT(RL/XL)
           IF(RELAT) THEN
             IF (LNG.EQ.1) WRITE(LU,100) M,TESTL
             IF (LNG.EQ.2) WRITE(LU,101) M,TESTL
           ELSE
             IF (LNG.EQ.1) WRITE(LU,200) M,TESTL
             IF (LNG.EQ.2) WRITE(LU,201) M,TESTL
           ENDIF
         ENDIF
         GOTO 1000
      ENDIF
C
C              M+1  0
C CALCUL DE ( G   , G  )
C
      GMP1G0 = P_DOTS(G,G0,MESH)
C
C CALCUL DE BETA
C
      BETA = GMP1G0 / GMG0
C
C            M
C CALCUL DE H   (ON AVAIT MIS H+K DANS H)
C
      CALL OS('X=X-Y   ' , H , K , K , C )
C
C            M+1
C CALCUL DE P
C
      CALL OS('X=CX    ' , P , P , P , BETA**2 )
      CALL OS('X=X+Y   ' , P , G , G , C       )
      CALL OS('X=X+CY  ' , P , H , H , 2*BETA  )
C
C            M+1
C CALCUL DE K
C
      CALL OS('X=Y     ' , K , G , G , C    )
      CALL OS('X=X+CY  ' , K , H , H , BETA )
C
      GOTO 20
C
1000  RETURN
C
C-----------------------------------------------------------------------
C
C   FORMATS
C
100   FORMAT(1X,'CGSQUA (BIEF) : ',1I8,' ITERATIONS',
     *          ' PRECISION RELATIVE: ',G16.7)
101   FORMAT(1X,'CGSQUA (BIEF) : ',1I8,' ITERATIONS',
     *          ' RELATIVE PRECISION: ',G16.7)
200   FORMAT(1X,'CGSQUA (BIEF) : ',1I8,' ITERATIONS',
     *          ' PRECISION ABSOLUE: ',G16.7)
201   FORMAT(1X,'CGSQUA (BIEF) : ',1I8,' ITERATIONS',
     *          ' ABSOLUTE PRECISION: ',G16.7)
102   FORMAT(1X,'CGSQUA (BIEF) : MAX D''ITERATIONS ATTEINT ',1I8,
     *          ' PRECISION RELATIVE: ',G16.7)
103   FORMAT(1X,'CGSQUA (BIEF) : EXCEEDING MAXIMUM ITERATIONS ',1I8,
     *          ' RELATIVE PRECISION: ',G16.7)
202   FORMAT(1X,'CGSQUA (BIEF) : MAX D''ITERATIONS ATTEINT ',1I8,
     *          ' PRECISION ABSOLUE: ',G16.7)
203   FORMAT(1X,'CGSQUA (BIEF) : EXCEEDING MAXIMUM ITERATIONS ',1I8,
     *          ' ABSOLUTE PRECISION:',G16.7)
C
C-----------------------------------------------------------------------
C
      END 
 
 
