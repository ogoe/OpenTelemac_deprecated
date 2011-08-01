C                       ****************
                        SUBROUTINE GMRES
C                       ****************
C
     * (X,A,B,MESH,R0,V,AV,CFG,INFOGR,AUX)
C
C***********************************************************************
C BIEF VERSION 5.6           26/08/93    C  MOULIN    (LNH) 30 87 83 81
C                            24/04/97  J-M  HERVOUET  (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION : RESOLUTION D'UN SYSTEME LINEAIRE A X = B
C             PAR LA METHODE GMRES (GENERALISED MINIMUM RESIDUAL)
C
C
C  ATTENTION : LES PRECONDITIONNEMENTS CROUT ET GSEBE SONT ICI TRAITES
C              COMME UN PRECONDITIONNEMENT LU. POUR CROUT CELA CONDUIT
C              A NE PAS PRENDRE EN COMPTE LA DIAGONALE DE LA MATRICE DE
C              PRECONDITIONNEMENT.
C
C
C  ALGORITHME:
C
C        |
C        |   INITIALISATION
C        |   ---------------
C        |
C        |   R0 = B - A X0
C        |
C        |   V1 = R0 / IIR0II
C        |
C        |   K = DIMENSION CHOISIE POUR L'ESPACE DE KRYLOV
C        |
C        |
C        |
C        |
C        |   ITERATIONS
C        |   ----------
C        |
C        | 1) CONSTRUCTION D'UNE BASE ORTHONORMALE  (VI) DE DIM K+1
C        |
C        |   I = 1,...,K
C        |
C        |       VI+1 = A * VI
C        |
C        |       J = 1,...,I
C        |
C        |       HI+1,J = ( VI+1 , VJ )   MATRICE DE HESSENBERG (K+1,K)
C        |
C        |       VI+1 <--- VI+1  -  HI+1,J * VJ
C        |
C        |       VI+1 = VI+1 / IIVI+1II
C        |
C        |
C        | 2) MULTIPLICATION DE H PAR Q  =  RK * RK-1 * ... * R1
C        |
C        |  RI MATRICE DE ROTATION DE GIVENS :
C        |                                    -                     -
C        |                                   I ID(J-1)               I
C        |                                   I                       I
C        |                                   I        CJ  SJ         I
C        |                              RI = I       -SJ  CJ         I
C        |                                   I                       I
C        |                                   I              ID(K-J)  I
C        |                                    -                     -
C        |
C        |                         2    2
C        |       AVEC      CJ + SJ  =  1
C        |
C        |                 CJ ET SJ TELS QUE H DEVIENT TRIANGULAIRE
C        |
C        |                 ID(J-1) :    MATRICE IDENTITE J-1*J-1
C        |
C        |
C        | 3) RESOLUTION DU SYSTEME H Y = Q E
C        |                             -      -
C        |       E VECTEUR (NR0,0,0,0,0,0,.....)
C        |
C        |       NR0 NORME DU RESIDU
C        |
C        |
C        |
C        | 4) CALCUL DE X(M+1) = X(M)  +  V * Y
C        |
C        |  V : MATRICE (3*NPOIN,K) DONT LES COLONNES SONT LES VJ
C        |
C        | 5) TEST SUR LE RESIDU...
C        |
C
C-----------------------------------------------------------------------
C                         PRECONDITIONNEMENT
C-----------------------------------------------------------------------
C  VALEUR DE PRECON   I                  SIGNIFICATION
C-----------------------------------------------------------------------
C        0 OU 1       I  RIEN
C                     I
C        2            I  PRECONDITIONNEMENT DIAGONAL AVEC LA DIAGONALE
C                     I  DE LA MATRICE.
C                     I
C        3            I  PRECONDITIONNEMENT BLOC-DIAGONAL.
C                     I
C        5            I  PRECONDITIONNEMENT DIAGONAL AVEC LA VALEUR
C                     I  ABSOLUE DE LA DIAGONALE DE LA MATRICE.
C                     I
C        7            I  PRECONDITIONNEMENT DE CROUT PAR ELEMENT
C                     I
C       11            I  PRECONDITIONNEMENT DE GAUSS-SEIDEL PAR ELEMENT
C                     I
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |      X         |<-- |  VALEUR INITIALE, PUIS SOLUTION
C |      A         | -->|  MATRICE DU SYSTEME
C |      B         | -->|  SECOND MEMBRE DU SYSTEME.
C |      MESH      | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE.
C |      R0        | -->|  TABLEAU DE TRAVAIL DE DIMENSION NPOIN
C |      V         | -->|  TABLEAU DE TRAVAIL DE DIMENSION NPOIN
C |      AV        | -->|  TABLEAU DE TRAVAIL DE DIMENSION NPOIN
C |      INFOGR    | -->|  SI OUI, ON IMPRIME UN COMPTE-RENDU
C |      AUX       | -->|  MATRICE DE TRAVAIL UTILISEE AVEC LES
C |                |    |  PRECONDITIONNEMENTS DE TYPE CROUT.
C |________________|____|______________________________________________-
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : MATRBL , OS , DOTS
C
C**********************************************************************
C
      USE BIEF, EX_GMRES => GMRES
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(SLVCFG), INTENT(INOUT)    :: CFG
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: B
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: X,V,AV,R0
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ), INTENT(IN)     :: A
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: AUX
      LOGICAL, INTENT(IN)            :: INFOGR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
C  ICI ON PREND COMME MAXIMUM CFG%KRYLOV=20
C
      DOUBLE PRECISION H(21,20),C(20),S(20),E(21),BID
C
      DOUBLE PRECISION R,ZZ,NB,PREC,NR0
C
      INTEGER I,J,K,L,M
C
      LOGICAL RELAT,CROUT,GSEB,PRECO
C
      INTRINSIC SQRT,ABS
C
C-----------------------------------------------------------------------
C
      K = CFG%KRYLOV
C
      CROUT=.FALSE.
      IF(MOD(CFG%PRECON,7).EQ.0) CROUT=.TRUE.
      GSEB=.FALSE.
      IF(MOD(CFG%PRECON,11).EQ.0) GSEB=.TRUE.
      PRECO=.FALSE.
      IF(CROUT.OR.GSEB.OR.MOD(CFG%PRECON,13).EQ.0) PRECO=.TRUE.
C
      IF(PRECO) THEN
C                  -1
C       CALCUL DE L   B
        CALL GODOWN(B, AUX,B,'D',MESH,.FALSE.)
      ENDIF
C
C INITIALISATIONS
C
      NB = P_DOTS(B,B,MESH)
      NB = SQRT(NB)
C
      RELAT = .TRUE.
      IF(NB.LT.1.D0) THEN
        NB = 1.D0
        RELAT = .FALSE.
      ENDIF
C
      M = 0
C
C INITIALISATION DU RESIDU R : A X0 - B
C
      CALL MATRBL( 'X=AY    ',R0,A,X,BID,MESH)
C
      IF(PRECO) THEN
        CALL GODOWN(R0, AUX,R0,'D',MESH,.FALSE.)
        CALL PUOG(X,AUX,X, 'D',MESH,.FALSE.)
      ENDIF
C
      CALL OS( 'X=X-Y   ' , R0 , B , B , BID )
C
C ON CONTROLE QUE LA PRECISION N'EST PAS DEJA ATTEINTE
C
      NR0=P_DOTS(R0,R0,MESH)
      NR0=SQRT(NR0)
      PREC = NR0/NB
C
      IF (PREC.LE.CFG%EPS) GO TO 3000
C
C-----------------------------------------------------------------------
C                   BOUCLE DES ITERATIONS
C-----------------------------------------------------------------------
C
20    CONTINUE
C
      M = M+1
C
C CALCUL DU VECTEUR V1 = - R0 / NORME(R0)
C (SIGNE - CAR R = AX - B AU LIEU DE B - AX)
C
      CALL OS('X=CY    ', V%ADR(1)%P , R0 , R0 , -1.D0/NR0 )
C
C-----------------------------------------------------------------------
C         CALCUL DE LA BASE ORTHONORMALE ET DE LA MATRICE H
C-----------------------------------------------------------------------
C
C     K-1 PREMIERES COLONNES
C
      DO 10 J=1,K-1
C
        IF(PRECO) THEN
          CALL GOUP(B , AUX , V%ADR(J)%P ,'D',MESH,.TRUE.)
          CALL MATRBL( 'X=AY    ',AV%ADR(J)%P,A,B,BID,MESH)
          CALL GODOWN(AV%ADR(J)%P,AUX,AV%ADR(J)%P,'D',MESH,.FALSE.)
        ELSE
          CALL MATRBL( 'X=AY    ',AV%ADR(J)%P,A,V%ADR(J)%P,BID,MESH)
        ENDIF
C
        CALL OS('X=Y     ', V%ADR(J+1)%P ,AV%ADR(J)%P , X , BID )
C
        DO 30 I = 1,J
C
          H(I,J) = P_DOTS( V%ADR(J+1)%P , V%ADR(I)%P , MESH )    
C
          CALL OS('X=X+CY  ',V%ADR(J+1)%P,V%ADR(I)%P,V%ADR(I)%P,-H(I,J))
C
30      CONTINUE
C
         H(J+1,J)=P_DOTS(V%ADR(J+1)%P,V%ADR(J+1)%P,MESH)
         H(J+1,J) = SQRT( H(J+1,J) )    
C
         CALL OS('X=CX    ',V%ADR(J+1)%P, B, B, 1.D0/H(J+1,J))
C
10    CONTINUE
C
C K-IEME COLONNE (ON NE CONSTRUIT PAS ENTIEREMENT LE VECTEUR V(K+1) )
C
        IF(PRECO) THEN
          CALL GOUP(B , AUX , V%ADR(K)%P , 'D' ,MESH,.TRUE.)
          CALL MATRBL( 'X=AY    ',AV%ADR(K)%P,A,B,BID,MESH)
          CALL GODOWN(AV%ADR(K)%P,AUX,AV%ADR(K)%P,'D',MESH,.FALSE.)
        ELSE
          CALL MATRBL( 'X=AY    ',AV%ADR(K)%P,A,V%ADR(K)%P,BID,MESH)
        ENDIF
C

        H(K+1,K) = P_DOTS( AV%ADR(K)%P , AV%ADR(K)%P , MESH )
C
        DO 31 I = 1,K
C
          H(I,K) = P_DOTS( AV%ADR(K)%P , V%ADR(I)%P , MESH )
          H(K+1,K) = H(K+1,K) - H(I,K)**2
C
31      CONTINUE
C       EN PRINCIPE H(K+1,K) EST POSITIF
C       A LA PRECISION DE LA MACHINE PRES
        H(K+1,K) = SQRT( ABS(H(K+1,K)) )
C
C-----------------------------------------------------------------------
C CONSTRUCTION DES ROTATIONS DE GIVENS ET APPLICATION A H ET E
C-----------------------------------------------------------------------
C
C     AUTRES COMPOSANTES DE E NULLES
      E(1) = NR0
C
C  ROTATIONS DE 1 A K
C
      DO 40 I = 1 , K
C
C     MODIFICATION DE LA COLONNE I DE H PAR LES ROTATIONS PRECEDENTES
      IF(I.GE.2) THEN
        DO 41 J = 1,I-1
          ZZ       =  C(J) * H(J,I) + S(J) * H(J+1,I)
          H(J+1,I) = -S(J) * H(J,I) + C(J) * H(J+1,I)
          H(J,I) = ZZ
41      CONTINUE
      ENDIF
C     MODIFICATION DE LA COLONNE I DE H PAR LA ROTATION I
      R = SQRT( H(I,I)**2 + H(I+1,I)**2 )
      IF(ABS(R).LT.1.D-6) THEN
        IF(INFOGR) THEN
          IF (LNG.EQ.1) WRITE(LU,91) R
          IF (LNG.EQ.2) WRITE(LU,92) R
        ENDIF
        GO TO 3000
      ENDIF
      C(I) =  H(I,I)   / R
      S(I) =  H(I+1,I) / R
      H(I,I) = R
C     H(I+1,I) = 0.D0    (ON NE S'EN RESERVIRA PAS)
C     MODIFICATION DU VECTEUR E
      E(I+1) = -S(I) * E(I)
      E(I  ) =  C(I) * E(I)
C
40    CONTINUE
C
C-----------------------------------------------------------------------
C RESOLUTION DU SYSTEME H*Y = E     (H TRIANGULAIRE SUP DE DIMENSION K)
C                                    Y CONFONDU AVEC E
C
C LE FAIT QUE H(I,I) EST NON NUL A ETE VERIFIE PLUS HAUT SUR R
C-----------------------------------------------------------------------
C
      E(K) = E(K) / H(K,K)
      DO 120 J = K-1,1,-1
      DO 130 L = J+1,K
        E(J) = E(J) - H(J,L) * E(L)
130   CONTINUE
      E(J) = E(J) / H(J,J)
120   CONTINUE
C
C-----------------------------------------------------------------------
C ON FORME LA SOLUTION POUR L'ETAPE M : X(M+1) = X(M) + VK * Y(K)
C-----------------------------------------------------------------------
C
      DO 150 L = 1,K
C
C  CALCUL DE LA NOUVELLE SOLUTION X
C
        CALL OS('X=X+CY  ', X , V%ADR(L)%P , X , E(L) )
C
C  CALCUL DU RESIDU : RM+1 = RM + A * VK * ZK
C
        CALL OS('X=X+CY  ', R0 , AV%ADR(L)%P , R0 , E(L) )
C
150   CONTINUE
C
C ON CONTROLE QUE LA PRECISION N'EST PAS ATTEINTE
C ON PEUT MONTRER QUE LA NORME DU RESIDU EST ABS(E(K+1))
C
C     NR0  =norme de R0
      NR0  = ABS(E(K+1))
      PREC = NR0/NB
C
      IF (PREC.GE.CFG%EPS.AND.M.LT.CFG%NITMAX)  GOTO 20
C
3000  CONTINUE
C
C-----------------------------------------------------------------------
C
      IF(PRECO) THEN
        CALL GOUP( X,AUX,X,'D',MESH,.FALSE.)
      ENDIF
C
C-----------------------------------------------------------------------
C
C     IF(INFOGR) THEN
        IF(M.LT.CFG%NITMAX) THEN
          IF(INFOGR) THEN
          IF(RELAT) THEN
            IF (LNG.EQ.1) WRITE(LU,101) M,PREC
            IF (LNG.EQ.2) WRITE(LU,102) M,PREC
          ELSE
            IF (LNG.EQ.1) WRITE(LU,201) M,PREC
            IF (LNG.EQ.2) WRITE(LU,202) M,PREC
          ENDIF
          ENDIF
        ELSE
          IF(RELAT) THEN
            IF (LNG.EQ.1) WRITE(LU,103) M,PREC
            IF (LNG.EQ.2) WRITE(LU,104) M,PREC
          ELSE
            IF (LNG.EQ.1) WRITE(LU,203) M,PREC
            IF (LNG.EQ.2) WRITE(LU,204) M,PREC
          ENDIF
        ENDIF
C     ENDIF
C
      RETURN
C
C-----------------------------------------------------------------------
C
 91   FORMAT(1X,'GMRES (BIEF) : ECHEC DE L''ALGORITHME  R= ',G16.7,/,
     *       1X,'               SI LA MATRICE EST DIAGONALE, PRENDRE',/,
     *       1X,'               COMME SOLVEUR LE GRADIENT CONJUGUE.')
 92   FORMAT(1X,'GMRES (BIEF) : ALGORITHM FAILED  R=',G16.7,/,
     *       1X,'               IF THE MATRIX IS DIAGONAL, CHOOSE',/,
     *       1X,'               THE CONJUGATE GRADIENT SOLVER.')
101   FORMAT(1X,'GMRES (BIEF) : ',
     *                     1I8,' ITERATIONS, PRECISION RELATIVE:',G16.7)
102   FORMAT(1X,'GMRES (BIEF) : ',
     *                     1I8,' ITERATIONS, RELATIVE PRECISION:',G16.7)
201   FORMAT(1X,'GMRES (BIEF) : ',
     *                     1I8,' ITERATIONS, PRECISION ABSOLUE :',G16.7)
202   FORMAT(1X,'GMRES (BIEF) : ',
     *                     1I8,' ITERATIONS, ABSOLUTE PRECISION:',G16.7)
103   FORMAT(1X,'GMRES (BIEF) : MAX D'' ITERATIONS ATTEINT:',
     *                     1I8,' PRECISION RELATIVE:',G16.7)
104   FORMAT(1X,'GMRES (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     *                     1I8,' RELATIVE PRECISION:',G16.7)
203   FORMAT(1X,'GMRES (BIEF) : MAX D'' ITERATIONS ATTEINT:',
     *                     1I8,' PRECISION ABSOLUE :',G16.7)
204   FORMAT(1X,'GMRES (BIEF) : EXCEEDING MAXIMUM ITERATIONS:',
     *                     1I8,' ABSOLUTE PRECISION:',G16.7)
C
C-----------------------------------------------------------------------
C
      END
