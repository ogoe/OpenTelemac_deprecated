C                       *****************
                        SUBROUTINE FILTER
C                       *****************
C
     *(VEC,BLDMAT,T1,T2,
     * A,FORMUL,
     * XMUL,F,G,H,U,V,W,
     * MESH,MSK,MASKEL,N)
C
C***********************************************************************
C BIEF VERSION 5.1           24/04/97    J-M HERVOUET (LNH) 30 87 80 18
C
C***********************************************************************
C
C  FONCTION  : FILTRE D'UN VECTEUR A L'AIDE D'UNE MATRICE
C
C              AVEC UNE MATRICE DE MASSE ON OBTIENT PAR EXEMPLE
C              UN LISSAGE.
C
C              NOTE : SI BLDMAT=.FALSE. ON SUPPOSE A DONNEE
C                     ET ELLE N'EST PAS RECONSTRUITE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |    VEC         |<-->| VECTEUR A FILTRER
C |    BLDMAT      | -->| LOGIQUE : ON CONSTRUIT LA MATRICE OU PAS.
C |    T1          | -->| TABLEAU DE TRAVAIL.
C |    T2          |<-->| TABLEAU DE TRAVAIL. MATRICE A MASS-LUMPEE
C |                |    | EN SORTIE (VOIR AUSSI XMUL)
C |    A           |<-->| MATRICE (DONNEE OU CONSTRUITE SUIVANT BLDMAT)
C |    FORMUL      | -->| FORMULE DECRIVANT LA MATRICE
C |                |    | (MEMES CONVENTIONS QUE DANS MATRIX)
C |    XMUL        | -->| FACTEUR MULTIPLICATIF NON NUL
C |                |    | N'A AUCUNE INFLUENCE SAUF SUR T2.
C |    F,G,H,U,V,W | -->| FONCTIONS INTERVENANT DANS LA MATRICE
C |    MESH,       | -->| BLOCS DU MAILLAGE.
C |    MSK,MASKEL  | -->| LOGIQUE ET TABLEAU POUR LE MASQUAGE
C |    N           | -->| NOMBRE DE FOIS OU ON FAIT L'OPERATION.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
C  SOUS-PROGRAMMES APPELES :  LUMP
C
C***********************************************************************
C
      USE BIEF, EX_FILTER => FILTER
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: N
      DOUBLE PRECISION, INTENT(IN)  :: XMUL
      LOGICAL, INTENT(IN)           :: BLDMAT,MSK
      CHARACTER(LEN=16), INTENT(IN) :: FORMUL
      TYPE(BIEF_MESH), INTENT(INOUT):: MESH
      TYPE(BIEF_OBJ), INTENT(INOUT) :: VEC,A,T1,T2
      TYPE(BIEF_OBJ), INTENT(IN)    :: F,G,H,U,V,W,MASKEL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
C
      DOUBLE PRECISION C
C
C-----------------------------------------------------------------------
C
      DO 10 I=1,N
C
C  CALCUL (EVENTUEL) DE LA MATRICE SUIVANT LA FORMULE DONNEE
C
      IF(BLDMAT.AND.I.EQ.1) THEN
C
          CALL MATRIX(A,'M=N     ',FORMUL,VEC%ELM,VEC%ELM,
     *                XMUL,F,G,H,U,V,W,MESH,MSK,MASKEL)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CALCUL DU PRODUIT A * VEC (AVEC ASSEMBLAGE)
C
      CALL MATVEC( 'X=AY    ',T1,A,VEC,C,MESH)
      IF(NCSIZE.GT.1) CALL PARCOM(T1,2,MESH)
C
C-----------------------------------------------------------------------
C
C  CONDENSATION DE A SUR SA DIAGONALE
C
      IF(I.EQ.1) THEN
        CALL LUMP(T2,A,MESH,XMUL)
        IF(NCSIZE.GT.1) CALL PARCOM(T2,2,MESH)
        CALL OS('X=1/Y   ',T2,T2,T2,C,IOPT=2,INFINI=0.D0,ZERO=1.D-20)
      ENDIF
C
C-----------------------------------------------------------------------
C
C  CALCUL DE F = A * F / (AGGLOMEREE DE A)
C
C  ON VERIFIE LES DIVISIONS PAR ZERO A CAUSE DES POINTS EXTERIEURS DU
C  FORMAT LEONARD, QUI PEUVENT AVOIR DES VALEURS NULLES.
C
      CALL OS('X=YZ    ',X=VEC,Y=T1,Z=T2)
C
C-----------------------------------------------------------------------
C
10    CONTINUE
C
C-----------------------------------------------------------------------
C
      RETURN
      END 
 
