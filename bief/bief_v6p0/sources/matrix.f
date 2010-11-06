C                       *****************
                        SUBROUTINE MATRIX
C                       *****************
C
     *(M,OP,FORMUL,IELM1,IELM2,XMUL,F,G,H,U,V,W,MESH,MSK,MASKEL)
C
C***********************************************************************
C BIEF VERSION 6.0      25/06/2008     JM HERVOUET (LNHE) 01 30 87 80 18
C
C
C 25/06/2008 JMH : PAS D'APPEL DE MATRIY SI NOMBRE D'ELEMENTS NUL
C 14/10/2009 JMH : ARGUMENTS ADDED TO ASSEX3
C
C***********************************************************************
C
C  FONCTION : CALCULS DE MATRICES
C
C             LA MATRICE EST IDENTIFIEE PAR LA FORMULE CONTENUE DANS
C             LA CHAINE DE CARACTERES FORMUL
C
C             OP EST UNE CHAINE DE 8 CARACTERES LA FACON DONT M EST
C             MODIFIEE. LA SYNTAXE EST LA MEME QUE OM, PAR EXEMPLE
C
C             OP='M=N     '
C             OP='M=TN    '
C             OP='M=M+N   '
C             OP='M=M+TN  '
C
C             TOUTES LES OPERATIONS DE OM QUI CONTIENNENT N A DROITE DU
C             SIGNE = SONT VALIDES.
C
C-----------------------------------------------------------------------
C
C  SIGNIFICATION DE IELM ET IELM2
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS
C
C  10 : TRIANGLE P0            1
C  11 : TRIANGLE P1            3
C  12 : TRIANGLE QUASI-BULLE   4
C
C  20 : QUADRILATERE Q0        1
C  21 : QUADRILATERE Q1        4
C
C  40 : PRISMES TELEMAC-3D P0  1
C  41 : PRISMES TELEMAC-3D P1  6
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  M             |<-->|  MATRICE A REMPLIR OU MODIFIER
C |  OP            | -->|  VOIR PLUS HAUT.
C |  FORMUL        | -->|  FORMULE DECRIVANT LA MATRICE
C |  IELM1         | -->|  TYPE D'ELEMENT POUR LES LIGNES
C |  IELM2         | -->|  TYPE D'ELEMENT POUR LES COLONNES
C |  XMUL          | -->|  COEFFICIENT MULTIPLICATEUR DU RESULTAT
C |  F,G,H         | -->|  FONCTIONS INTERVENANT DANS LA FORMULE
C |  U,V,W         | -->|  COMPOSANTES D'UN VECTEUR U DANS LA FORMULE
C |  MESH          | -->|  BLOC DES TABLEAUX D'ENTIERS DU MAILLAGE
C |  MSK           | -->|  SI OUI, PRESENCE D'ELEMENTS MASQUES.
C |  MASKEL        | -->|  TABLEAU DE MASQUAGE DES ELEMENTS
C |                |    |  =1. : NORMAL   =0. : ELEMENT MASQUE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : MATRIY
C
C**********************************************************************
C
      USE BIEF, EX_MATRIX => MATRIX
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)            :: IELM1,IELM2
C
      DOUBLE PRECISION, INTENT(IN)   :: XMUL
C
      LOGICAL, INTENT(IN)            :: MSK
C
      CHARACTER(LEN=16), INTENT(IN)  :: FORMUL
      CHARACTER(LEN=8), INTENT(IN)   :: OP
C
      TYPE(BIEF_OBJ), INTENT(IN)     :: F,G,H,U,V,W,MASKEL
      TYPE(BIEF_OBJ), INTENT(INOUT)  :: M
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER NELMAX,NELEM,NPT,SS
      INTEGER, DIMENSION(:), POINTER :: IKLE
      DOUBLE PRECISION, DIMENSION(:), POINTER :: SURFAC,XX,YY,ZZ
      DOUBLE PRECISION C
      LOGICAL LEGO
C
C-----------------------------------------------------------------------
C
C     STOCKAGE 1 POUR LA MATRICE DE TRAVAIL.
C     ASSEXT PEUT ENSUITE LE MODIFIER
C
      MESH%M%STO = 1
C
C-----------------------------------------------------------------------
C  APPEL DU SOUS-PROGRAMME D'AIGUILLAGE ET D'ASSEMBLAGE : MATRIY
C-----------------------------------------------------------------------
C
C     LEGO CAN BE MODIFIED BY MATRIY
      LEGO = .TRUE.
C
      IF(DIMENS(IELM1).EQ.MESH%DIM) THEN
C       MATRICE NORMALE
        NELEM  = MESH%NELEM
        NELMAX = MESH%NELMAX
        IKLE   =>MESH%IKLE%I
        SURFAC =>MESH%SURFAC%R
        XX=>MESH%XEL%R
        YY=>MESH%YEL%R
        ZZ=>MESH%ZEL%R
      ELSE
C       MATRICE DE BORD
        NELEM  = MESH%NELEB
        NELMAX = MESH%NELEBX
        IKLE   =>MESH%IKLBOR%I
        SURFAC =>MESH%LGSEG%R
        XX=>MESH%X%R
        YY=>MESH%Y%R
        ZZ=>MESH%Z%R
      ENDIF
C
C     MATRIY FILLS THE DIAGONAL AND EXTRA DIAGONAL TERMS
C
C     CHOICE OF PRE-ASSEMBLY STORAGE
      IF(M%STO.EQ.1.OR.M%STO.EQ.3) THEN
        SS = 1
      ELSEIF(M%STO.EQ.2) THEN
        SS = 2
      ELSE
        WRITE(LU,*) 'UNKNOWN STORAGE IN MATRIX'
        CALL PLANTE(1)
        STOP
      ENDIF
C
      IF(NELEM.GT.0) THEN
        CALL MATRIY(FORMUL,MESH%M%X%R,
     *              MESH%M%TYPDIA,MESH%M%TYPEXT,
     *              XMUL,F,G,H,U,V,W,
     *              F%R,G%R,H%R,U%R,V%R,W%R,
     *              MESH%W%R,LEGO,XX,YY,ZZ,
     *              SURFAC,IKLE,MESH%NBOR%I ,
     *              NELEM,NELMAX,IELM1,IELM2,SS)
      ENDIF
C
C  MISE A JOUR DES INFORMATIONS DE LA MATRICE
C    
      NPT = NBPTS(MIN(IELM1,IELM2))
C
      IF(NPT.GT.MESH%M%D%MAXDIM1) THEN
        IF (LNG.EQ.1) WRITE(LU,500) MESH%M%NAME
        IF (LNG.EQ.2) WRITE(LU,501) MESH%M%NAME
        IF (LNG.EQ.1) WRITE(LU,2000) IELM1
        IF (LNG.EQ.2) WRITE(LU,2001) IELM1
        IF (LNG.EQ.1) WRITE(LU,3000) IELM2
        IF (LNG.EQ.2) WRITE(LU,3001) IELM2
        CALL PLANTE(1)
        STOP
      ENDIF
C
      MESH%M%D%ELM  = MIN(IELM1,IELM2)
      MESH%M%D%DIM1 = NPT
      MESH%M%ELMLIN = IELM1
      MESH%M%ELMCOL = IELM2
C
C  ASSEMBLAGE EVENTUEL DE LA DIAGONALE
C
      IF(LEGO.AND.MESH%M%TYPDIA.EQ.'Q') THEN
C
            CALL ASSVEC(MESH%M%D%R,
     *                  IKLE,NPT,NELEM,NELMAX,MIN(IELM1,IELM2),
     *                  MESH%W%R,LEGO,MESH%LV,MSK,MASKEL%R)
C
      ENDIF
C
C  MASQUAGE EVENTUEL DES TERMES EXTRA-DIAGONAUX
C
C     NOTE : POUR LE STOCKAGE 2, LES TERMES EXTRA-DIAGONAUX NE SONT
C            PAS A LEUR PLACE MAIS CA NE CHANGE RIEN POUR LA
C            MULTIPLICATION PAR MASKEL.
C
      IF(MSK) CALL OM( 'M=MSK(M)',MESH%M,MESH%M,MASKEL,C,MESH)
C
C  ASSEMBLAGE EVENTUEL DES TERMES EXTRA-DIAGONAUX
C
      IF(M%STO.EQ.3) THEN
C       COPIE DE LA DIAGONALE
        CALL OS('X=Y     ',MESH%MSEG%D,MESH%M%D,MESH%M%D,0.D0)
        MESH%MSEG%TYPDIA(1:1)='Q'
C       ASSEMBLAGE DES TERMES EXTRA-DIAGONAUX
        IF(MESH%M%TYPEXT.EQ.'Q'.OR.MESH%M%TYPEXT.EQ.'S') THEN
C       CASE OF MATRICES WITH INVERTED STORAGE OF OFF-DIAGONAL TERMS
C       SO FAR ONLY MAMURD. SEE INVERSION OF DIM1_EXT AND DIM2_EXT
C       AND 2 INSTEAD OF 1 FOR STOXMT
          IF(FORMUL(1:6).EQ.'MAMURD') THEN
            CALL ASSEX3(MESH%MSEG%X%R,MESH%M%STO,MESH%M%NAME,
     *                  MESH%M%ELMLIN,MESH%M%ELMCOL,
     *                  MESH%M%TYPEXT,MESH%M%X%R,
     *                  DIM2_EXT(IELM1,IELM2,MESH%M%STO,MESH%M%TYPEXT),
     *                  DIM1_EXT(IELM1,IELM2,MESH%M%STO,MESH%M%TYPEXT),
     *                  2,
     *                  MESH,MESH%NELMAX,MESH%ELTSEG%I,MESH%ORISEG%I)
          ELSE
            CALL ASSEX3(MESH%MSEG%X%R,MESH%M%STO,MESH%M%NAME,
     *                  MESH%M%ELMLIN,MESH%M%ELMCOL,
     *                  MESH%M%TYPEXT,MESH%M%X%R,
     *                  DIM1_EXT(IELM1,IELM2,MESH%M%STO,MESH%M%TYPEXT),
     *                  DIM2_EXT(IELM1,IELM2,MESH%M%STO,MESH%M%TYPEXT),
     *                  1,
     *                  MESH,MESH%NELMAX,MESH%ELTSEG%I,MESH%ORISEG%I)        
          ENDIF
        ENDIF
        MESH%MSEG%TYPEXT=MESH%M%TYPEXT
        MESH%MSEG%ELMLIN = IELM1
        MESH%MSEG%ELMCOL = IELM2
        MESH%MSEG%D%ELM  = MIN(IELM1,IELM2)
        MESH%MSEG%D%DIM1 = NPT
        MESH%MSEG%X%DIM1 = DIM1_EXT(IELM1,IELM2,M%STO,MESH%MSEG%TYPEXT)
        MESH%MSEG%X%DIM2 = DIM2_EXT(IELM1,IELM2,M%STO,MESH%MSEG%TYPEXT)
      ENDIF
C
C     DIMENSIONS DU TABLEAU DES TERMES EXTRADIAGONAUX
C     ATTENTION AU M%STO (ON NE MET PAS MESH%M%STO CAR IL VAUT 1)
C                             VOIR DEBUT DU SOUS-PROGRAMME
C
      MESH%M%X%DIM1 = DIM1_EXT(IELM1,IELM2,M%STO,MESH%M%TYPEXT)
      MESH%M%X%DIM2 = DIM2_EXT(IELM1,IELM2,M%STO,MESH%M%TYPEXT)
C
C-----------------------------------------------------------------------
C  APRES TRAVAIL SUR MESH%M, ACTUALISATION FINALE DE M
C-----------------------------------------------------------------------
C
      IF(M%STO.EQ.1) THEN
        CALL OM( OP , M , MESH%M , F , C , MESH )
      ELSEIF(M%STO.EQ.3) THEN
        CALL OM( OP , M , MESH%MSEG , F , C , MESH )
      ELSE
        WRITE(LU,*) 'MATRIX, STOCKAGE INCONNU : ',M%STO
        STOP        
      ENDIF
C
C-----------------------------------------------------------------------
C
500   FORMAT(1X,'MATRIX (BIEF) : MATRICE ',A6,' TROP PETITE')
501   FORMAT(1X,'MATRIX (BIEF) : MATRIX ',A6,' TOO SMALL')
2000  FORMAT(1X,'                POUR IELM1 = ',1I6)
2001  FORMAT(1X,'                FOR IELM1 = ',1I6)
3000  FORMAT(1X,'                ET IELM2 = ',1I6)
3001  FORMAT(1X,'                AND IELM2 = ',1I6)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
