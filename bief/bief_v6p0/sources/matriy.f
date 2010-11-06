C                       *****************
                        SUBROUTINE MATRIY
C                       *****************
C
     *(FORMUL,XM,TYPDIA,TYPEXT,
     * XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,T,LEGO,
     * XEL,YEL,ZEL,SURFAC,IKLE,NBOR,
     * NELEM,NELMAX,IELM1,IELM2,S)
C
C***********************************************************************
C BIEF VERSION 5.9      21/07/2008     JM HERVOUET (LNHE) 01 30 87 80 18
C
C 16/09/2005: MT04PP NOW IN 2 OPTIONS, AND WITHOUT VERTICAL UPWIND
C
C***********************************************************************
C 
C FONCTION : CALCULS DE MATRICES
C
C             LA MATRICE EST IDENTIFIEE PAR LA FORMULE CONTENUE DANS
C             LA CHAINE DE CARACTERES FORMUL
C
C-----------------------------------------------------------------------
C
C  SIGNIFICATION DE IELM ET IELM2
C
C  TYPE D'ELEMENT      NOMBRE DE POINTS
C
C  A OU 11 : TRIANGLE P1            3
C  B OU 12 : TRIANGLE QUASI-BULLE   4
C  C OU 13 : TRIANGLE P1-ISO P1     6
C  D OU 14 : TRIANGLE P2            7
C  E       : RIEN POUR L'INSTANT
C
C  F OU 21 : QUADRILATERE Q1        4
C  G OU 22 : QUADRILATERE Q2        8
C  H OU 24 : QUADRILATERE Q2        9
C
C  T OU 31 : TETRAEDRE P1           4
C
C  P OU 41 : PRISMES TELEMAC-3D     6
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  FORMUL        | -->|  FORMULE DECRIVANT LA MATRICE
C |  XM            |<-->|  TERMES EXTRA-DIAGONAUX
C |  XMUL          | -->|  COEFFICIENT MULTIPLICATEUR DU RESULTAT
C |  SF,SG,SH      | -->|  STRUCTURES DES FONCTIONS F,G,H ET
C |  SU,SV,SW      | -->|  U,V,H CI-DESSOUS.
C |  F,G,H         | -->|  FONCTIONS INTERVENANT DANS LA FORMULE
C |  U,V,W         | -->|  COMPOSANTES D'UN VECTEUR U DANS LA FORMULE
C |  T             | -->|  TABLEAU DE TRAVAIL DE DIMENSION QUI
C |                |    |  CONTIENDRA LA DIAGONALE NON ASSEMBLEE
C |  LEGO          | -->|  LOGIQUE : POUR ASSEMBLER LA DIAGONALE
C |  XEL,YEL,ZEL   | -->|  COORDONNEES DES POINTS DANS L'ELEMENT
C |  SURFAC        | -->|  SURFACE DES ELEMENTS.
C |  IKLE          | -->|  PASSAGE DE LA NUMEROTATION LOCALE A GLOBALE
C |  NELEM         | -->|  NOMBRE D'ELEMENTS DU MAILLAGE
C |  NELMAX        | -->|  NOMBRE MAXIMUM D'ELEMENTS DU MAILLAGE
C |                |    |  (CAS DU MAILLAGE ADAPTATIF)
C |  IELM1         | -->|  TYPE D'ELEMENT POUR LES LIGNES
C |  IELM2         | -->|  TYPE D'ELEMENT POUR LES COLONNES
C |  S             | -->|  TYPE DE STOCKAGE DE LA MATRICE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : MATRIX
C PROGRAMMES APPELES : TOUS LES SOUS-PROGRAMMES DE MATRICES ELEMENTAIRES
C                      ASSVEC
C
C**********************************************************************
C
      USE BIEF, EX_MATRIY => MATRIY
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NELMAX,NELEM,IELM1,IELM2,S
      INTEGER, INTENT(IN)             :: IKLE(NELMAX,*),NBOR(*)
      LOGICAL, INTENT(INOUT)          :: LEGO
      TYPE(BIEF_OBJ), INTENT(IN)      :: SF,SG,SH,SU,SV,SW
      DOUBLE PRECISION, INTENT(IN)    :: F(*),G(*),H(*),U(*),V(*),W(*)
      DOUBLE PRECISION, INTENT(IN)    :: SURFAC(NELMAX)
      DOUBLE PRECISION, INTENT(IN)    :: XEL(*),YEL(*),ZEL(*)
      DOUBLE PRECISION, INTENT(IN)    :: XMUL 
      DOUBLE PRECISION, INTENT(INOUT) :: XM(NELMAX,*),T(NELMAX,*)
      CHARACTER(LEN=16), INTENT(IN)   :: FORMUL 
      CHARACTER(LEN=1), INTENT(INOUT) :: TYPDIA,TYPEXT    
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+     
C           
      LOGICAL SIGMAG,INCHYD,SPECAD
C
      CHARACTER(LEN=1)  TDIA,TEXT      
C
C-----------------------------------------------------------------------
C
      INTEGER ICOORD
C
      INTEGER AAQ(3,3,2),BBQ(4,4,2),ABQ(3,4,2),BAQ(4,3,2),PPQ(6,6,2)
      INTEGER AAS(3,3,2),BBS(4,4,2),OOS(2,2,2),FFS(4,4,2),PPS(6,6,2)
      INTEGER ACQ(3,6,2),CAQ(6,3,2)
C
C
C  ATTENTION LES MATRICES CI-DESSOUS SONT A TRANSPOSER
C  DANS LE CAS QUELCONQUE, A CAUSE DE LA NOTATION FORTRAN DES DATA
C
C  ATTENTION : OM N'A PAS ETE PARAMETRE AVEC CES TABLEAUX
C              CES DATA FIGURENT AUSSI DANS MATVCT.
C
      DATA OOS/  0 ,  1 ,
     *           1 ,  0 ,
C S=2 NON PREVU
     *           0 ,  0 ,
     *           0 ,  0 /
C
C     P1-P1 EBE SYMETRIQUE (S=1)
      DATA AAS/  0 ,  1 ,  2 ,
     *           1 ,  0 ,  3 ,
     *           2 ,  3 ,  0 ,
C     P1-P1 EBE PRE-ASSEMBLE SYMETRIQUE (S=2)
     *           0 ,  1 ,  3 ,
     *           1 ,  0 ,  2 ,
     *           3 ,  2 ,  0 /
C
C     P1-P1 EBE NON SYMETRIQUE (S=1)
      DATA AAQ/  0 ,  4 ,  5 ,
     *           1 ,  0 ,  6 ,
     *           2 ,  3 ,  0 ,
C     P1-P1 EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
     *           0 ,  4 ,  3 ,
     *           1 ,  0 ,  5 ,
     *           6 ,  2 ,  0 /
C
C     QUASI-BULLE QUASI-BULLE EBE SYMETRIQUE (S=1)
      DATA BBS/  0 ,  1 ,  2 ,  3 ,
     *           1 ,  0 ,  4 ,  5 ,
     *           2 ,  4 ,  0 ,  6 ,
     *           3 ,  5 ,  6 ,  0 ,
C     QUASI-BULLE QUASI-BULLE EBE PRE-ASSEMBLE SYMETRIQUE (S=2)
     *           0 ,  4 ,  6 ,  1 ,
     *           4 ,  0 ,  5 ,  2 ,
     *           6 ,  5 ,  0 ,  3 ,
     *           1 ,  2 ,  3 ,  0 /
C
C     QUASI-BULLE QUASI-BULLE EBE NON SYMETRIQUE (S=1)
      DATA BBQ/  0 ,  7 ,  8 ,  9 ,
     *           1 ,  0 , 10 , 11 ,
     *           2 ,  4 ,  0 , 12 ,
     *           3 ,  5 ,  6 ,  0 ,
C     QUASI-BULLE QUASI-BULLE EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
     *           0 , 10 ,  6 ,  7 ,
     *           4 ,  0 , 11 ,  8 ,
     *          12 ,  5 ,  0 ,  9 ,
     *           1 ,  2 ,  3 ,  0 /
C
C     P1 QUASI-BULLE EBE NON SYMETRIQUE (S=1)
      DATA ABQ/  0 ,  4 ,  7 ,
     *           1 ,  0 ,  8 ,
     *           2 ,  5 ,  0 ,
     *           3 ,  6 ,  9 ,
C     P1 QUASI-BULLE EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
     *           0 ,  7 ,  3 ,
     *           1 ,  0 ,  8 ,
     *           9 ,  2 ,  0 ,
     *           4 ,  5 ,  6 /
C
C     QUASI-BULLE P1 EBE NON SYMETRIQUE (S=1)
      DATA BAQ/  0 ,  3 ,  5 ,  7 ,
     *           1 ,  0 ,  6 ,  8 ,
     *           2 ,  4 ,  0 ,  9 ,
C     QUASI-BULLE P1 EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
     *           0 ,  7 ,  3 ,  4 ,
     *           1 ,  0 ,  8 ,  5 ,
     *           9 ,  2 ,  0 ,  6 /
C     P1 P2 EBE NON SYMETRIQUE (S=1)
      DATA ACQ/   0 ,  6  , 11 ,
     *            1 ,  0  , 12 ,
     *            2 ,  7  ,  0 ,
     *            3 ,  8  , 13 ,
     *            4 ,  9  , 14 ,         
     *            5 ,  10 , 15 , 
C S=2 NON PREVU  
     *            0 ,  0 ,  0  ,
     *            0 ,  0 ,  0  ,
     *            0 ,  0 ,  0  ,
     *            0 ,  0 ,  0  ,
     *            0 ,  0 ,  0  ,         
     *            0 ,  0 ,  0  /
C     P2 P1 EBE NON SYMETRIQUE (S=1)
      DATA CAQ/  0 ,  3 ,  5 ,  7 , 10 , 13 ,
     *           1 ,  0 ,  6 ,  8 , 11 , 14 ,
     *           2 ,  4 ,  0 ,  9 , 12 , 15 ,
C     P2 P1 EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
C     -NON PREVU
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     AJOUTE PAR JMJ MAIS PAS UTILISE
C     PRISMES P1-P1 ET TRIANGLES P2 EBE SYMETRIQUE (S=1)
      DATA PPS/  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,
     *           1 ,  0 ,  6 ,  7 ,  8 ,  9 ,
     *           2 ,  6 ,  0 , 10 , 11 , 12 ,
     *           3 ,  7 , 10 ,  0 , 13 , 14 ,
     *           4 ,  8 , 11 , 13 ,  0 , 15 ,
     *           5 ,  9 , 12 , 14 , 15 ,  0 ,
C     PRISMES P1-P1 EBE PRE-ASSEMBLE SYMETRIQUE (S=2) - NON PREVU
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     PRISMES P1-P1 EBE NON SYMETRIQUE (S=1)
      DATA PPQ/  0 , 16 , 17 , 18 , 19 , 20 ,
     *           1 ,  0 , 21 , 22 , 23 , 24 ,
     *           2 ,  6 ,  0 , 25 , 26 , 27 ,
     *           3 ,  7 , 10 ,  0 , 28 , 29 ,
     *           4 ,  8 , 11 , 13 ,  0 , 30 ,
     *           5 ,  9 , 12 , 14 , 15 ,  0 ,
C     PRISMES P1-P1 EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2) - NON PREVU
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     QUADRANGLES Q1-Q1 EBE SYMETRIQUE (S=1)
      DATA FFS/  0 ,  1 ,  2 ,  3 ,
     *           1 ,  0 ,  4 ,  5 ,
     *           2 ,  4 ,  0 ,  6 ,
     *           3 ,  5 ,  6 ,  0 ,
C     QUADRANGLES Q1-Q1 EBE PRE-ASSEMBLE SYMETRIQUE (S=2) - NON PREVU
     *           0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 /
C
C-----------------------------------------------------------------------
C
C  TEST DU TYPE DE MATRICE
C
C=======================================================================
C     MATRICE DE MASSE
C=======================================================================
C
      IF(FORMUL(1:16).EQ.'MATMAS          ') THEN
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
             CALL MT01AA(   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                                       T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                                        T(1,3)   ,
     *                   XMUL,SURFAC,NELEM,NELMAX)
C
            TYPDIA='Q'
            TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C.......................................................................
C         ERREUR SUR L'ELEMENT DE COLONNE
C.......................................................................
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
          CALL MT01BB
     * (   T(1,1)   ,XM(1,BBS(1,2,S)),XM(1,BBS(1,3,S)),XM(1,BBS(1,4,S)),
     *                      T(1,2)   ,XM(1,BBS(2,3,S)),XM(1,BBS(2,4,S)),
     *                                       T(1,3)   ,XM(1,BBS(3,4,S)),
     *                                                        T(1,4)   ,
     *      XMUL,SURFAC,NELEM,NELMAX)
C
            TYPDIA='Q'
            TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C.......................................................................
C         ERREUR SUR L'ELEMENT DE COLONNE
C.......................................................................
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       TRIANGLE P2
C-----------------------------------------------------------------------
C
        ELSEIF(IELM1.EQ.13) THEN
C
C       TEST DE L'ELEMENT DE COLONNE
C
C.......................................................................
C         TRIANGLE P2
C.......................................................................
C
          IF(IELM2.EQ.13) THEN
             CALL MT01CC
     * (   T(1,1)   ,XM(1,PPS(1,2,S)),XM(1,PPS(1,3,S)),
     *   XM(1,PPS(1,4,S)),XM(1,PPS(1,5,S)),XM(1,PPS(1,6,S)),
     *     T(1,2)   ,XM(1,PPS(2,3,S)),XM(1,PPS(2,4,S)),
     *   XM(1,PPS(2,5,S)),XM(1,PPS(2,6,S)),
     *     T(1,3)   ,XM(1,PPS(3,4,S)),XM(1,PPS(3,5,S)),
     *   XM(1,PPS(3,6,S)),
     *     T(1,4)   ,XM(1,PPS(4,5,S)),XM(1,PPS(4,6,S)),
     *     T(1,5)   ,XM(1,PPS(5,6,S)),T(1,6),
     *     XMUL,SURFAC,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE SEGMENT P1
        ELSEIF(IELM1.EQ.1) THEN
C.......................................................................
C         ELEMENT DE COLONNE SEGMENT P1
          IF(IELM2.EQ.1.AND.S.EQ.1) THEN
             CALL MT01OO(   T(1,1)   ,XM(1,OOS(1,2,S)),
     *                                       T(1,2)   ,
     *                   XMUL,SURFAC,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE ELEMENT DE COLONNE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C>>>>
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
        ELSEIF(IELM1.EQ.41) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.41) THEN
             CALL MT01PP(T,XM,XMUL,ZEL,SURFAC,IKLE,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TETRAEDRE T1
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.31.OR.IELM2.EQ.51) THEN
             CALL MT01TT(T,XM,XMUL,XEL,YEL,ZEL,IKLE,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE DE DIFFUSION
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'MATDIF') THEN
C
C     LE CARACTERE 7 DIT SI INCHYD OU NON
C
      INCHYD = .FALSE.
      IF(FORMUL(7:7).EQ.'2') INCHYD = .TRUE.
C
C     TEST DE L'ELEMENT DE LIGNE
C
C-----------------------------------------------------------------------
C       TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
C
C     LE CARACTERE 7 DIT AUSSI SI ON VEUT LE TERME
C     DE DIFFUSION POUR ESTEL.
C
      IF(FORMUL(7:7).NE.'3') THEN
C
             CALL MT02AA(  T(1,1),XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                                   T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                                      T(1,3) ,
     *                         XMUL,SU,U,SV,V,XEL,YEL,SURFAC,
     *             IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX,FORMUL)
C
      ELSE
C
             CALL MT02AA_2( T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                                       T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                                        T(1,3)   ,
     *                   XMUL,SU,SV,U,V,XEL,YEL,SURFAC,NELEM,NELMAX)

C
      ENDIF
C
      TYPDIA='Q'
      TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT02BB
     * (   T(1,1)   ,XM(1,BBS(1,2,S)),XM(1,BBS(1,3,S)),XM(1,BBS(1,4,S)),
     *                      T(1,2)   ,XM(1,BBS(2,3,S)),XM(1,BBS(2,4,S)),
     *                                       T(1,3)   ,XM(1,BBS(3,4,S)),
     *                                                        T(1,4)   ,
     *          XMUL,SU,U,XEL,YEL,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P2
C      
        ELSEIF(IELM1.EQ.13) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P2
C         UTILISATION DE LA MATRICE PPS car ICI AUSSI ON AURA UNE MATRICE 
C         6x6 SYMETRIQUE         

          IF(IELM2.EQ.13) THEN
             CALL MT02CC
     * (   T(1,1)   ,XM(1,PPS(1,2,S)),XM(1,PPS(1,3,S)),
     *   XM(1,PPS(1,4,S)),XM(1,PPS(1,5,S)),XM(1,PPS(1,6,S)),
     *     T(1,2)   ,XM(1,PPS(2,3,S)),XM(1,PPS(2,4,S)),
     *   XM(1,PPS(2,5,S)),XM(1,PPS(2,6,S)),
     *     T(1,3)   ,XM(1,PPS(3,4,S)),XM(1,PPS(3,5,S)),
     *   XM(1,PPS(3,6,S)),
     *     T(1,4)   ,XM(1,PPS(4,5,S)),XM(1,PPS(4,6,S)),
     *     T(1,5)   ,XM(1,PPS(5,6,S)),T(1,6),
     *          XMUL,SU,U,XEL,YEL,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
        ELSEIF(IELM1.EQ.41) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.41) THEN
C
C         CALL MT02PT(T,XM,XMUL,SF,SG,SH,F,G,H,
C    *                XEL,YEL,ZEL,IKLE,NELEM,NELMAX,INCHYD)
          CALL MT02PP(T,XM,XMUL,SF,SG,SH,F,G,H,
     *                XEL,YEL,ZEL,SURFAC,IKLE,NELEM,NELMAX,INCHYD,
     *                FORMUL)
C
          TYPDIA='Q'
          TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TETRAEDRE
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.31.OR.IELM2.EQ.51) THEN
C
          CALL MT02TT(T,XM,XMUL,SF,SG,SH,F,G,H,
     *                XEL,YEL,ZEL,IKLE,NELEM,NELMAX,INCHYD)
C
          TYPDIA='Q'
          TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     CONTRIBUTION DE SUPG A LA MATRICE DE MASSE
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'MASUPG          ') THEN
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
             CALL MT03AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,SF,SG,SU,SV,F,G,U,V,XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                   NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT03BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *          XMUL,SF,SG,SU,SV,F,G,U,V,
     *          XEL,YEL,IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P2
        ELSEIF(IELM1.EQ.13) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P2
          IF(IELM2.EQ.13) THEN
             CALL MT03CC
     * (   T(1,1)   ,XM(1,PPQ(1,2,S)),XM(1,PPQ(1,3,S)),
     *   XM(1,PPQ(1,4,S)),XM(1,PPQ(1,5,S)),XM(1,PPQ(1,6,S)),
     *   XM(1,PPQ(2,1,S)),T(1,2),XM(1,PPQ(2,3,S)),XM(1,PPQ(2,4,S)),
     *   XM(1,PPQ(2,5,S)),XM(1,PPQ(2,6,S)),XM(1,PPQ(3,1,S)),
     *   XM(1,PPQ(3,2,S)),T(1,3),XM(1,PPQ(3,4,S)),XM(1,PPQ(3,5,S)),
     *   XM(1,PPQ(3,6,S)),XM(1,PPQ(4,1,S)),XM(1,PPQ(4,2,S)),
     *   XM(1,PPQ(4,3,S)),T(1,4),XM(1,PPQ(4,5,S)),XM(1,PPQ(4,6,S)),
     *   XM(1,PPQ(5,1,S)),XM(1,PPQ(5,2,S)),XM(1,PPQ(5,3,S)),
     *   XM(1,PPQ(5,4,S)),  T(1,5)   ,XM(1,PPQ(5,6,S)),
     *   XM(1,PPQ(6,1,S)),XM(1,PPQ(6,2,S)) ,XM(1,PPQ(6,3,S)),
     *   XM(1,PPQ(6,4,S)),XM(1,PPQ(6,5,S)),T(1,6),
     *          XMUL,SF,SG,SU,SV,F,G,U,V,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *          NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='Q'
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE U.GRAD U.GRAD
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'MAUGUG') THEN
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
             CALL MT04AA(   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                                       T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                                        T(1,3)   ,
     *                   XMUL,SU,SV,U,V,XEL,YEL,SURFAC,IKLE,
     *                   NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT04BB
     * (   T(1,1)   ,XM(1,BBS(1,2,S)),XM(1,BBS(1,3,S)),XM(1,BBS(1,4,S)),
     *                      T(1,2)   ,XM(1,BBS(2,3,S)),XM(1,BBS(2,4,S)),
     *                                       T(1,3)   ,XM(1,BBS(3,4,S)),
     *                                                        T(1,4)   ,
     *          XMUL,SU,SV,U,V,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P2
        ELSEIF(IELM1.EQ.13) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P2
          IF(IELM2.EQ.13) THEN
             CALL MT04CC
     * (   T(1,1)   ,XM(1,PPS(1,2,S)),XM(1,PPS(1,3,S)),
     *   XM(1,PPS(1,4,S)),XM(1,PPS(1,5,S)),XM(1,PPS(1,6,S)),
     *     T(1,2)   ,XM(1,PPS(2,3,S)),XM(1,PPS(2,4,S)),
     *   XM(1,PPS(2,5,S)),XM(1,PPS(2,6,S)),
     *     T(1,3)   ,XM(1,PPS(3,4,S)),XM(1,PPS(3,5,S)),
     *   XM(1,PPS(3,6,S)),
     *     T(1,4)   ,XM(1,PPS(4,5,S)),XM(1,PPS(4,6,S)),
     *     T(1,5)   ,XM(1,PPS(5,6,S)),T(1,6),
     *          XMUL,SU,SV,U,V,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          IKLE(1,4),IKLE(1,5),IKLE(1,6),NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF          
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
        ELSEIF(IELM1.EQ.41) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.41) THEN
             CALL MT04PP(T,XM,XMUL,SU,SV,SW,U,V,W,
     *                   XEL,YEL,ZEL,SURFAC,IKLE,NELEM,NELMAX,FORMUL)
C           
             TYPDIA='Q'
             IF(FORMUL(7:7).EQ.'2') THEN
               TYPEXT='S'
             ELSE
               TYPEXT='Q'
             ENDIF
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TETRAEDRE
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.31.OR.IELM2.EQ.51) THEN
             CALL MT04TT(T,XM,XMUL,SU,SV,SW,U,V,W,
     *                   XEL,YEL,ZEL,IKLE,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE U.GRAD
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'MATVGR') THEN
C
      SIGMAG = .FALSE.
      IF(FORMUL(7:7).EQ.'2') SIGMAG = .TRUE.
      SPECAD = .FALSE.
      IF(FORMUL(8:8).EQ.'2') SPECAD = .TRUE.
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
             CALL MT05AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,SU,SV,U,V,XEL,YEL,IKLE,
     *                   NELEM,NELMAX,FORMUL)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT05BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *          XMUL,SU,SV,U,V,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,FORMUL)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P2
        ELSEIF(IELM1.EQ.13) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P2
          IF(IELM2.EQ.13) THEN
             CALL MT05CC
     * (   T(1,1)   ,XM(1,PPQ(1,2,S)),XM(1,PPQ(1,3,S)),
     *   XM(1,PPQ(1,4,S)),XM(1,PPQ(1,5,S)),XM(1,PPQ(1,6,S)),
     *   XM(1,PPQ(2,1,S)),T(1,2),XM(1,PPQ(2,3,S)),XM(1,PPQ(2,4,S)),
     *   XM(1,PPQ(2,5,S)),XM(1,PPQ(2,6,S)),XM(1,PPQ(3,1,S)),
     *   XM(1,PPQ(3,2,S)),T(1,3),XM(1,PPQ(3,4,S)),XM(1,PPQ(3,5,S)),
     *   XM(1,PPQ(3,6,S)),XM(1,PPQ(4,1,S)),XM(1,PPQ(4,2,S)),
     *   XM(1,PPQ(4,3,S)),T(1,4),XM(1,PPQ(4,5,S)),XM(1,PPQ(4,6,S)),
     *   XM(1,PPQ(5,1,S)),XM(1,PPQ(5,2,S)),XM(1,PPQ(5,3,S)),
     *   XM(1,PPQ(5,4,S)),  T(1,5)   ,XM(1,PPQ(5,6,S)),
     *   XM(1,PPQ(6,1,S)),XM(1,PPQ(6,2,S)) ,XM(1,PPQ(6,3,S)),
     *   XM(1,PPQ(6,4,S)),XM(1,PPQ(6,5,S)),T(1,6),
     *          XMUL,SU,SV,U,V,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *          NELEM,NELMAX,FORMUL)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
        ELSEIF(IELM1.EQ.41) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.41) THEN
             CALL MT05PP(T,XM,XMUL,SU,SV,SW,U,V,W,SF,SG,SH,F,G,H,
     *                   XEL,YEL,ZEL,IKLE,NELEM,NELMAX,SURFAC,SIGMAG,
     *                   SPECAD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TETRAEDRE
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.31.OR.IELM2.EQ.51) THEN
             CALL MT05TT(T,XM,XMUL,SU,SV,SW,U,V,W,
     *                   XEL,YEL,ZEL,IKLE,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE F PSI PSJ
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'FMATMA') THEN
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
C
             CALL MT06AA(   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                                       T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                                        T(1,3)   ,
     *                   XMUL,SF,F,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),NELEM,NELMAX)          
C
            TYPDIA='Q'
            TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1 DE FACE LATERALE DE PRISMES DECOUPES
C       EN TETRAEDRES
C
        ELSEIF(IELM1.EQ.61.OR.IELM1.EQ.81) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.61.OR.IELM1.EQ.81) THEN
C
C     LE CARACTERE 7 DIT AUSSI SI ON VEUT LE TERME
C     DE PLUS POUR ESTEL-3D.
C
      IF(FORMUL(7:7).NE.'2') THEN
C
                CALL MT06FT
     * (   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                      T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                       T(1,3)   ,
     *          XMUL,SF,F,XEL,YEL,ZEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          NBOR,NELEM,NELMAX)
C
      ELSE
                CALL MT06FT2
     * (   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                      T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                       T(1,3)   ,
     *          XMUL,SF,F,SU,U,XEL,YEL,ZEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          NBOR,NELEM,NELMAX)
      ENDIF
            TYPDIA='Q'
            TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT06BB
     * (   T(1,1)   ,XM(1,BBS(1,2,S)),XM(1,BBS(1,3,S)),XM(1,BBS(1,4,S)),
     *                      T(1,2)   ,XM(1,BBS(2,3,S)),XM(1,BBS(2,4,S)),
     *                                       T(1,3)   ,XM(1,BBS(3,4,S)),
     *                                                        T(1,4)   ,
     *          XMUL,SF,F,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUADRATIQUE
C
        ELSEIF(IELM1.EQ.13) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUADRATIQUE
          IF(IELM2.EQ.13) THEN
             CALL MT06CC
     * ( T(1,1)   ,XM(1,PPS(1,2,S)),XM(1,PPS(1,3,S)),
     *   XM(1,PPS(1,4,S)),XM(1,PPS(1,5,S)),XM(1,PPS(1,6,S)),
     *   T(1,2)   ,XM(1,PPS(2,3,S)),XM(1,PPS(2,4,S)),
     *   XM(1,PPS(2,5,S)),XM(1,PPS(2,6,S)),
     *   T(1,3)   ,XM(1,PPS(3,4,S)),XM(1,PPS(3,5,S)),XM(1,PPS(3,6,S)),
     *   T(1,4)   ,XM(1,PPS(4,5,S)),XM(1,PPS(4,6,S)),
     *   T(1,5)   ,XM(1,PPS(5,6,S)),T(1,6)          ,
     *          XMUL,SF,F,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          IKLE(1,4),IKLE(1,5),IKLE(1,6),NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF          
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE QUADRILATERE DE FACE LATERALE DE PRISME
        ELSEIF(IELM1.EQ.71) THEN
C
C.......................................................................
C
C         ATTENTION !!!!!!!!!!
C         ELEMENT DE COLONNE QUADRANGLE POUR LES PRISMES DE TELEMAC-3D
          IF(IELM2.EQ.71) THEN
C
               CALL MT06FF
     * (   T(1,1)   ,XM(1,FFS(1,2,S)),XM(1,FFS(1,3,S)),XM(1,FFS(1,4,S)),
     *                      T(1,2)   ,XM(1,FFS(2,3,S)),XM(1,FFS(2,4,S)),
     *                                       T(1,3)   ,XM(1,FFS(3,4,S)),
     *                                                        T(1,4)   ,
     *         XMUL,SF,F,XEL,YEL,ZEL,
     *         IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *         NBOR,NELEM,NELMAX)
C
               TYPDIA='Q'
               TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C<<<<
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE SEGMENT P1
        ELSEIF(IELM1.EQ.1) THEN
C.......................................................................
C         ELEMENT DE COLONNE SEGMENT P1
          IF(IELM2.EQ.1.AND.S.EQ.1) THEN
             CALL MT06OO(   T(1,1)   ,XM(1,OOS(1,2,S)),
     *                                       T(1,2)   ,
     *                   XMUL,SF,F,SURFAC,IKLE(1,1),IKLE(1,2),
     *                   NBOR,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C<<<<
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE SEGMENT P2
        ELSEIF(IELM1.EQ.2) THEN
C.......................................................................
C         ELEMENT DE COLONNE SEGMENT P2
          IF(IELM2.EQ.2.AND.S.EQ.1) THEN
               CALL MT06OC
     * (   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                      T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                       T(1,3)   ,
     *        XMUL,SF,F,SURFAC,IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NBOR,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
C
        ELSE IF (IELM1.EQ.41) THEN
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF (IELM2.EQ.41) THEN
             CALL MT06PP(T,XM,
     *                   XMUL,SF,F,ZEL,SURFAC,IKLE,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
c             CALL SETDIA(M,'Q')
c             CALL SETEXT(M,'S')
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TETRAEDRE
C
        ELSE IF (IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF (IELM2.EQ.31.OR.IELM2.EQ.51) THEN
             CALL MT06TT(T,XM,
     *                   XMUL,SF,F,
     *                   XEL,YEL,ZEL,IKLE,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
c             CALL SETDIA(M,'Q')
c             CALL SETEXT(M,'S')
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE DE MASSE MASS-LUMPEE
C=======================================================================
C
      ELSEIF(FORMUL(1:16).EQ.'MSLUMP          ') THEN
C
C     TEST DE L'ELEMENT DE LIGNE
C
C-----------------------------------------------------------------------
C       TRIANGLE P1
C-----------------------------------------------------------------------
C
        IF(IELM1.EQ.11) THEN
C
C       TEST DE L'ELEMENT DE COLONNE
C
C.......................................................................
C         TRIANGLE P1
C.......................................................................
C
          IF(IELM2.EQ.11) THEN
             CALL MT07AA(   T(1,1)   ,XM(1,AAS(1,2,S)),XM(1,AAS(1,3,S)),
     *                                       T(1,2)   ,XM(1,AAS(2,3,S)),
     *                                                        T(1,3)   ,
     *                   XMUL,SF,F,SURFAC,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       TRIANGLE QUASI-BULLE
C-----------------------------------------------------------------------
C
        ELSEIF(IELM1.EQ.12) THEN
C
C       TEST DE L'ELEMENT DE COLONNE
C
C.......................................................................
C         TRIANGLE QUASI-BULLE
C.......................................................................
C
          IF(IELM2.EQ.12) THEN
             CALL MT07BB
     * (   T(1,1)   ,XM(1,BBS(1,2,S)),XM(1,BBS(1,3,S)),XM(1,BBS(1,4,S)),
     *                      T(1,2)   ,XM(1,BBS(2,3,S)),XM(1,BBS(2,4,S)),
     *                                       T(1,3)   ,XM(1,BBS(3,4,S)),
     *                                                        T(1,4)   ,
     *          XMUL,SF,F,SURFAC,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       TRIANGLE P2
C-----------------------------------------------------------------------
C
        ELSEIF(IELM1.EQ.13) THEN
C
C       TEST DE L'ELEMENT DE COLONNE
C
C.......................................................................
C         TRIANGLE P2
C.......................................................................
C
          IF(IELM2.EQ.13) THEN
             CALL MT07CC
     * (   T(1,1)   ,XM(1,PPS(1,2,S)),XM(1,PPS(1,3,S)),
     *   XM(1,PPS(1,4,S)),XM(1,PPS(1,5,S)),XM(1,PPS(1,6,S)),
     *     T(1,2)   ,XM(1,PPS(2,3,S)),XM(1,PPS(2,4,S)),
     *   XM(1,PPS(2,5,S)),XM(1,PPS(2,6,S)),
     *     T(1,3)   ,XM(1,PPS(3,4,S)),XM(1,PPS(3,5,S)),
     *   XM(1,PPS(3,6,S)),
     *     T(1,4)   ,XM(1,PPS(4,5,S)),XM(1,PPS(4,6,S)),
     *     T(1,5)   ,XM(1,PPS(5,6,S)),T(1,6),
     *          XMUL,SF,F,SURFAC,NELEM,NELMAX)
C
             TYPDIA='Q'
             TYPEXT='S'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE U GRADIENT
C=======================================================================
C
      ELSEIF(FORMUL(1:15).EQ.'MATFGR         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
      IF(FORMUL(16:16).EQ.'X') THEN
        ICOORD=1
      ELSEIF(FORMUL(16:16).EQ.'Y') THEN
        ICOORD=2
      ELSEIF(FORMUL(16:16).EQ.'Z') THEN
        ICOORD=3
      ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
C-----------------------------------------------------------------------
C
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
C.......................................................................
C
          IF(IELM2.EQ.11) THEN
             CALL MT08AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,SF,F,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C
          ELSEIF(IELM2.EQ.12) THEN
             CALL MT08AB
     * (   T(1,1)   ,XM(1,ABQ(1,2,S)),XM(1,ABQ(1,3,S)),XM(1,ABQ(1,4,S)),
     *  XM(1,ABQ(2,1,S)),   T(1,2)   ,XM(1,ABQ(2,3,S)),XM(1,ABQ(2,4,S)),
     *  XM(1,ABQ(3,1,S)),XM(1,ABQ(3,2,S)),   T(1,3)   ,XM(1,ABQ(3,4,S)),
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C 
           ELSEIF(IELM2.EQ.13) THEN
             CALL MT08AC(   T(1,1)   ,XM(1,ACQ(1,2,S)),XM(1,ACQ(1,3,S)),
     *                   XM(1,ACQ(1,4,S)),XM(1,ACQ(1,5,S)),
     *                   XM(1,ACQ(1,6,S)),XM(1,ACQ(2,1,S)),
     *                   T(1,2)  ,XM(1,ACQ(2,3,S)),
     *                   XM(1,ACQ(2,4,S)),XM(1,ACQ(2,5,S)),
     *                   XM(1,ACQ(2,6,S)),XM(1,ACQ(3,1,S)),
     *                   XM(1,ACQ(3,2,S)),T(1,3)  ,XM(1,ACQ(3,4,S)),
     *                   XM(1,ACQ(3,5,S)),XM(1,ACQ(3,6,S)),
     *                   XMUL,SF,F,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *                   NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT08BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C         ELEMENT DE COLONNE TRIANGLE LINEAIRE
          ELSEIF(IELM2.EQ.11) THEN
             CALL MT08BA
     *         (     T(1,1)     ,XM(1,BAQ(1,2,S)),XM(1,BAQ(1,3,S)),
     *          XM(1,BAQ(2,1,S)),     T(1,2)     ,XM(1,BAQ(2,3,S)),
     *          XM(1,BAQ(3,1,S)),XM(1,BAQ(3,2,S)),     T(1,3)     ,
     *          XM(1,BAQ(4,1,S)),XM(1,BAQ(4,2,S)),XM(1,BAQ(4,3,S)),
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C
C       ELEMENT DE LIGNE PRISME P1
        ELSEIF(IELM1.EQ.41) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.41.AND.ICOORD.EQ.3) THEN
             CALL MT08PP(T,XM,XMUL,SF,F,SURFAC,IKLE,NELEM,NELMAX)
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
        ELSEIF(IELM1.EQ.31.OR.IELM1.EQ.51) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.31.OR.IELM2.EQ.51) THEN
             CALL MT08TT(T,XM,XMUL,XEL,YEL,ZEL,SF,F,IKLE,NELEM,NELMAX)
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C     MATRICE F U GRADIENT (NON PROGRAMMEE)
C=======================================================================
C
C     ELSEIF(FORMUL(1:15).EQ.'MATQGR         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
C     IF(FORMUL(16:16).EQ.'X') THEN
C       ICOORD=1
C     ELSEIF(FORMUL(16:16).EQ.'Y') THEN
C       ICOORD=2
C     ELSEIF(FORMUL(16:16).EQ.'Z') THEN
C       ICOORD=3
C     ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
C-----------------------------------------------------------------------
C
C       IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
C.......................................................................
C
C         IF(IELM1.EQ.11) THEN
C            CALL MT09AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
C    *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
C    *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
C    *                   XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
C    *                   NELEM,NELMAX,ICOORD)
C
C            TYPDIA='Q'
C            TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
C         ELSE
C           IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
C           IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
C           IF (LNG.EQ.1) WRITE(LU,2000) IELM1
C           IF (LNG.EQ.2) WRITE(LU,2001) IELM1
C           IF (LNG.EQ.1) WRITE(LU,3000) IELM2
C           IF (LNG.EQ.2) WRITE(LU,3001) IELM2
C           CALL PLANTE(1)
C           STOP
C         ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSE
C         IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
C         IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
C         IF (LNG.EQ.1) WRITE(LU,2000) IELM1
C         IF (LNG.EQ.2) WRITE(LU,2001) IELM1
C         CALL PLANTE(1)
C         STOP
C       ENDIF
C
C=======================================================================
C     MATRICE  U.N (NON PROGRAMMEE)
C=======================================================================
C
C     ELSEIF(FORMUL(1:15).EQ.'??????         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
C     IF(FORMUL(16:16).EQ.'X') THEN
C       ICOORD=1
C     ELSEIF(FORMUL(16:16).EQ.'Y') THEN
C       ICOORD=2
C     ELSEIF(FORMUL(16:16).EQ.'Z') THEN
C       ICOORD=3
C     ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
C-----------------------------------------------------------------------
C
C       IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
C.......................................................................
C
C         IF(IELM1.EQ.11) THEN
C            CALL MT10AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
C    *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
C    *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
C    *                   XMUL,SF,SG,SH,SU,SV,SW,F,G,H,U,V,W,
C    *                   XEL,YEL,ZEL,SURFAC,
C    *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
C    *                   NELEM,NELMAX,ICOORD)
C
C            TYPDIA='Q'
C            TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
C         ELSE
C           IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
C           IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
C           IF (LNG.EQ.1) WRITE(LU,2000) IELM1
C           IF (LNG.EQ.2) WRITE(LU,2001) IELM1
C           IF (LNG.EQ.1) WRITE(LU,3000) IELM2
C           IF (LNG.EQ.2) WRITE(LU,3001) IELM2
C           CALL PLANTE(1)
C           STOP
C         ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSE
C         IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
C         IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
C         IF (LNG.EQ.1) WRITE(LU,2000) IELM1
C         IF (LNG.EQ.2) WRITE(LU,2001) IELM1
C         CALL PLANTE(1)
C         STOP
C       ENDIF
C
C=======================================================================
C     MATRICE   - PSIJ GRAD(F PSII)
C=======================================================================
C
      ELSEIF(FORMUL(1:15).EQ.'MATGRF         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
      IF(FORMUL(16:16).EQ.'X') THEN
        ICOORD=1
      ELSEIF(FORMUL(16:16).EQ.'Y') THEN
        ICOORD=2
      ELSEIF(FORMUL(16:16).EQ.'Z') THEN
        ICOORD=3
      ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
C-----------------------------------------------------------------------
C
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
C.......................................................................
C
          IF(IELM2.EQ.11) THEN
             CALL MT11AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,SF,F,XEL,YEL,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C
          ELSEIF(IELM2.EQ.12) THEN
             CALL MT11AB
     * (   T(1,1)   ,XM(1,ABQ(1,2,S)),XM(1,ABQ(1,3,S)),XM(1,ABQ(1,4,S)),
     *  XM(1,ABQ(2,1,S)),   T(1,2)   ,XM(1,ABQ(2,3,S)),XM(1,ABQ(2,4,S)),
     *  XM(1,ABQ(3,1,S)),XM(1,ABQ(3,2,S)),   T(1,3)   ,XM(1,ABQ(3,4,S)),
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C
          ELSEIF(IELM2.EQ.13) THEN
             CALL MT11AC(
     *       T(1,1)   ,XM(1,ACQ(1,2,S)),XM(1,ACQ(1,3,S)),
     *       XM(1,ACQ(1,4,S)),XM(1,ACQ(1,5,S)),
     *       XM(1,ACQ(1,6,S)),XM(1,ACQ(2,1,S)),
     *       T(1,2)  ,XM(1,ACQ(2,3,S)),
     *       XM(1,ACQ(2,4,S)),XM(1,ACQ(2,5,S)),
     *       XM(1,ACQ(2,6,S)),XM(1,ACQ(3,1,S)),
     *       XM(1,ACQ(3,2,S)),T(1,3)  ,XM(1,ACQ(3,4,S)),
     *       XM(1,ACQ(3,5,S)),XM(1,ACQ(3,6,S)),
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *          IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT11BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C         ELEMENT DE COLONNE TRIANGLE LINEAIRE
          ELSEIF(IELM2.EQ.11) THEN
             CALL MT11BA
     *         (     T(1,1)     ,XM(1,BAQ(1,2,S)),XM(1,BAQ(1,3,S)),
     *          XM(1,BAQ(2,1,S)),     T(1,2)     ,XM(1,BAQ(2,3,S)),
     *          XM(1,BAQ(3,1,S)),XM(1,BAQ(3,2,S)),     T(1,3)     ,
     *          XM(1,BAQ(4,1,S)),XM(1,BAQ(4,2,S)),XM(1,BAQ(4,3,S)),
     *          XMUL,SF,F,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C   MATRICE  PSIJ GRAD(F)   U .GRAD(PSII)
C=======================================================================
C
      ELSEIF(FORMUL(1:15).EQ.'MATUGH         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
      IF(FORMUL(16:16).EQ.'X') THEN
        ICOORD=1
      ELSEIF(FORMUL(16:16).EQ.'Y') THEN
        ICOORD=2
      ELSEIF(FORMUL(16:16).EQ.'Z') THEN
        ICOORD=3
      ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
C-----------------------------------------------------------------------
C
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
C.......................................................................
C
          IF(IELM2.EQ.11) THEN
             CALL MT12AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,SF,SU,SV,F,U,V,XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C
          ELSEIF(IELM2.EQ.12) THEN
             CALL MT12AB
     * (   T(1,1)   ,XM(1,ABQ(1,2,S)),XM(1,ABQ(1,3,S)),XM(1,ABQ(1,4,S)),
     *  XM(1,ABQ(2,1,S)),   T(1,2)   ,XM(1,ABQ(2,3,S)),XM(1,ABQ(2,4,S)),
     *  XM(1,ABQ(3,1,S)),XM(1,ABQ(3,2,S)),   T(1,3)   ,XM(1,ABQ(3,4,S)),
     *          XMUL,SF,SU,SV,F,U,V,
     *          XEL,YEL,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P2
C.......................................................................
C
          ELSEIF(IELM2.EQ.13) THEN
             CALL MT12AC( 
     *       T(1,1)   ,XM(1,ACQ(1,2,S)),XM(1,ACQ(1,3,S)),
     *       XM(1,ACQ(1,4,S)),XM(1,ACQ(1,5,S)),
     *       XM(1,ACQ(1,6,S)),XM(1,ACQ(2,1,S)),
     *       T(1,2)  ,XM(1,ACQ(2,3,S)),
     *       XM(1,ACQ(2,4,S)),XM(1,ACQ(2,5,S)),
     *       XM(1,ACQ(2,6,S)),XM(1,ACQ(3,1,S)),
     *       XM(1,ACQ(3,2,S)),T(1,3)  ,XM(1,ACQ(3,4,S)),
     *       XM(1,ACQ(3,5,S)),XM(1,ACQ(3,6,S)),
     *          XMUL,SF,SU,SV,F,U,V,
     *          XEL,YEL,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          IKLE(1,5),IKLE(1,6),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
             CALL MT12BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *          XMUL,SF,SU,SV,F,U,V,XEL,YEL,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
          ELSEIF(IELM2.EQ.11) THEN
             CALL MT12BA
     *         (     T(1,1)     ,XM(1,BAQ(1,2,S)),XM(1,BAQ(1,3,S)),
     *          XM(1,BAQ(2,1,S)),     T(1,2)     ,XM(1,BAQ(2,3,S)),
     *          XM(1,BAQ(3,1,S)),XM(1,BAQ(3,2,S)),     T(1,3)     ,
     *          XM(1,BAQ(4,1,S)),XM(1,BAQ(4,2,S)),XM(1,BAQ(4,3,S)),
     *          XMUL,SF,SU,SV,F,U,V,
     *          XEL,YEL,SURFAC,
     *          IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *          NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C=======================================================================
C   MATRICE  PSIJ GRAD(PSII) (SIGNE CHANGE PAR RAPPORT A 3.0)
C   (PSIJ GRAD(PSII) DANS LE CAS B/A)
C=======================================================================
C
      ELSEIF(FORMUL(1:15).EQ.'MATGRA         ') THEN
C
C     LE CARACTERE 16 EST LA COORDONNEE CHOISIE
C
      IF(FORMUL(16:16).EQ.'X') THEN
        ICOORD=1
      ELSEIF(FORMUL(16:16).EQ.'Y') THEN
        ICOORD=2
      ELSEIF(FORMUL(16:16).EQ.'Z') THEN
        ICOORD=3
      ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
C-----------------------------------------------------------------------
C
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
C.......................................................................
C
          IF(IELM2.EQ.11) THEN
             CALL MT13AA(   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
          ELSEIF(IELM2.EQ.12) THEN
             CALL MT13AB(   T(1,1)   ,XM(1,ABQ(1,2,S)),XM(1,ABQ(1,3,S)),
     *                   XM(1,ABQ(1,4,S)),
     *                   XM(1,ABQ(2,1,S)),   T(1,2)   ,XM(1,ABQ(2,3,S)),
     *                   XM(1,ABQ(2,4,S)),
     *                   XM(1,ABQ(3,1,S)),XM(1,ABQ(3,2,S)),   T(1,3)   ,
     *                   XM(1,ABQ(3,4,S)),
     *                   XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
             CALL MT13BA
     *         (     T(1,1)     ,XM(1,BAQ(1,2,S)),XM(1,BAQ(1,3,S)),
     *          XM(1,BAQ(2,1,S)),     T(1,2)     ,XM(1,BAQ(2,3,S)),
     *          XM(1,BAQ(3,1,S)),XM(1,BAQ(3,2,S)),     T(1,3)     ,
     *          XM(1,BAQ(4,1,S)),XM(1,BAQ(4,2,S)),XM(1,BAQ(4,3,S)),
     *          XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          ELSEIF(IELM2.EQ.12) THEN
             CALL MT13BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *          XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P2
        ELSEIF(IELM1.EQ.13) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
             CALL MT13CA
     *         (     T(1,1)     ,XM(1,CAQ(1,2,S)),XM(1,CAQ(1,3,S)),
     *          XM(1,CAQ(2,1,S)),     T(1,2)     ,XM(1,CAQ(2,3,S)),
     *          XM(1,CAQ(3,1,S)),XM(1,CAQ(3,2,S)),     T(1,3)     ,
     *          XM(1,CAQ(4,1,S)),XM(1,CAQ(4,2,S)),XM(1,CAQ(4,3,S)),
     *          XM(1,CAQ(5,1,S)),XM(1,CAQ(5,2,S)),XM(1,CAQ(5,3,S)),
     *          XM(1,CAQ(6,1,S)),XM(1,CAQ(6,2,S)),XM(1,CAQ(6,3,S)),
     *          XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P2
          ELSEIF(IELM2.EQ.13) THEN
             CALL MT13CC
     * (   T(1,1)   ,XM(1,PPQ(1,2,S)),XM(1,PPQ(1,3,S)),
     *   XM(1,PPQ(1,4,S)),XM(1,PPQ(1,5,S)),XM(1,PPQ(1,6,S)),
     *   XM(1,PPQ(2,1,S)),T(1,2),XM(1,PPQ(2,3,S)),XM(1,PPQ(2,4,S)),
     *   XM(1,PPQ(2,5,S)),XM(1,PPQ(2,6,S)),XM(1,PPQ(3,1,S)),
     *   XM(1,PPQ(3,2,S)),T(1,3),XM(1,PPQ(3,4,S)),XM(1,PPQ(3,5,S)),
     *   XM(1,PPQ(3,6,S)),XM(1,PPQ(4,1,S)),XM(1,PPQ(4,2,S)),
     *   XM(1,PPQ(4,3,S)),T(1,4),XM(1,PPQ(4,5,S)),XM(1,PPQ(4,6,S)),
     *   XM(1,PPQ(5,1,S)),XM(1,PPQ(5,2,S)),XM(1,PPQ(5,3,S)),
     *   XM(1,PPQ(5,4,S)),  T(1,5)   ,XM(1,PPQ(5,6,S)),
     *   XM(1,PPQ(6,1,S)),XM(1,PPQ(6,2,S)) ,XM(1,PPQ(6,3,S)),
     *   XM(1,PPQ(6,4,S)),XM(1,PPQ(6,5,S)),T(1,6),
     *          XMUL,XEL,YEL,NELEM,NELMAX,ICOORD)
C
             TYPDIA='Q'
             TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C-----------------------------------------------------------------------
C         ERREUR SUR L'ELEMENT DE COLONNE
C-----------------------------------------------------------------------
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
C>>>>
C=======================================================================
C     MATRICE MURD
C=======================================================================
C
      ELSEIF(FORMUL(1:6).EQ.'MAMURD') THEN
C
C     LE CARACTERE 7 DIT SI SIGMAG OU NON
C
      SIGMAG = .FALSE.
      IF(FORMUL(7:7).EQ.'2') SIGMAG = .TRUE.
C
C     LE CARACTERE 8 DONNE LES DETAILS POUR APPEL DE VC04PP
C
      SPECAD = .FALSE.
      IF(FORMUL(8:8).EQ.'2') SPECAD = .TRUE.
C
C     LES CARACTERES 14 A 16 DONNENT L'OPTION DU SCHEMA
C
      IF(FORMUL(14:16).EQ.'PSI') LEGO = .FALSE.
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE PRISME P1
        IF(IELM1.EQ.41) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE PRISME P1
          IF(IELM2.EQ.41) THEN
             CALL MT14PP(T,XM,PPQ(1,1,S),LEGO,
     *                   XMUL,SU,SV,SW,U,V,W,SF,SG,SH,F,G,H,
     *                   XEL,YEL,ZEL,SURFAC,IKLE,NELEM,NELMAX,SIGMAG,
     *                   SPECAD)
C
            TYPDIA='Q'
            TYPEXT='Q'
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C.......................................................................
C         ERREUR SUR L'ELEMENT DE COLONNE
C.......................................................................
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       AUTRE ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C       ELSEIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C<<<<
C
C=======================================================================
C     MATRICE DE BOUSSINESQ
C=======================================================================
C
      ELSEIF(FORMUL(1:7).EQ.'FFBT   ') THEN
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE P1
        IF(IELM1.EQ.11) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE P1
          IF(IELM2.EQ.11) THEN
            CALL MT99AA
     *                  (   T(1,1)   ,XM(1,AAQ(1,2,S)),XM(1,AAQ(1,3,S)),
     *                   XM(1,AAQ(2,1,S)),   T(1,2)   ,XM(1,AAQ(2,3,S)),
     *                   XM(1,AAQ(3,1,S)),XM(1,AAQ(3,2,S)),   T(1,3)   ,
     *                   XMUL,SF,F,XEL,YEL,SURFAC,
     *                   IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                   NELEM,NELMAX,FORMUL,TDIA,TEXT)
C
            TYPDIA = TDIA
            TYPEXT = TEXT
C.......................................................................
C         AUTRE
C.......................................................................
C
C         ELSEIF
C
C.......................................................................
C         ERREUR SUR L'ELEMENT DE COLONNE
C.......................................................................
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C-----------------------------------------------------------------------
C       AUTRE ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C       ELEMENT DE LIGNE TRIANGLE QUASI-BULLE
C
        ELSEIF(IELM1.EQ.12) THEN
C
C.......................................................................
C         ELEMENT DE COLONNE TRIANGLE QUASI-BULLE
          IF(IELM2.EQ.12) THEN
            CALL MT99BB
     * (   T(1,1)   ,XM(1,BBQ(1,2,S)),XM(1,BBQ(1,3,S)),XM(1,BBQ(1,4,S)),
     *  XM(1,BBQ(2,1,S)),   T(1,2)   ,XM(1,BBQ(2,3,S)),XM(1,BBQ(2,4,S)),
     *  XM(1,BBQ(3,1,S)),XM(1,BBQ(3,2,S)),   T(1,3)   ,XM(1,BBQ(3,4,S)),
     *  XM(1,BBQ(4,1,S)),XM(1,BBQ(4,2,S)),XM(1,BBQ(4,3,S)),   T(1,4)   ,
     *         XMUL,SF,F,XEL,YEL,SURFAC,
     *         IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *         NELEM,NELMAX,FORMUL,TDIA,TEXT)
C
C
            TYPDIA = TDIA
            TYPEXT = TEXT
C
C.......................................................................
C         AUTRE
C.......................................................................
C
C
C.......................................................................
C         ERREUR SUR L'ELEMENT DE COLONNE
C.......................................................................
C
          ELSE
            IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
            IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
            IF (LNG.EQ.1) WRITE(LU,2000) IELM1
            IF (LNG.EQ.2) WRITE(LU,2001) IELM1
            IF (LNG.EQ.1) WRITE(LU,3000) IELM2
            IF (LNG.EQ.2) WRITE(LU,3001) IELM2
            CALL PLANTE(1)
            STOP
          ENDIF
C
C-----------------------------------------------------------------------
C       ERREUR SUR L'ELEMENT DE LIGNE
C-----------------------------------------------------------------------
C
C
        ELSE
          IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
          IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
          IF (LNG.EQ.1) WRITE(LU,2000) IELM1
          IF (LNG.EQ.2) WRITE(LU,2001) IELM1
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSE
C
C=======================================================================
C     ERREUR : TYPE DE MATRICE NON PROGRAMME
C=======================================================================
C
        IF (LNG.EQ.1) WRITE(LU,1000) FORMUL
        IF (LNG.EQ.2) WRITE(LU,1001) FORMUL
        CALL PLANTE(1)
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
1000  FORMAT(1X,'MATRIY (BIEF) : MATRICE NON PREVUE : ',A16)
1001  FORMAT(1X,'MATRIY (BIEF) : MATRIX NOT IMPLEMENTED:',A16)
2000  FORMAT(1X,'                POUR IELM1 = ',1I6)
2001  FORMAT(1X,'                FOR IELM1 = ',1I6)
3000  FORMAT(1X,'                ET IELM2 = ',1I6)
3001  FORMAT(1X,'                AND IELM2 = ',1I6)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
