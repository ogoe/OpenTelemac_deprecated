C                       *****************
                        SUBROUTINE MATVCT
C                       *****************
C
     *(OP, X , DA,TYPDIA,XA,TYPEXT, Y ,
     * C,IKLE,NPT,NELEM,NELMAX,W,LEGO,IELM1,IELM2,IELMX,LV,
     * S,P,IKLEM1,DIMIKM,LIMVOI,MXPTVS,NPMAX,NPOIN,NPTFR,
     * GLOSEG,SIZGLO,SIZXA)
C
C***********************************************************************
C BIEF VERSION 6.0        05/02/2010  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C 31/03/2008 : ALGIANE FROEHLY (MATMECA) : TRIANGLE QUADRATIQUE
C 05/02/2010 : JMH NSEG1 AND NSEG2 MODIFIED BEFORE CALLING MVSEG
C***********************************************************************
C
C   FONCTION : OPERATIONS MATRICE VECTEUR
C
C   OP EST UNE CHAINE DE 8 CARACTERES QUI INDIQUE L'OPERATION QUI SERA
C   EFFECTUEE SUR LES VECTEURS X,Y ET LA MATRICE M. LE RESULTAT
C   EST LE VECTEUR X (DONT UNE PARTIE NON ASSEMBLEE PEUT ETRE DANS
C   LE TABLEAU W SI LEGO = .FALSE.)
C
C   CES OPERATIONS SONT DIFFERENTES SUIVANT LE TYPE DE DIAGONALE
C   ET LE TYPE DES TERMES EXTRADIAGONAUX.
C
C   OPERATIONS PROGRAMMEES :
C
C      OP = 'X=AY    '  : X = AY
C      OP = 'X=X+AY  '  : X = X + AY
C      OP = 'X=X-AY  '  : X = X - AY
C      OP = 'X=X+CAY '  : X = X + C AY
C      OP = 'X=TAY   '  : X = TA Y (TRANSPOSEE DE A)
C      OP = 'X=X+TAY '  : X = X + TA Y
C      OP = 'X=X-TAY '  : X = X - TA Y
C      OP = 'X=X+CTAY'  : X = X + C TA Y
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      OP        | -->|OPERATION A EFFECTUER
C |      X         |<-- |VECTEUR IMAGE
C |      DA        | -->| DIAGONALE DE LA MATRICE.
C |      TYPDIA    | -->| TYPE DE LA DIAGONALE (CHAINE DE CARACTERES)
C |                |    | TYPDIA = 'Q' : DIAGONALE QUELCONQUE
C |                |    | TYPDIA = 'I' : DIAGONALE IDENTITE.
C |                |    | TYPDIA = '0' : DIAGONALE NULLE.
C |      XA        | -->| TERMES EXTRA-DIAGONAUX DE LA MATRICE
C |      TYPEXT    | -->| TYPE DES TERMES EXTRADIAGONAUX
C |                |    | TYPEXT = 'Q' : QUELCONQUES.
C |                |    | TYPEXT = 'S' : SYMETRIQUES.
C |                |    | TYPEXT = '0' : NULS.
C |      Y         | -->| VECTEUR OPERANDE
C |      C         | -->| CONSTANTE DONNEE
C |      IKLE      | -->| CORRESPONDANCE NUMEROTATIONS LOCALE ET GLOBALE
C |      NPT       | -->| DIMENSION DE LA DIAGONALE
C |      NELEM     | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |      NELMAX    | -->| NOMBRE D'ELEMENTS DU MAILLAGE
C |      W         |<-- | TABLEAUX DE TRAVAIL QUI CONTIENNENT UNE PARTIE
C |                |    | DU RESULTAT SI L'ON CHOISIT LE MODE NON
C |                |    | ASSEMBLE ( LEGO = .FALSE. )
C |      LEGO      | -->| = .TRUE. W1,2,... SONT ASSEMBLES SUR X
C |                |    | =.FALSE. W1,2,... SONT LAISSES TELS QUELS.
C |      IELM1     | -->| TYPE D'ELEMENT (LIGNES DE LA MATRICE).
C |      IELM2     | -->| TYPE D'ELEMENT (COLONNES DE LA MATRICE).
C |      IELMX     | -->| TYPE D'ELEMENT DU RESULTAT
C |                |    | EGAL A IELM1 OU IELM2 SUIVANT L'OPERATION.
C |      LV        | -->| LONGUEUR DU VECTEUR POUR LA VECTORISATION.
C |      S         | -->| TYPE DE STOCKAGE.
C |      P         | -->| TYPE DE PRODUIT MATRICE X VECTEUR.
C |      IKLEM1    | -->| TABLEAU DE CONNECTIVITE UTILISE POUR LE
C |                | -->| STOCKAGE DE TYPE 2.
C |      DIMIKM    | -->| PREMIERE DIMENSION DE IKLEM1.
C |      LIMVOI    | -->| TABLEAU UTILISE POUR LE STOCKAGE 2.
C |      MXPTVS    | -->| PREMIERE DIMENSION DE LIMVOI.
C |                |    | NOMBRE MAXIMUM DE VOISINS D'UN POINT
C |      NPMAX     | -->| NOMBRE MAXIMUM DE POINTS DU MAILLAGE.
C |      NPOIN     | -->| NOMBRE DE POINTS DU MAILLAGE.
C |      NPTFR     | -->| NOMBRE DE POINTS DE BORD.
C |      NELBOR    | -->| NUMEROS DES ELEMENTS DE BORD.
C |      NULONE    | -->| NUMEROS LOCAUX DANS L'ELEMENT DES POINTS DE
C |                |    | BORD.
C |      SIZXA     | -->| FIRST DIMENSION OF ARRAY XA
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   : MV0303 , MV0404 , MV0606 , PLANTE
C                        MV0304 , MV0305 , MV0403 , MV0202 
C                        ASSVEC
C
C***********************************************************************
C
      USE BIEF, EX_MATVCT => MATVCT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: IELM1,IELM2,IELMX,NPOIN,NPMAX,S,P,SIZXA
      INTEGER, INTENT(INOUT) :: NPT
      INTEGER, INTENT(IN) :: NELEM,NELMAX,LV,DIMIKM,MXPTVS,NPTFR,SIZGLO
      INTEGER, INTENT(IN) :: IKLE(NELMAX,*),IKLEM1(*),LIMVOI(*)
      INTEGER, INTENT(IN) :: GLOSEG(SIZGLO,2)
      CHARACTER(LEN=8), INTENT(IN)    :: OP
      CHARACTER(LEN=1),INTENT(IN)     :: TYPDIA,TYPEXT
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      DOUBLE PRECISION, INTENT(IN)    :: Y(*),DA(*),XA(SIZXA,*),C
      DOUBLE PRECISION, INTENT(INOUT) :: W(NELMAX,*)
      LOGICAL, INTENT(IN)             :: LEGO      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER NSEG1,NSEG2,SYM,NPT2      
C
      INTEGER AAQ(3,3,2),ABQ(3,4,2),BAQ(4,3,2)
      INTEGER ACQ(3,6,2),BBQ(4,4,2),CAQ(6,3,2),PPQ(6,6,2)
      INTEGER AAS(3,3,2),BBS(4,4,2),PPS(6,6,2)
      INTEGER OOS(2,2,2)
C
      DOUBLE PRECISION Z(1)    
C
      INTRINSIC MIN
C
C     CES DATA FIGURENT AUSSI DANS MATVCT
C
C     DATA OOS/  0 ,  1 ,
C    *           1 ,  0 ,
C S=2 NON PREVU
C    *           0 ,  0 ,
C    *           0 ,  0 /
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
C     P1 P2 EBE NON SYMETRIQUE (S=1)
      DATA ACQ/   0 ,  6 ,  11,
     *            1 ,  0 ,  12,
     *            2 ,  7 ,  0 ,
     *            3 ,  8 ,  13,
     *         	  4 ,  9 ,  14, 	
     *		  5 ,  10,  15, 
C S=2 NON PREVU  
     *            0 ,  0 ,  0 ,
     *            0 ,  0 ,  0 ,
     *            0 ,  0 ,  0 ,
     *            0 ,  0 ,  0 ,
     *         	  0 ,  0 ,  0 , 	
     *		  0 ,  0 ,  0 /
C     QUASI-BULLE P1 EBE NON SYMETRIQUE (S=1)
      DATA BAQ/  0 ,  3 ,  5 ,  7 ,
     *           1 ,  0 ,  6 ,  8 ,
     *           2 ,  4 ,  0 ,  9 ,
C     QUASI-BULLE P1 EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
     *           0 ,  7 ,  3 ,  4 ,
     *           1 ,  0 ,  8 ,  5 ,
     *           9 ,  2 ,  0 ,  6 /

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
C     PRISMES P1-P1 ET TRIANGLES P2 EBE SYMETRIQUE (S=1)
      DATA PPS/  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,
     *           1 ,  0 ,  6 ,  7 ,  8 ,  9 ,
     *           2 ,  6 ,  0 , 10 , 11 , 12 ,
     *           3 ,  7 , 10 ,  0 , 13 , 14 ,
     *           4 ,  8 , 11 , 13 ,  0 , 15 ,
     *           5 ,  9 , 12 , 14 , 15 ,  0 ,
C     PRISMES P1-P1 ET TRIANGLES P2 EBE PRE-ASSEMBLE SYMETRIQUE (S=2) 
C     - NON PREVU
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C     PRISMES P1-P1 ET TRIANGLES P2 EBE NON SYMETRIQUE (S=1)
      DATA PPQ/  0 , 16 , 17 , 18 , 19 , 20 ,
     *           1 ,  0 , 21 , 22 , 23 , 24 ,
     *           2 ,  6 ,  0 , 25 , 26 , 27 ,
     *           3 ,  7 , 10 ,  0 , 28 , 29 ,
     *           4 ,  8 , 11 , 13 ,  0 , 30 ,
     *           5 ,  9 , 12 , 14 , 15 ,  0 ,
C    PRISMES P1-P1 ET TRIANGLES P2 EBE PRE-ASSEMBLE NON SYMETRIQUE (S=2)
C     - NON PREVU
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
     *           0 ,  0 ,  0 ,  0 ,  0 ,  0 /
C
C-----------------------------------------------------------------------
C
      IF(S.EQ.1) THEN
C
C-----------------------------------------------------------------------
C
C     STOCKAGE CLASSIQUE EBE ET PRODUIT MATRICE X VECTEUR EBE CLASSIQUE
C
C-----------------------------------------------------------------------
C
      IF(IELM1.EQ.1) THEN
C
        IF(IELM2.EQ.1) THEN
          CALL MV0202(OP, X , DA,TYPDIA,
     *                XA(1,1),XA(1,2),TYPEXT, Y,C,
     *                IKLE(1,1),IKLE(1,2),
     *                NPT,NELEM,W(1,1),W(1,2))
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSEIF(IELM1.EQ.11) THEN
C
        IF(IELM2.EQ.11) THEN
          IF(TYPEXT(1:1).EQ.'S') THEN
            CALL MV0303(OP, X , DA,TYPDIA,
     *                  XA(1,AAS(1,2,S)),
     *                  XA(1,AAS(1,3,S)),
     *                  XA(1,AAS(2,1,S)),
     *                  XA(1,AAS(2,3,S)),
     *                  XA(1,AAS(3,1,S)),
     *                  XA(1,AAS(3,2,S)),
     *                  TYPEXT,Y,C,
     *                  IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                  NPT,NELEM,
     *                  W(1,1),W(1,2),W(1,3))
          ELSE
            CALL MV0303(OP, X , DA,TYPDIA,
     *                  XA(1,AAQ(1,2,S)),
     *                  XA(1,AAQ(1,3,S)),
     *                  XA(1,AAQ(2,1,S)),
     *                  XA(1,AAQ(2,3,S)),
     *                  XA(1,AAQ(3,1,S)),
     *                  XA(1,AAQ(3,2,S)),
     *                  TYPEXT,Y,C,
     *                  IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                  NPT,NELEM,
     *                  W(1,1),W(1,2),W(1,3))
          ENDIF
C
        ELSEIF(IELM2.EQ.12) THEN
C
          CALL MV0304(OP, X , DA,TYPDIA,
     *                XA(1,ABQ(1,2,S)),
     *                XA(1,ABQ(1,3,S)),
     *                XA(1,ABQ(1,4,S)),
     *                XA(1,ABQ(2,1,S)),
     *                XA(1,ABQ(2,3,S)),
     *                XA(1,ABQ(2,4,S)),
     *                XA(1,ABQ(3,1,S)),
     *                XA(1,ABQ(3,2,S)),
     *                XA(1,ABQ(3,4,S)),
     *                TYPEXT, Y,C,
     *                IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                NPT,NELEM,
     *                W(1,1),W(1,2),W(1,3),W(1,4))
C 
         ELSEIF(IELM2.EQ.13) THEN
C
          NPT2=NBPTS(IELM2)
          CALL MV0306(OP, X , DA,TYPDIA,
     *                XA(1,ACQ(1,2,S)), XA(1,ACQ(1,3,S)),              
     *                XA(1,ACQ(1,4,S)), XA(1,ACQ(1,5,S)),
     *                XA(1,ACQ(1,6,S)), XA(1,ACQ(2,1,S)),
     *                XA(1,ACQ(2,3,S)), XA(1,ACQ(2,4,S)),
     *                XA(1,ACQ(2,5,S)), XA(1,ACQ(2,6,S)),
     *                XA(1,ACQ(3,1,S)), XA(1,ACQ(3,2,S)),
     *                XA(1,ACQ(3,4,S)), XA(1,ACQ(3,5,S)),
     *                XA(1,ACQ(3,6,S)), 
     *                TYPEXT, Y,C,
     *                IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *                NPT,NPT2,NELEM,
     *                W(1,1),W(1,2),W(1,3),
     *                W(1,4),W(1,5),W(1,6))
C 
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSEIF(IELM1.EQ.12.OR.IELM2.EQ.31.OR.IELM2.EQ.51) THEN
C
        IF(IELM2.EQ.12.OR.IELM2.EQ.31.OR.IELM2.EQ.51) THEN
          IF(TYPEXT(1:1).EQ.'S') THEN
            CALL MV0404(OP, X , DA,TYPDIA,
     *                  XA(1,BBS(1,2,S)),
     *                  XA(1,BBS(1,3,S)),
     *                  XA(1,BBS(1,4,S)),
     *                  XA(1,BBS(2,1,S)),
     *                  XA(1,BBS(2,3,S)),
     *                  XA(1,BBS(2,4,S)),
     *                  XA(1,BBS(3,1,S)),
     *                  XA(1,BBS(3,2,S)),
     *                  XA(1,BBS(3,4,S)),
     *                  XA(1,BBS(4,1,S)),
     *                  XA(1,BBS(4,2,S)),
     *                  XA(1,BBS(4,3,S)),
     *                  TYPEXT, Y,C,
     *                  IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                  NPT,NELEM,
     *                  W(1,1),W(1,2),W(1,3),W(1,4))
          ELSE
            CALL MV0404(OP, X , DA,TYPDIA,
     *                  XA(1,BBQ(1,2,S)),
     *                  XA(1,BBQ(1,3,S)),
     *                  XA(1,BBQ(1,4,S)),
     *                  XA(1,BBQ(2,1,S)),
     *                  XA(1,BBQ(2,3,S)),
     *                  XA(1,BBQ(2,4,S)),
     *                  XA(1,BBQ(3,1,S)),
     *                  XA(1,BBQ(3,2,S)),
     *                  XA(1,BBQ(3,4,S)),
     *                  XA(1,BBQ(4,1,S)),
     *                  XA(1,BBQ(4,2,S)),
     *                  XA(1,BBQ(4,3,S)),
     *                  TYPEXT, Y,C,
     *                  IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                  NPT,NELEM,
     *                  W(1,1),W(1,2),W(1,3),W(1,4))
          ENDIF
        ELSEIF(IELM2.EQ.11) THEN
          CALL MV0403(OP, X , DA,TYPDIA,
     *                XA(1,BAQ(1,2,S)),
     *                XA(1,BAQ(1,3,S)),
     *                XA(1,BAQ(2,1,S)),
     *                XA(1,BAQ(2,3,S)),
     *                XA(1,BAQ(3,1,S)),
     *                XA(1,BAQ(3,2,S)),
     *                XA(1,BAQ(4,1,S)),
     *                XA(1,BAQ(4,2,S)),
     *                XA(1,BAQ(4,3,S)),
     *                TYPEXT, Y,C,
     *                IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
     *                NPT,NELEM,
     *                W(1,1),W(1,2),W(1,3),W(1,4))
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSEIF(IELM1.EQ.41.OR.IELM1.EQ.13) THEN
C
        IF(IELM2.EQ.41.OR.IELM2.EQ.13) THEN
C	   
          CALL MV0606(OP, X , DA,TYPDIA,XA,TYPEXT, Y,C,
     *                IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *                NPT,NELEM,NELMAX,
     *                W(1,1),W(1,2),W(1,3),W(1,4),W(1,5),W(1,6))
C	
        ELSEIF(IELM2.EQ.11) THEN
C
C         HERE IELM1=13
          NPT2=NBPTS(IELM1)
          CALL MV0603(OP, X , DA,TYPDIA,
     *                XA(1,CAQ(1,2,S)),XA(1,CAQ(1,3,S)),
     *                XA(1,CAQ(2,1,S)),XA(1,CAQ(2,3,S)),
     *                XA(1,CAQ(3,1,S)),XA(1,CAQ(3,2,S)),
     *                XA(1,CAQ(4,1,S)),XA(1,CAQ(4,2,S)),
     *                XA(1,CAQ(4,3,S)),XA(1,CAQ(5,1,S)),
     *                XA(1,CAQ(5,2,S)),XA(1,CAQ(5,3,S)),
     *                XA(1,CAQ(6,1,S)),XA(1,CAQ(6,2,S)),
     *                XA(1,CAQ(6,3,S)),
     *                TYPEXT, Y,C,
     *                IKLE(1,1),IKLE(1,2),IKLE(1,3),
     *                IKLE(1,4),IKLE(1,5),IKLE(1,6),
     *                NPT,NPT2,NELEM,
     *                W(1,1),W(1,2),W(1,3),
     *                W(1,4),W(1,5),W(1,6))
C        
        ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
C  IELM1 NON PREVU : ERREUR
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
        IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
100     FORMAT(1X,'MATVCT (BIEF) : ELEMENTS ',1I2,' ET ',1I2,/,1X,
     *            'ET STOCKAGE ',1I2,'   CAS NON PREVU')
101     FORMAT(1X,'MATVCT (BIEF) : ELEMENTS ',1I2,' AND ',1I2,/,1X,
     *            'AND STORAGE ',1I2,'   CASE NOT IMPLEMENTED')
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C     ASSEMBLAGE FINAL EVENTUEL DE X
C
C     ICI COMME INIT = FALSE IL N'Y A PEUT-ETRE PAS BESOIN DE NPT
      NPT = NBPTS(IELMX)
      IF(LEGO) CALL ASSVEC(X,IKLE,NPT,NELEM,NELMAX,IELMX,W,
     *                                 .FALSE.,LV,.FALSE.,Z)
C
      ELSEIF(S.EQ.3.AND.P.EQ.2) THEN
C
C-----------------------------------------------------------------------
C
C  STOCKAGE SEGMENT ET PRODUIT MATRICE X VECTEUR FRONTAL
C
C-----------------------------------------------------------------------
C
      IF(IELM1.EQ.1) THEN
C
        IF(IELM2.EQ.1) THEN
C         CALL MW0202(OP, X , DA,TYPDIA,
C    *                XA(1,1),XA(1,2),TYPEXT, Y,C,
C    *                IKLE(1,1),IKLE(1,2),
C    *                NPT,NELEM,
C    *                W(1,1),W(1,2))
C       ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSEIF(IELM1.EQ.11) THEN
C
        IF(IELM2.EQ.11) THEN
C
          CALL MW0303(OP, X , DA,TYPDIA,XA,TYPEXT, Y,C,
     *                IKLEM1,DIMIKM,LIMVOI,MXPTVS,NPMAX,NPOIN,W)
C
        ELSEIF(IELM2.EQ.12) THEN
C
C         CALL MW0304(OP, X , DA,TYPDIA,
C    *                XA(1,1),XA(1,2),XA(1,3),
C    *                XA(1,4),XA(1,5),XA(1,6),
C    *                XA(1,7),XA(1,8),XA(1,9),
C    *                TYPEXT, Y,C,
C    *                IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
C    *                NPT,NELEM,
C    *                W(1,1),W(1,2),W(1,3),W(1,4))
C       ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSEIF(IELM1.EQ.12) THEN
C
        IF(IELM2.EQ.12) THEN
C           CALL MW0404(OP, X , DA,TYPDIA,
C    *                  XA(1,1),XA(1,2),XA(1,3),XA(1,1),
C    *                  XA(1,4),XA(1,5),XA(1,2),XA(1,4),
C    *                  XA(1,6),XA(1,3),XA(1,5),XA(1,6),
C    *                  TYPEXT, Y,C,
C    *                  IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
C    *                  NPT,NELEM,
C    *                  W(1,1),W(1,2),W(1,3),W(1,4))
C       ELSEIF(IELM2.EQ.11) THEN
C         CALL MW0403(OP, X , DA,TYPDIA,
C    *                XA(1,1),XA(1,2),XA(1,3),
C    *                XA(1,4),XA(1,5),XA(1,6),
C    *                XA(1,7),XA(1,8),XA(1,9),
C    *                TYPEXT, Y,C,
C    *                IKLE(1,1),IKLE(1,2),IKLE(1,3),IKLE(1,4),
C    *                NPT,NELEM,
C    *                W(1,1),W(1,2),W(1,3),W(1,4))
C       ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
      ELSEIF(IELM1.EQ.41) THEN
C
        IF(IELM2.EQ.41) THEN
C
C         CALL MW0606(OP, X , DA,TYPDIA,XA,TYPEXT, Y,C,
C    *                IKLE(1,1),IKLE(1,2),IKLE(1,3),
C    *                IKLE(1,4),IKLE(1,5),IKLE(1,6),
C    *                NPT,NELEM,NELMAX,
C    *                W(1,1),W(1,2),W(1,3),W(1,4),W(1,5),W(1,6))
C       ELSE
          IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
          IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
          CALL PLANTE(1)
          STOP
        ENDIF
C
C  IELM1 NON PREVU : ERREUR
C
      ELSE
C
        IF (LNG.EQ.1) WRITE(LU,100) IELM1,IELM2,S
        IF (LNG.EQ.2) WRITE(LU,101) IELM1,IELM2,S
        CALL PLANTE(1)
        STOP
C
      ENDIF
C
C  STOCKAGE PAR SEGMENTS
C
      ELSEIF(S.EQ.3.AND.P.EQ.1) THEN
C
C-----------------------------------------------------------------------
C
C  STOCKAGE SEGMENT ET PRODUIT MATRICE X VECTEUR CLASSIQUE
C
C-----------------------------------------------------------------------
C
      NSEG1 = NBSEG(IELM1)
      NSEG2 = NBSEG(IELM2)
C
C     IN LINEAR-QUADRATIC RECTANGULAR MATRICES, PURELY QUADRATIC
C     SEGMENTS ARE NOT CONSIDERED (NUMBER 13,14 AND 15, SO 3 PER ELEMENT)
C
      IF(IELM1.EQ.11.AND.IELM2.EQ.13) THEN
        NSEG2=NSEG2-3*NELEM
      ELSEIF(IELM1.EQ.13.AND.IELM2.EQ.11) THEN
        NSEG1=NSEG1-3*NELEM
      ENDIF
C
      IF(TYPEXT(1:1).EQ.'Q') THEN
        SYM = MIN(NSEG1,NSEG2)
      ELSE
        SYM = 0
      ENDIF
      CALL MVSEG (OP, X , DA,TYPDIA,XA(1,1),XA(SYM+1,1),
     *            TYPEXT,Y,C,NPT,NELEM,NSEG1,NSEG2,
     *            GLOSEG(1,1),GLOSEG(1,2),IELM1,IELM2)
C
C-----------------------------------------------------------------------
C
C  STOCKAGE NON PREVU
C
C-----------------------------------------------------------------------
C
      ELSE
        IF (LNG.EQ.1) WRITE(LU,102) S,P
        IF (LNG.EQ.2) WRITE(LU,103) S,P
102     FORMAT(1X,'MATVCT (BIEF) : ',1I2,' ET ',1I2,/,1X,
     *            'STOCKAGE ET PRODUIT MATRICE-VECTEUR INCOMPATIBLES')
103     FORMAT(1X,'MATVCT (BIEF) : ',1I2,' AND ',1I2,/,1X,
     *            'STORAGE AND MATRIX-VECTOR PRODUCT INCOMPATIBLE')
        CALL PLANTE(1)
        STOP
      ENDIF
C
C=======================================================================
C
      RETURN
      END
