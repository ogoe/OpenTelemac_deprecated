C                       *****************
                        SUBROUTINE ALLMAT
C                       *****************
C
     *( MAT , NOM , IELM1 , IELM2 , CFG , TYPDIA , TYPEXT )
C
C***********************************************************************
C BIEF VERSION 5.1            01/03/95    J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C  FONCTION  : ALLOCATION EN MEMOIRE D'UNE STRUCTURE DE MATRICE REELLE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   MAT          | -->| MATRICE A ALLOUER
C |   NOM          | -->| NOM FORTRAN DU TABLEAU
C |   IELM1        | -->| TYPE D'ELEMENT PAR LIGNES
C |   IELM2        | -->| TYPE D'ELEMENT PAR COLONNES
C |   CFG(1)       | -->| TYPE DE STOCKAGE.
C |   CFG(2)       | -->| CHOIX DU PRODUIT MATRICE X VECTEUR
C |   TYPDIA       | -->| TYPE DE DIAGONALE ('Q', 'I' OU '0')
C |   TYPEXT       | -->| TYPE DE TERMES EXTRA-DIAGONAUX:'Q','S' OU '0'
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_ALLMAT => ALLMAT
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_OBJ)  , INTENT(INOUT) :: MAT
      CHARACTER(LEN=6), INTENT(IN)    :: NOM
      INTEGER         , INTENT(IN)    :: IELM1,IELM2,CFG(2)
      CHARACTER(LEN=1), INTENT(IN)    :: TYPDIA,TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELMD
C
      CHARACTER(LEN=6) :: NAME
C
C-----------------------------------------------------------------------
C  HEADER COMMON TO ALL OBJECTS
C-----------------------------------------------------------------------
C
C     TO CHECK MEMORY CRASHES
C
      MAT%KEY = 123456
C
C     TYPE OF OBJECT (HERE MATRIX)
C
      MAT%TYPE = 3
C
C     NAME OF OBJECT
C
      MAT%NAME = NOM
C
C-----------------------------------------------------------------------
C  SPECIFIC PART FOR MATRICES
C-----------------------------------------------------------------------
C
C     ELEMENT OF DIAGONAL (SMALLEST ELEMENT IF MATRIX IS RECTANGULAR)
      IELMD = IELM1
      IF(NBPTS(IELM2).LT.NBPTS(IELM1)) IELMD = IELM2
C
C  TYPE OF STORAGE
C
      MAT%STO = CFG(1)
C
C  TYPES OF ELEMENTS FOR LINE AND COLUMN
C
      MAT%ELMLIN = IELM1
      MAT%ELMCOL = IELM2
C
C  ALLOCATION OF THE DIAGONAL (TILL THERE MAT%D IS ONLY A POINTER)
C
C     MAT%D WILL POINT TO AN EXISTING BIEF_OBJ
      ALLOCATE(MAT%D)
C
      NAME = 'D' // NOM(1:5)
      IF(TYPDIA(1:1).EQ.'Q') THEN
C        ONLY CASE WHERE THE DIAGONAL DOES EXIST
         CALL ALLVEC(1,MAT%D,NAME,IELMD,1,2)
      ELSE
         CALL ALLVEC(1,MAT%D,NAME,0    ,1,0)
      ENDIF
C     TYPE IS FORGOTTEN UNTIL INITIALISATION OF MATRIX
C     MAT%TYPDIA = TYPDIA(1:1)
      MAT%TYPDIA = '?'
C
C     ALLOCATION OF OFF-DIAGONAL TERMS (AS FOR DIAGONAL)
C
      ALLOCATE(MAT%X)
C
      NAME = 'X' // NOM(1:5)
C
         CALL ALLVEC(1,MAT%X,NAME,
     *               DIM1_EXT(IELM1,IELM2,CFG(1),TYPEXT),
     *               DIM2_EXT(IELM1,IELM2,CFG(1),TYPEXT),0)
C
C     TYPE IS FORGOTTEN UNTIL INITIALISATION OF MATRIX
C     MAT%TYPEXT = TYPEXT(1:1)
      MAT%TYPEXT = '?'
C
C     MATRIX x VECTOR PRODUCT
      MAT%PRO = CFG(2)
C
C-----------------------------------------------------------------------
C
C     IF(LNG.EQ.1) WRITE(LU,*) 'MATRICE : ',NOM,' ALLOUEE'
C     IF(LNG.EQ.2) WRITE(LU,*) 'MATRIX: ',NOM,' ALLOCATED'
C
C-----------------------------------------------------------------------
C
      RETURN
      END
