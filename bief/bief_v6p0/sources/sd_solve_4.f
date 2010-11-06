C		       ********************* 
                       SUBROUTINE SD_SOLVE_4
C		       *********************
C 
     &(NPOIN,NSEGB,GLOSEGB,DAB1,DAB2,DAB3,DAB4,XAB1,XAB2,XAB3,XAB4,
     & XX1,XX2,CVB1,CVB2,INFOGR,TYPEXT)
C
C***********************************************************************
C BIEF VERSION 5.9     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION: RESOLUTION DIRECTE D'UN SYSTEME BLOCK 2 X 2 PAR
C            PERMUTATION MINIMUM DEGREE ET DECOMPOSITION LDLT 
C
C            TRANSFORMATION STOCKAGE SEGMENT EN STOCKAGE COMPACT (MORSE)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPOIN        | -->| NOMBRE D'INCONNUES
C |   NSEGB        | -->| NOMBRE DE SEGMENTS 
C |   GLOSEG       | -->| NUMEROS GLOBAUX DES POINTS DES SEGMENTS
C |   DA,XA        | -->| DIAGONALES ET TERMES EXTRA-DIAGONAUX DES
C |                |    | MATRICES
C |   XX1,XX2      |<-- | SOLUTIONS
C |   CVB1,CVB2    | -->| SECONDS MEMBRES
C |   INFOGR       | -->| IF, YES INFORMATIONS ON LISTING
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_SOLVE_4 => SD_SOLVE_4
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,NSEGB
      INTEGER, INTENT(IN) :: GLOSEGB(NSEGB*2)
      LOGICAL, INTENT(IN) :: INFOGR
      DOUBLE PRECISION, INTENT(IN)    :: DAB1(NPOIN),DAB2(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: DAB3(NPOIN),DAB4(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: XAB1(NSEGB),XAB2(NSEGB)
      DOUBLE PRECISION, INTENT(IN)    :: XAB3(NSEGB),XAB4(NSEGB)
      DOUBLE PRECISION, INTENT(INOUT) :: XX1(NPOIN),XX2(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: CVB1(NPOIN),CVB2(NPOIN)    
      CHARACTER(LEN=1), INTENT(IN)    :: TYPEXT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C              
      INTEGER NPBLK,NSEGBLK,I
C
      INTEGER, ALLOCATABLE          :: GLOSEG4(:)
      DOUBLE PRECISION, ALLOCATABLE :: XA(:),DA(:)
      DOUBLE PRECISION, ALLOCATABLE :: RHS(:),XINC(:)     
C
      INTEGER SIZE_GLOSEG4,SIZE_DA,SIZE_XA,SIZE_RHS,SIZE_XINC
C
      DATA SIZE_GLOSEG4/0/
      DATA SIZE_DA     /0/
      DATA SIZE_XA     /0/
      DATA SIZE_RHS    /0/
      DATA SIZE_XINC   /0/
C
      SAVE
C
C-----------------------------------------------------------------------
C
      NPBLK=NPOIN*2
      NSEGBLK=4*NSEGB+NPOIN 
C 
      IF(SIZE_GLOSEG4.EQ.0) THEN      
        ALLOCATE(GLOSEG4(2*NSEGBLK))
        SIZE_GLOSEG4=    2*NSEGBLK
      ELSEIF(            2*NSEGBLK.GT.SIZE_GLOSEG4) THEN
        DEALLOCATE(GLOSEG4)
        ALLOCATE(GLOSEG4(2*NSEGBLK))
        SIZE_GLOSEG4=    2*NSEGBLK
      ENDIF 
      IF(SIZE_DA.EQ.0) THEN
        ALLOCATE(DA(NPBLK))
        SIZE_DA=    NPBLK
      ELSEIF(       NPBLK.GT.SIZE_DA) THEN
        DEALLOCATE(DA)
        ALLOCATE(DA(NPBLK))
        SIZE_DA=    NPBLK
      ENDIF                
      IF(SIZE_XA.EQ.0) THEN      
        ALLOCATE(XA(2*NSEGBLK))
        SIZE_XA=    2*NSEGBLK
      ELSEIF(       2*NSEGBLK.GT.SIZE_XA) THEN
        DEALLOCATE(XA)
        ALLOCATE(XA(2*NSEGBLK))
        SIZE_XA=    2*NSEGBLK
      ENDIF              
      IF(SIZE_RHS.EQ.0) THEN
        ALLOCATE(RHS(NPBLK))
        SIZE_RHS=    NPBLK
      ELSEIF(        NPBLK.GT.SIZE_RHS) THEN
        DEALLOCATE(RHS)
        ALLOCATE(RHS(NPBLK))
        SIZE_RHS=    NPBLK
      ENDIF 
      IF(SIZE_XINC.EQ.0) THEN
        ALLOCATE(XINC(NPBLK))
        SIZE_XINC=    NPBLK
      ELSEIF(         NPBLK.GT.SIZE_XINC) THEN
        DEALLOCATE(XINC)
        ALLOCATE(XINC(NPBLK))
        SIZE_XINC=    NPBLK
      ENDIF             
C
C-----------------------------------------------------------------------
C       
C     1. SECOND MEMBRE DU SYSTEME
C     ===========================
C
      DO I=1,NPOIN
        RHS(I)      = CVB1(I)
        RHS(I+NPOIN)= CVB2(I)
      ENDDO       
C              
C     2. CONSTRUCTION STOCKAGE SEGMENT MATRICE BLOCK (DE 4)
C     ===================================================== 
C	    
      CALL SD_STRSG4(NPOIN,NSEGB,GLOSEGB,NPBLK,NSEGBLK,GLOSEG4)
C     
      CALL SD_FABSG4(NPOIN,NSEGB,DAB1,DAB2,DAB3,DAB4,
     *               XAB1,XAB2,XAB3,XAB4,NPBLK,NSEGBLK,DA,XA)
C
C     3. RESOLUTION COMME UNE MATRICE SYMETRIQUE NORMALE
C     ================================================== 
C 
      CALL SD_SOLVE_1(NPBLK,NSEGBLK,GLOSEG4,NSEGBLK,DA,XA,
     *                XINC,RHS,INFOGR,TYPEXT)
C
C     4. RECUPERATION DES INCONNUES
C     ============================= 
C 
      DO I=1,NPOIN
        XX1(I)= XINC(I)
        XX2(I)= XINC(I+NPOIN)
      ENDDO
C
C-----------------------------------------------------------------------
C 
      RETURN              
      END
