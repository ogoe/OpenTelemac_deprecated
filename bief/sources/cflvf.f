C                       ****************
                        SUBROUTINE CFLVF
C                       ****************
C
     *(DTMAX,HSTART,H,FXMAT,FXMATPAR,MAS,DT,FXBOR,SMH,YASMH,TAB1,NSEG,
     * NPOIN,NPTFR,GLOSEG,SIZGLO,MESH,MSK,MASKPT)
C
C***********************************************************************
C BIEF VERSION 6.0           30/11/2009   C-T PHAM (LNHE) 01 30 87 85 93
C
C JMH LE 11/04/2008 : YASMH AJOUTE
C JMH LE 10/06/2008 : SIZGLO AJOUTE
C JMH LE 02/10/2008 : PARALLELISME (FXMATPAR AJOUTE, ETC.)
C JMH LE 30/11/2009 : REFINED COMPUTATION OF DTMAX (AS IN 3D)
C***********************************************************************
C
C  FONCTION  : COMPUTING THE MAXIMUM TIME STEP ENABLING MONOTONICITY
C              IN THE ADVECTION STEP              
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    DTMAX       |<-- | MAXIMUM TIME STEP FOR STABILITY
C |    HSTART      | -->| H AT BEGINNING OF SUB TIME STEP
C |    H           | -->| H AT THE END OF FULL TIME STEP
C |    FXMAT       | -->| FLUXES
C |    FXMATPAR    | -->| FLUXES ASSEMBLED IN PARALLEL
C |    MAS         | -->| INTEGRAL OF TEST FUNCTIONS (=AREA AROUND POINTS)
C |    DT          | -->| TIME STEP
C |    FXBOR       | -->| BOUNDARY FLUXES
C |    SMH         | -->| RIGHT HAND SIDE OF CONTINUITY EQUATION
C |    YASMH       | -->| IF YES, TAKE SHM INTO ACCOUNT
C |    TAB1        | -->| WORK ARRAY
C |    NSEG        | -->| NUMBER OF SEGMENTS
C |    NPOIN       | -->| NUMBER OF POINTS IN THE MESH
C |    NPTFR       | -->| NUMBER OF BOUNDARY POINTS
C |    GLOSEG      | -->| GLOBAL NUMBER OF THE 2 POINTS OF A SEGMENT
C |    SIZGLO      | -->| FIRST DIMENSION OF GLOSEG
C |    MESH        | -->| MESH STRUCTURE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_CFLVF => CFLVF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NSEG,NPOIN,NPTFR,SIZGLO
      INTEGER, INTENT(IN)             :: GLOSEG(SIZGLO,2)
      DOUBLE PRECISION, INTENT(INOUT) :: DTMAX
      DOUBLE PRECISION, INTENT(IN)    :: DT,HSTART(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: H(NPOIN),MAS(NPOIN),SMH(NPOIN)
C                                              NOT NPTFR, SEE TRACVF                                            
      DOUBLE PRECISION, INTENT(IN)    :: FXBOR(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FXMAT(NSEG),FXMATPAR(NSEG)
      LOGICAL, INTENT(IN)             :: YASMH,MSK
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: TAB1
      TYPE(BIEF_MESH), INTENT(INOUT)  :: MESH
      TYPE(BIEF_OBJ), INTENT(IN)      :: MASKPT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
      DOUBLE PRECISION DENOM,A,B
C
C-----------------------------------------------------------------------
C
C COMPUTATION OF CRITERION FOR COURANT NUMBER
C
      DO I = 1,NPOIN
        TAB1%R(I) = 0.D0
      ENDDO
C
C WE USE HERE FXMAT ASSEMBLED IN PARALLEL FOR UPWINDING
C
      DO I = 1,NSEG
        IF(FXMATPAR(I).LT.0.D0) THEN
          TAB1%R(GLOSEG(I,1)) = TAB1%R(GLOSEG(I,1)) + FXMAT(I)
        ELSEIF(FXMATPAR(I).GT.0.D0) THEN
          TAB1%R(GLOSEG(I,2)) = TAB1%R(GLOSEG(I,2)) - FXMAT(I)             
        ENDIF  
      ENDDO
C
      IF(NCSIZE.GT.1) CALL PARCOM(TAB1,2,MESH)
C
C     MASKING TAB1
C
      IF(MSK) THEN
        CALL OS('X=XY    ',X=TAB1,Y=MASKPT)
      ENDIF
C
C STABILITY (AND MONOTONICITY) CRITERION
C
C NOTE THAT TAB1(I)<0 MIN(FXBOR(I),0.D0)<0 AND -MAX(SMH(I),0.D0)<0
C           SO ABS(TAB1(I)+MIN(FXBOR(I),0.D0)-MAX(SMH(I),0.D0))=
c              -(TAB1(I)+MIN(FXBOR(I),0.D0)-MAX(SMH(I),0.D0))      
C
C     ANY TIME LARGER THAN THE REMAINING DT
      DTMAX = 2.D0*DT
C
C     SEE RELEASE NOTES 5.7, CRITERION AT THE END OF 4.4 PAGE 33
C     BUT HERE THE FINAL H IS NOT H(N+1) BUT A FUNCTION OF DTMAX ITSELF
C     H FINAL = HSTART + DTMAX/DT *(H-HSTART)
C
      IF(YASMH) THEN
        DO I = 1,NPOIN
          DENOM=TAB1%R(I)+MIN(FXBOR(I),0.D0)-MAX(SMH(I),0.D0)
          A=-MAS(I)/MIN(DENOM,-1.D-12)
          B=DT+A*(HSTART(I)-H(I))
          IF(B.GT.0.D0) THEN
            DTMAX = MIN(DTMAX,A*HSTART(I)*DT/B)
          ENDIF
        ENDDO
      ELSE
        DO I = 1,NPOIN
          DENOM=TAB1%R(I)+MIN(FXBOR(I),0.D0)
          A=-MAS(I)/MIN(DENOM,-1.D-12)
          B=DT+A*(HSTART(I)-H(I))
          IF(B.GT.0.D0) THEN
            DTMAX = MIN(DTMAX,A*HSTART(I)*DT/B)
          ENDIF
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
