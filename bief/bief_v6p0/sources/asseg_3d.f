C                       *******************
                        SUBROUTINE ASSEG_3D
C                       *******************
C
     *(FLOW,F,NPOIN3,NPLAN,NSEG2D,GLOSEG,SIZGLO,INIFLO)
C
C***********************************************************************
C BIEF VERSION 6.0        18/05/09    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION  : ASSEMBLING HORIZONTAL EDGE BY EDGE FLUXES ON POINTS
C                         ----------
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    FLOW        | -->| FLUXES (SIZE OF FLOW MAY NOT EXCEED 
C |                |    | NSEG2D*NPLAN, THOUGH THE TOTAL NUMBER OF
C |                |    | SEGMENTS IS LARGER)
C |    F           |<---| F WHERE THE FLUXES ARE ASSEMBLED
C |    NPOIN3      | -->| NUMBER OF POINTS
C |    NPLAN       | -->| NUMBER OF PLANES
C |    NSEG2D      | -->| NUMBER OF SEGMENTS IN 2D
C |    GLOSEG      | -->| GLOBAL NUMBER OF THE 2 POINTS OF A SEGMENT
C |    SIZGLO      | -->| FIRST DIMENSION OF GLOSEG
C |    INIFLO      | -->| IF(YES) F WILL BE INITIALISED AT 0.
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : CVDFTR
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER, INTENT(IN)             :: NSEG2D,NPOIN3,SIZGLO,NPLAN
      INTEGER, INTENT(IN)             :: GLOSEG(SIZGLO,2)
      DOUBLE PRECISION, INTENT(IN)    :: FLOW(*)
C                                        HERE * = NESG2D*NPLAN < NSEG3D
      DOUBLE PRECISION, INTENT(INOUT) :: F(NPOIN3)
      LOGICAL, INTENT(IN)             :: INIFLO
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER ISEG,I
C
C-----------------------------------------------------------------------
C
C     INITIALISATION OF FLOW TO 0.D0
C
      IF(INIFLO) THEN
        DO I = 1,NPOIN3
          F(I) = 0.D0
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C
C     ASSEMBLING THE FLUXES OF HORIZONTAL SEGMENTS
C
      DO ISEG = 1,NSEG2D*NPLAN
        F(GLOSEG(ISEG,1))=F(GLOSEG(ISEG,1))+FLOW(ISEG)
        F(GLOSEG(ISEG,2))=F(GLOSEG(ISEG,2))-FLOW(ISEG)  
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
