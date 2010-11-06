C                       ********************
                        SUBROUTINE FLUX3DLIM
C                       ********************
C
     *(FLOW,FLULIM,NPLAN,NSEG2D)
C
C***********************************************************************
C BIEF VERSION 6.0        19/05/09    J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C  FONCTION  : LIMITING 3D HORIZONTAL EDGE BY EDGE FLUXES ON POINTS
C                          ----------
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|_______________________________________________
C |    FLOW        | -->| FLUXES (SIZE OF FLOW MAY NOT EXCEED 
C |                |    | NSEG2D*NPLAN, THOUGH THE TOTAL NUMBER OF
C |                |    | SEGMENTS IS LARGER)
C |    FLULIM      | -->| LIMITING FACTOR OF 2D SEGMENTS
C |    NPLAN       | -->| NUMBER OF PLANES
C |    NSEG2D      | -->| NUMBER OF SEGMENTS IN 2D
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
      INTEGER, INTENT(IN)             :: NSEG2D,NPLAN
      DOUBLE PRECISION, INTENT(INOUT) :: FLOW(*)
C                                        HERE * = NESG2D*NPLAN < NSEG3D
      DOUBLE PRECISION, INTENT(IN)    :: FLULIM(NSEG2D)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER ISEG,ISEG3D,IPLAN
C
C-----------------------------------------------------------------------
C
C     LIMITATION OF 3D FLUXES BY COEFFICIENT OF 2D FLUXES
C
      DO IPLAN=1,NPLAN
        DO ISEG=1,NSEG2D
          ISEG3D=ISEG+(IPLAN-1)*NSEG2D
          FLOW(ISEG3D)=FLOW(ISEG3D)*FLULIM(ISEG)
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN
      END
