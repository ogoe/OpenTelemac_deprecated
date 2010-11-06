C                    *************************************
                     DOUBLE PRECISION FUNCTION STA_DIS_CUR
C                    *************************************
C
     *(IFRLIQ,FLUX,PTS,QZ,NFRLIQ,ZN)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.9   17/08/94  J-M HERVOUET (LNHE) 01 30 87 80 18
C
C***********************************************************************
C
C FONCTION  : DONNE LA VALEUR DE LA COTE DE LA SURFACE LIBRE
C             EN FONCTION DU DEBIT EN INTERPOLANT DANS UNE
C             COURBE DE TARAGE
C
C-----------------------------------------------------------------------
C
C FUNCTION  : GIVES THE PRESCRIBED VALUE OF FREE SURFACE AT
C             A LIQUID BOUNDARY BY INTERPOLATING IN A STAGE-DISCHARGE
C             CURVE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   IFRLIQ       | -->| LIQUID BOUNDARY NUMBER
C |   FLUX         | -->| ACTUAL FLUX AT THIS BOUNDARY
C |   PTS          | -->| NUMBER OF POINTS IN THE STAGE-DISCHARGE CURVE
C |   QZ           | -->| ARRAY WITH STAGE-DISCHARGE CURVES
C |   NFRLIQ       | -->| NUMBER OF LIQUID BOUNDARIES
C |   ZN           | -->| PREVIOUS ELEVATION (FOR RELAXATION)
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C  APPELE PAR : BORD
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: IFRLIQ,NFRLIQ,PTS
      DOUBLE PRECISION, INTENT(IN) :: ZN,FLUX,QZ(2,NFRLIQ,PTS)
C                                                         PTS AT LEAST
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER I
      DOUBLE PRECISION GOAL,TETA,Q1,Q2
C
C-----------------------------------------------------------------------
C
      Q1=QZ(1,IFRLIQ,1)
      Q2=QZ(1,IFRLIQ,PTS)
      IF(FLUX.LE.Q1) THEN
C       OUTSIDE THE CURVE WITH LOWER DISCHARGE
        GOAL=QZ(2,IFRLIQ,1)
      ELSEIF(FLUX.GE.Q2) THEN
C       OUTSIDE THE CURVE WITH HIGHER DISCHARGE
        GOAL=QZ(2,IFRLIQ,PTS)
      ELSE
C       IN BETWEEN: CASE WITH INTERPOLATION
        I=1
1       CONTINUE
        Q2=QZ(1,IFRLIQ,I+1)
        IF(FLUX.GE.Q1.AND.FLUX.LE.Q2) THEN
          TETA=(FLUX-Q1)/MAX(Q2-Q1,1.D-8)
          GOAL=QZ(2,IFRLIQ,I)+TETA*(QZ(2,IFRLIQ,I+1)-QZ(2,IFRLIQ,I))
        ELSE
          I=I+1
          Q1=Q2
          GO TO 1
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C
C     RELAXATION OF RESULT
C
      STA_DIS_CUR=ZN+0.02D0*(GOAL-ZN)                
C
C-----------------------------------------------------------------------
C
      RETURN
      END
