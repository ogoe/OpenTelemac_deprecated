C                       ***************************
                        SUBROUTINE CALC_NEIGHB_GRAD
C                       ***************************
C
     * ( GRAD2 , X , Y , NPOIN ,ZF , TPUNKT )
C
C***********************************************************************
C SISYPHE VERSION 6.0    04/02/2011               O. GOETHEL    
C
C
C***********************************************************************
C
C  FONCTION : CALCULATES THE GRADIENT BETWEEN NEIGHBOURPOINTS
C
C
C***********************************************************************
C
C
      USE DECLARATIONS_SISYPHE, ONLY : NEIGHB_DIFF, ELEM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
c
      INTEGER NPOIN,I,K,TPUNKT(NPOIN)
      DOUBLE PRECISION X(NPOIN),Y(NPOIN),ZF(NPOIN),GRAD2(NPOIN)
      DOUBLE PRECISION DEV(20),KLEIN(NPOIN),GROSS(NPOIN)
C
c
c
      DO I=1,NPOIN
c
c
c
        GRAD2(I)=0.D0
c
        DO K=1,ELEM(I)%ANZ
         NEIGHB_DIFF(I)%HZ(K)=ZF(I)-ZF(ELEM(I)%NR(K))
         IF(NEIGHB_DIFF(I)%HZ(K).LT.0.D0)NEIGHB_DIFF(I)%HZ(K)=0.D0 
c
         NEIGHB_DIFF(I)%L(K)=DSQRT( (X(ELEM(I)%NR(K))-X(I))**2.D0 +
     *                       (Y(ELEM(I)%NR(K))-Y(I))**2.D0 )
C
c
C
         IF(DABS(GRAD2(I)).LT.(NEIGHB_DIFF(I)%HZ(K)/
     *                         NEIGHB_DIFF(I)%L(K)))THEN
            GRAD2(I)=NEIGHB_DIFF(I)%HZ(K)/NEIGHB_DIFF(I)%L(K)
            TPUNKT(I)=ELEM(I)%NR(K)
         ENDIF
c
        ENDDO
C
c
      ENDDO
C
C
      RETURN
      END
C