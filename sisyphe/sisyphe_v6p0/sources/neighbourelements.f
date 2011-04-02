C                         ****************************
                          SUBROUTINE NEIGHBOURELEMENTS
C                         ****************************
C
     *( NPOIN , NELEM , IKLE)
C
C***********************************************************************
C SISYPHE VERSION 6.0    04/02/2011               O. GOETHEL    
C
C
C***********************************************************************
C
C  FONCTION : CALCULATES THE SEDIMENT VOLUME BENEATH A MESHPOINT
C
C
C***********************************************************************
C
C
C
C
      USE DECLARATIONS_SISYPHE, ONLY : NEIGHB_ELEM, NEIGHB_DIFF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
C
C
       INTEGER I,NPOIN,NELEM,K,L,IKLE(NELEM,3) 
C
C
c
      ALLOCATE(NEIGHB_ELEM(NPOIN))
      ALLOCATE(NEIGHB_DIFF(NPOIN))
C
C
       FORALL (I=1:NPOIN) NEIGHB_ELEM(I)%ANZ=0
       FORALL (I=1:NPOIN,K=1:20) NEIGHB_ELEM(I)%NR(K)=0
c
C
      DO I=1,NPOIN
C
C
       DO K=1,NELEM
C
C
        DO L=1,3
C
         IF(I.EQ.IKLE(K,L)) THEN
C
          NEIGHB_ELEM(I)%ANZ=NEIGHB_ELEM(I)%ANZ+1
          NEIGHB_ELEM(I)%NR(NEIGHB_ELEM(I)%ANZ)=K
C
         ENDIF
C
        ENDDO
C
C
       ENDDO
C
C
      ENDDO
C
C
C
      RETURN
      END
C