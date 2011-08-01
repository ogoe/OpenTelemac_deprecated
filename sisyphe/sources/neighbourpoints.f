c                 **************************
                  SUBROUTINE NEIGHBOURPOINTS
C                 **************************
C
     *      ( IKLE , NELEM , NPOIN )
C
C***********************************************************************
C SISYPHE VERSION 6.0    04/02/2011               O. GOETHEL    
C
C
C***********************************************************************
C
C  FONCTION : FIND THE NEIGHBOURS OF A MESHPOINT
C
C
C***********************************************************************
C
C
      USE DECLARATIONS_SISYPHE, ONLY: ELEM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      LOGICAL TEST1,TEST2,TEST3
      INTEGER NPOIN,NELEM,I,K,IELEM,IKLE(NELEM,3)
C
C
C
      ALLOCATE(ELEM(NPOIN))
C
C
C
      FORALL (I=1:NPOIN) ELEM(I)%ANZ=0
      FORALL (I=1:NPOIN,K=1:20) ELEM(I)%NR(K)=0
      FORALL (I=1:NPOIN,K=1:20) ELEM(I)%DIFFX(K)=0.D0
      FORALL (I=1:NPOIN,K=1:20) ELEM(I)%DIFFY(K)=0.D0
      FORALL (I=1:NPOIN,K=1:20) ELEM(I)%DIFFZF(K)=0.D0
C
C
C
C
      DO I=1,NPOIN
       DO IELEM=1,NELEM
C
        IF(I.EQ.IKLE(IELEM,1)) THEN
C
         IF(ELEM(I)%ANZ.NE.0) THEN
          TEST2=.TRUE.
          TEST3=.TRUE.
          DO K=1,ELEM(I)%ANZ
           IF(ELEM(I)%NR(K).EQ.IKLE(IELEM,2)) TEST2=.FALSE.
          ENDDO
          DO K=1,ELEM(I)%ANZ
           IF(ELEM(I)%NR(K).EQ.IKLE(IELEM,3)) TEST3=.FALSE.
          ENDDO
C
          IF(TEST2) THEN
           ELEM(I)%ANZ=ELEM(I)%ANZ+1
           ELEM(I)%NR(ELEM(I)%ANZ)=IKLE(IELEM,2)
          ENDIF
          IF(TEST3) THEN
           ELEM(I)%ANZ=ELEM(I)%ANZ+1
           ELEM(I)%NR(ELEM(I)%ANZ)=IKLE(IELEM,3)
          ENDIF
C
         ELSE
C
         ELEM(I)%NR(1)=IKLE(IELEM,2)
         ELEM(I)%NR(2)=IKLE(IELEM,3)
         ELEM(I)%ANZ=2
C
        ENDIF
        ENDIF
C
C
        IF(I.EQ.IKLE(IELEM,2)) THEN
C
         IF(ELEM(I)%ANZ.NE.0) THEN
          TEST1=.TRUE.
          TEST3=.TRUE.
          DO K=1,ELEM(I)%ANZ
           IF(ELEM(I)%NR(K).EQ.IKLE(IELEM,1)) TEST1=.FALSE.
          ENDDO
          DO K=1,ELEM(I)%ANZ
           IF(ELEM(I)%NR(K).EQ.IKLE(IELEM,3)) TEST3=.FALSE.
          ENDDO
C
          IF(TEST1) THEN
           ELEM(I)%ANZ=ELEM(I)%ANZ+1
           ELEM(I)%NR(ELEM(I)%ANZ)=IKLE(IELEM,1)
          ENDIF
          IF(TEST3) THEN
           ELEM(I)%ANZ=ELEM(I)%ANZ+1
           ELEM(I)%NR(ELEM(I)%ANZ)=IKLE(IELEM,3)
          ENDIF
C
         ELSE
C
         ELEM(I)%NR(1)=IKLE(IELEM,1)
         ELEM(I)%NR(2)=IKLE(IELEM,3)
         ELEM(I)%ANZ=2
C
        ENDIF
        ENDIF
C
C
        IF(I.EQ.IKLE(IELEM,3)) THEN
C
         IF(ELEM(I)%ANZ.NE.0) THEN
          TEST1=.TRUE.
          TEST2=.TRUE.
          DO K=1,ELEM(I)%ANZ
           IF(ELEM(I)%NR(K).EQ.IKLE(IELEM,1)) TEST1=.FALSE.
          ENDDO
          DO K=1,ELEM(I)%ANZ
           IF(ELEM(I)%NR(K).EQ.IKLE(IELEM,2)) TEST2=.FALSE.
          ENDDO
C
          IF(TEST1) THEN
           ELEM(I)%ANZ=ELEM(I)%ANZ+1
           ELEM(I)%NR(ELEM(I)%ANZ)=IKLE(IELEM,1)
          ENDIF
          IF(TEST2) THEN
           ELEM(I)%ANZ=ELEM(I)%ANZ+1
           ELEM(I)%NR(ELEM(I)%ANZ)=IKLE(IELEM,2)
          ENDIF
C
         ELSE
C
         ELEM(I)%NR(1)=IKLE(IELEM,1)
         ELEM(I)%NR(2)=IKLE(IELEM,2)
         ELEM(I)%ANZ=2
C
        ENDIF
        ENDIF
C
C
       ENDDO
      ENDDO
C
      RETURN
      END
C