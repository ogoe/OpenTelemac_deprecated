C                    ********************
                     SUBROUTINE SANDSLIDE
C                    ********************
C
     *( ZF     , ZF1  )
C
C***********************************************************************
C SISYPHE VERSION 6.0    04/02/2011               O. GOETHEL    
C
C
C***********************************************************************
C
C  FONCTION : CALCULATES A SANDSLIDE BASED ON A CALCULATED FRICTION ANGLE
C
C
C***********************************************************************
C
C
C
C
      USE BIEF
      USE DECLARATIONS_SISYPHE, ONLY: MESH, NPOIN
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
      INTEGER I,ISCHL,I2,I3,K,DEBUG
      LOGICAL BETA1,TEST
C
      INTEGER TPUNKT(NPOIN)
C
C
      DOUBLE PRECISION  BETA(NPOIN),BETA_CRIT(NPOIN)
      DOUBLE PRECISION  VOL,VOL2,DIFVOL,VOL3
      DOUBLE PRECISION  GRADR(NPOIN)
C
      DOUBLE PRECISION  ISTHDIF,SOLLHDIF,DIF, ZFNEU2
C
C
      TYPE(BIEF_OBJ) :: ZF1,ZF
C
C
      DEBUG=1
C
      IF(DEBUG>0) PRINT *,'SANDSLIDE: FIRST CALCULATION OF GRADIENTS'
C
       CALL CALC_NEIGHB_GRAD(GRADR,MESH%X%R,MESH%Y%R ,NPOIN ,ZF1%R,
     *           TPUNKT)
C
C
C 
       IF(DEBUG>0) PRINT *,'SANDSLIDE: CALCULATING BOTTOM ANGLE'
C    
      DO I=1,NPOIN
       BETA(I)=DATAN(DABS(GRADR(I)))
      ENDDO
C

      CALL CALC_CRIT_BOTTOM_ANGLE(BETA_CRIT)
C
      IF(DEBUG>0) PRINT *,
     * 'SANDSLIDE: FIRST COMPARISON OF ANGLE AND CRITICAL ANGLE'
C     Check BETA_CRIT
      K=0
      BETA1=.FALSE.
C
      DO I=1,NPOIN
C
        IF (BETA(I).GT.BETA_CRIT(I)) THEN
          BETA1=.TRUE.
C        GOTO 10
          K=K+1
        END IF
c     
      END DO
C
10    CONTINUE
C
      ISCHL=1
C
      IF (BETA1) THEN
         print *,'SANDSLIDE: SLIDE AT ',K,' POINTS INITIATED'
      ENDIF
C
C
C    
      DO WHILE (BETA1)
      IF(DEBUG>1) print *,'ENTERED WHILE LOOP'
C
C     
      DO I=1,NPOIN
C
C
       IF(BETA(I).GT.BETA_CRIT(I))THEN
C
C
C
       ISTHDIF = ZF1%R(I) - ZF1%R(TPUNKT(I))
c
       SOLLHDIF = DTAN(BETA_CRIT(I) - 0.035D0) * !BETA_CRIT(I) - 2. degrees -> threshold
     *            DABS(MESH%X%R(I) - MESH%X%R(TPUNKT(I)))
C
c
       IF(ISTHDIF.LE.SOLLHDIF)THEN
c        print*,'DIF negativ'
        GOTO 110
c       print*,'test' 
       ENDIF
c
c
       DIF = ISTHDIF - SOLLHDIF
c       IF(ISCHL.LT.12)THEN
c       DIF = DIF /(20.D0 - DBLE(ISCHL))
c       ELSE
       DIF=DIF/4.D0
c       ENDIF
C
C
        CALL POINT_VOLUME(VOL,I,NPOIN,MESH%X%R,MESH%Y%R,ZF1%R,
     *               MESH%IKLE%I,MESH%NELEM)
C
C
        ZF1%R(I) = ZF1%R(I) - DIF
C
C
        CALL POINT_VOLUME(VOL2,I,NPOIN,MESH%X%R,MESH%Y%R,ZF1%R,
     *               MESH%IKLE%I,MESH%NELEM)
C
C
        DIFVOL = VOL - VOL2
C
        CALL POINT_VOLUME(VOL3,TPUNKT(I),NPOIN,MESH%X%R,MESH%Y%R,
     *               ZF1%R,MESH%IKLE%I,MESH%NELEM)
C        
C
        CALL POINT_HEIGHT(ZFNEU2 , VOL3+DIFVOL , TPUNKT(I) , 
     *             NPOIN , MESH%X%R , MESH%Y%R , ZF1%R ,
     *             MESH%IKLE%I , MESH%NELEM )
C
C
      ZF1%R(TPUNKT(I))=ZFNEU2   
C
c        CALL GRAD(GRADX%R,GRADY%R,NKR,MESH%X%R,MESH%Y%R , NPOIN ,ZF1%R
c     *            ,TPUNKTX , TPUNKTY )
C
c      DO K=1,NPOIN
c       BETAX(K)=DATAN(DABS(GRADX%R(K)))
c       BETAY(K)=DATAN(DABS(GRADY%R(K)))
c      END DO
C
        ENDIF
C
C 
 110    CONTINUE
C
        ENDDO
C
C
        CALL CALC_NEIGHB_GRAD(GRADR,MESH%X%R,MESH%Y%R , NPOIN ,ZF1%R
     *            ,TPUNKT )
C
      DO K=1,NPOIN
       BETA(K)=DATAN(DABS(GRADR(K)))
      END DO    
C  
C
C
      BETA1=.FALSE.
      DO I=1,NPOIN
       IF (BETA(I).GT.BETA_CRIT(I)) THEN
          BETA1=.TRUE.
          GOTO 120
       END IF
      END DO
C
120    CONTINUE
C
c      DO I=1,NPOIN
c         print*,'punkt',i,'zf',zf1%r(I)
c      enddo
C      PRINT *,'SANDSLIDE ITERATION NR.',ISCHL
c      IF (ISCHL.EQ.30) THEN
c         BETA1=.FALSE.
c      END IF
      ISCHL=ISCHL+1
      END DO  
Cend do while
      IF(ISCHL.GT.1) print*,'SANDSLIDE: NUMBER OF LOOPS:',ISCHL-1
      RETURN
      END
C