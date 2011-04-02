C                         ***********************
                          SUBROUTINE POINT_HEIGHT
C                         ***********************
C
     *( ZFNEU , VOL , SPOINT , NPOIN , X , Y , ZF , IKLE ,
     *  NELEM )
C
C***********************************************************************
C SISYPHE VERSION 6.0    04/02/2011               O. GOETHEL    
C
C
C***********************************************************************
C
C  FONCTION : CALCULATES THE NEW HEIGHT DEPENDING ON THE VOLUME
C
C
C***********************************************************************
C
C
C
C
      USE DECLARATIONS_SISYPHE, ONLY: NEIGHB_ELEM, ELEM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
C
C
C
      INTEGER SPOINT,NPOIN,I,NELEM,IKLE(NELEM,3)
      DOUBLE PRECISION X(NPOIN),Y(NPOIN),ZF(NPOIN),VOL,ZF2(NPOIN)
      DOUBLE PRECISION ZFNEU,ZAEHLER,NENNER,KLAMMER
      DOUBLE PRECISION a,b,c,F,s
C
C
      DO I=1,NPOIN
         ZF2(I)=ZF(I)+1.D5
      ENDDO
C
      ZAEHLER = 3.D0 * VOL
      NENNER = 0.D0
C
C
      DO I=1, NEIGHB_ELEM(SPOINT)%ANZ
C
C
      a=DSQRT((X(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1))-
     *         X(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2)))**2.D0 + 
     *        (Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1))-
     *         Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2)))**2.D0)
C
C
      b=DSQRT((X(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1))-
     *         X(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3)))**2.D0 + 
     *        (Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1))-
     *         Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3)))**2.D0)
C 
C
      c=DSQRT((X(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2))-
     *         X(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3)))**2.D0 + 
     *         (Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2))-
     *         Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3)))**2.D0)
C
C
      s=(a+b+c)/2.D0
C
      F=DSQRT(s*(s-a)*(s-b)*(s-c))
C
C
      NENNER = NENNER + F
C
C
C
      IF(SPOINT.EQ.IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1))THEN
C
       KLAMMER = ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2)) + 
     *           ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3))
C
      ELSE IF(SPOINT.EQ.IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2))THEN
C
       KLAMMER = ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1)) + 
     *           ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3))
C
      ELSE IF(SPOINT.EQ.IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3))THEN
C
       KLAMMER = ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1)) + 
     *           ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2))
      ENDIF
C
C
C
      ZAEHLER = ZAEHLER - F * KLAMMER
C
C
      ENDDO
C
C
C     NEUE HOEHE
C
      ZFNEU = ZAEHLER / NENNER
C
      ZFNEU = ZFNEU - 1.D5
C
C
      RETURN
      END
C