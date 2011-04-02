C                         ***********************
                          SUBROUTINE POINT_VOLUME
C                         ***********************
C
     *( VOL , SPOINT , NPOIN , X , Y , ZF , IKLE , NELEM)
C
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
      USE DECLARATIONS_SISYPHE, ONLY: NEIGHB_ELEM, ELEM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C
C
C
      INTEGER SPOINT,NPOIN,I,NELEM,IKLE(NELEM,3)
      DOUBLE PRECISiON X(NPOIN),Y(NPOIN),ZF(NPOIN),VOL,ZF2(NPOIN)
      DOUBLE PRECISION a,b,c,F,s
C
C
C
      DO I=1,NPOIN
         ZF2(I)=ZF(I)+1.D5 !shift to positive value
      ENDDO
C
C
C
      VOL=0.D0
      DO I=1,NEIGHB_ELEM(SPOINT)%ANZ
C
C     sites of a triangle
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
     *        (Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2))-
     *         Y(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3)))**2.D0)
C
C     area of a triangle
C
      s=(a+b+c)/2.D0
C
      F=DSQRT(s*(s-a)*(s-b)*(s-c))
C
C     prism volume
C
      VOL=VOL + ( ( ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),1)) + 
     *              ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),2)) +
     *              ZF2(IKLE(NEIGHB_ELEM(SPOINT)%NR(I),3))  ) * F/3.D0)
C
C
C
C
      ENDDO
C
C
      RETURN
      END
C