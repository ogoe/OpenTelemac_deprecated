C                       *****************
                        SUBROUTINE FDNRST 
C                       *****************
C
     *(IFRM,ITO,X,Y,NODENRS,NPOIN2,IFRM1,ITOP1)
C
C***********************************************************************
C  TELEMAC 2D VERSION 5.7    03/04/07    LEO POSTMA (DELFT HYDRAULICS)
C
C***********************************************************************
C
C     FONCTION  : FINDS THE NEAREST FROM -1 AND TO +1 POINTER
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   IFRM         | -->|  
C |   ITO          | -->|  
C |   X,Y          | -->|  NODE COORDINATES
C |   NODENRS      | -->|  IF > 0 : NODE NUMBER
C |                |    |  IF < 0 : - OPEN BOUNDARY NODE NUMBER
C |   NPOIN2       | -->|  NUMBER OF POINTS IN 2D
C |   NPTFR        | -->|  NOMBRE DE POINTS FRONTIERES 
C |--------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)          :: IFRM,ITO,NPOIN2
      INTEGER, INTENT(IN)          :: NODENRS(NPOIN2)
      INTEGER, INTENT(INOUT)       :: IFRM1,ITOP1
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN2), Y(NPOIN2)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPOIN
      DOUBLE PRECISION XFRM1,XTOP1,YFRM1,YTOP1,DISFRM,DISTO,DX,DY
C
C-----------------------------------------------------------------------
C
      XFRM1  = X(IFRM)
      YFRM1  = Y(IFRM)
      XTOP1  = X(ITO )
      YTOP1  = Y(ITO )
      DISFRM = XFRM1-XTOP1
      XFRM1  = XFRM1 + DISFRM
      XTOP1  = XTOP1 - DISFRM
      DISFRM = YFRM1-YTOP1
      YFRM1  = YFRM1 + DISFRM
      YTOP1  = YTOP1 - DISFRM
C
      DX     = XFRM1-X(1)
      DY     = YFRM1-Y(1)
      DISFRM = DX*DX + DY*DY
      DX     = XTOP1-X(1)
      DY     = YTOP1-Y(1)
      DISTO  = DX*DX + DY*DY
      IFRM1  = 1
      ITOP1  = 1
      DO IPOIN = 2, NPOIN2
         DX     = XFRM1-X(IPOIN)
         DY     = YFRM1-Y(IPOIN)
         DX     = DX*DX + DY*DY
         IF(DX.LT.DISFRM) THEN
           DISFRM = DX
           IFRM1  = IPOIN
         ENDIF
         DX     = XTOP1-X(IPOIN)
         DY     = YTOP1-Y(IPOIN)
         DX     = DX*DX + DY*DY
         IF(DX.LT.DISTO) THEN
           DISTO  = DX
           ITOP1  = IPOIN
         ENDIF
      ENDDO
      IFRM1 = NODENRS(IFRM1)
      ITOP1 = NODENRS(ITOP1)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
