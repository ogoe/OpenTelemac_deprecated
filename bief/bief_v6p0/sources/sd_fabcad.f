C                       ********************
                        SUBROUTINE SD_FABCAD
C                       ********************
C
     *(NPBLK,NSEGBLK,IN,IP,ISEGIP,
     * INDTRI,ISTRI,INX,IPX,ACTRI,XA1,XA2,DA,AC)    
C
C***********************************************************************
C BIEF VERSION 5.7     20/11/06   E. RAZAFINDRAKOTO (LNH) 01 30 87 74 03
C
C***********************************************************************
C
C  IMPORTANT NOTICE: INSPIRED FROM PACKAGE CMLIB3 - YALE UNIVERSITE-YSMP
C
C  FONCTION : CONSTRUCTION D'UN STOCKAGE COMPACT 
C             (INX,IPX) STRUCTURE AVEC LA DIAGONALE
C             VIA (IN,IP) = (XADJ, ADJNCY) DES TERMES EXTRADIAGONAUX 
C             ET LE STOCKAGE SEGMENT (ISEGIP, XA, DA)
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NPBLK        | -->|
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF, EX_SD_FABCAD => SD_FABCAD
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NPBLK,NSEGBLK
      INTEGER, INTENT(IN)             :: IN(NPBLK+1),IP(NSEGBLK*2+1)
      INTEGER, INTENT(IN)             :: ISEGIP(NSEGBLK*2+1)
      INTEGER, INTENT(INOUT)          :: INDTRI(NPBLK)
      INTEGER, INTENT(INOUT)          :: ISTRI(NPBLK)
      INTEGER, INTENT(INOUT)          :: INX(NPBLK+1)
      INTEGER, INTENT(INOUT)          :: IPX(NSEGBLK*2+NPBLK+1)
      DOUBLE PRECISION, INTENT(INOUT) :: ACTRI(NPBLK)
      DOUBLE PRECISION, INTENT(IN)    :: XA1(NSEGBLK),XA2(NSEGBLK)
      DOUBLE PRECISION, INTENT(IN)    :: DA(NPBLK)
      DOUBLE PRECISION, INTENT(INOUT) :: AC(NSEGBLK*2+NPBLK+1)     
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C      
      INTEGER I,J,J1,J2,JN,ISEG,ND
C       
C-----------------------------------------------------------------------
C      
C---> STOCKAGE COMPACT AVEC LA DIAGONALE : (XADJ, ADJNCY) = (INX,IPX)
C
      DO I = 1, NPBLK+1
        INX(I) = IN(I)+I-1
      ENDDO
      J2=1
      DO I = 1, NPBLK
         IPX(INX(I)) = I
         AC(INX(I)) = DA(I)
         DO J1 = INX(I), INX(I+1)-1
            JN = J1-I+1
            J = IP(JN)
            J2=J2+1
            ISEG = ISEGIP(JN)
            IPX(J2) = J
            IF(ISEG.LT.0) AC(J2) = XA1(-ISEG)
            IF(ISEG.GT.0) AC(J2) = XA2(ISEG)
         ENDDO
      ENDDO
      DO I = 1, NPBLK
         ND = INX(I+1)-INX(I)
         DO J = 1,ND
            ISTRI(J) = IPX(INX(I)+J-1)
            ACTRI(J) = AC(INX(I)+J-1)
         ENDDO
         CALL SD_STRTRI(ISTRI,ND,INDTRI)
         DO J = 1,ND
            J1 = INDTRI(J)    
            IPX(INX(I)+J-1) = ISTRI(J1)
            AC(INX(I)+J-1) = ACTRI(J1)
         ENDDO	   
      ENDDO
C
C-----------------------------------------------------------------------
C         
      RETURN
      END
