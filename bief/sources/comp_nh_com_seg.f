C                       **************************
                        SUBROUTINE COMP_NH_COM_SEG
C                       **************************
C
     *(ELTSEG,NELEM,NH_COM_SEG,DIM1NHCOM,NB_NEIGHB_SEG,NB_NEIGHB_PT_SEG,
     * GLOSEG,DIMGLO,KNOLG,NPOIN)
C
C***********************************************************************
C BIEF VERSION 6.0      22/07/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT EDF 2009                 
C 
C  22/07/2009 MODIFIED BY C. DENIS 
C  THIS VERSION AVOIDS THE (TOO) LARGE INTEGERS OF THE 5.9 VERSION                                     
C***********************************************************************
C
C    FONCTION : COMPLETING THE REAL ADDRESS OF SEGMENTS IN NH_COM_SEG
C               SEE PARINI WHERE NH_COM_SEG IS INITIALISED AT -999999
C               AND THEN FILLED WITH 4*IELEM+IFACE TO STORE IELEM AND
C               IFACE.
C
C               THEN THE ADDRESSES ARE ORDERED WITH RESPECT TO THE
C               GLOBAL NUMBER OF THE FIRST AND SECOND POINT OF 
C               EVERY SEGMENT, SO THAT THE PROCESSORS SHARE THE
C               INFORMATION ON THE SAME SEGMENTS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  ELTSEG        | -->| GIVES THE SEGMENT NUMBER OF EDGES OF ELEMENTS
C |  NELEM         | -->| NUMBER OF ELEMENTS
C |  NH_COM_SEG    | -->| ADDRESSES OF INTERFACE SEGMENTS
C |  DIM1NHCOM     | -->| FIRST DIMENSION OF NH_COM_SEG
C |  NB_NEIGHB_SEG | -->| NUMBER OF NEIGHBOUR PROCESSOR (FOR SEGMENTS)
C |NB_NEIGHB_PT_SEG| -->| NUMBER OF INTERFACE SEGMENTS FOR EVERY 
C |                |    | NEIGHBOUR PROCESSOR
C |  GLOSEG        | -->| GLOBAL NUMBERS (IN SUB-DOMAIN) OF POINTS 
C |                |    | OF A SEGMENT
C |  DIMGLO        | -->| FIRST DIMENSION OF GLOSEG
C |  KNOLG         | -->| GLOBAL NUMBERS (WHOLE MESH) FUNCTION OF 
C |                |    | LOCAL NUMBERS OF POINTS
C |  NPOIN         | -->| NUMBER OF POINTS
C |________________|____|______________________________________________|
C  MODE: -->(DONNEE NON MODIFIEE),<--(RESULTAT),<-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : INBIEF
C
C SOUS-PROGRAMME APPELE :
C
C***********************************************************************
C     
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)    :: NELEM,DIM1NHCOM,NB_NEIGHB_SEG,DIMGLO
      INTEGER, INTENT(IN)    :: NPOIN
      INTEGER, INTENT(INOUT) :: NH_COM_SEG(DIM1NHCOM,NB_NEIGHB_SEG)
      INTEGER, INTENT(IN)    :: ELTSEG(NELEM,3),GLOSEG(DIMGLO,2)
      INTEGER, INTENT(IN)    :: NB_NEIGHB_PT_SEG(NB_NEIGHB_SEG)
      INTEGER, INTENT(IN)    :: KNOLG(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IELEM,IFACE,ISEG,IPROC,IKA,I,J,B,NUMSEG
      INTEGER I11,I12,I21,I22
      LOGICAL IS_LE_THAN
C     
C-----------------------------------------------------------------------
C
      DO IPROC=1,NB_NEIGHB_SEG
        IKA = NB_NEIGHB_PT_SEG(IPROC)   
        DO ISEG=1,IKA
          IFACE=MOD(NH_COM_SEG(ISEG,IPROC),4)
          IELEM=(NH_COM_SEG(ISEG,IPROC)-IFACE)/4
          NUMSEG=ELTSEG(IELEM,IFACE)
          NH_COM_SEG(ISEG,IPROC)=NUMSEG
        ENDDO 
        IF(IKA.GT.1) THEN
          DO J=2,IKA
            B=NH_COM_SEG(J,IPROC)
            DO I=J-1,1,-1
              NUMSEG=NH_COM_SEG(I,IPROC)
              I11=KNOLG(GLOSEG(NUMSEG,1))
              I12=KNOLG(GLOSEG(NUMSEG,2))
              I21=KNOLG(GLOSEG(B     ,1))
              I22=KNOLG(GLOSEG(B     ,2))
              IF(I11.GT.I21) then
                IS_LE_THAN=.FALSE.
              ELSEIF(I11.LT.I21) THEN
                IS_LE_THAN=.TRUE.
              ELSEIF(I11.EQ.I21.AND.I12.GT.I22) THEN
                IS_LE_THAN=.FALSE.
              ELSEIF(I11.EQ.I21.AND.I12.LT.I22) THEN
                IS_LE_THAN=.TRUE.
              ELSEIF(I11.EQ.I21.AND.I12.EQ.I22) THEN
                IS_LE_THAN=.TRUE.
              ELSE
                WRITE(LU,*) 'UNEXPECTED CASE IN COMP_NH_COM_SEG'
                CALL PLANTE(1)
                STOP
              ENDIF
              IF(IS_LE_THAN) GO TO 10
              NH_COM_SEG(I+1,IPROC)=NH_COM_SEG(I,IPROC)
            ENDDO
 10         CONTINUE
            NH_COM_SEG(I+1,IPROC)=B
          ENDDO
        ENDIF              
      ENDDO
C
C-----------------------------------------------------------------------
C
      RETURN     
      END
