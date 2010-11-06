C                       *****************************
                        SUBROUTINE MULT_INTERFACE_SEG
C                       *****************************
C
     *(FSEG,NH_COM_SEG,DIM1NHCOM,NB_NEIGHB_SEG,
     * NB_NEIGHB_PT_SEG,XMUL,NSEG)
C
C***********************************************************************
C BIEF VERSION 5.9      27/02/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C COPYRIGHT EDF 2009                 
C                                     
C***********************************************************************
C
C    FONCTION : MULTIPLYING BY A CONSTANT THE INTERFACE VALUES OF A
C               FUNCTION DEFINED ON SEGMENTS
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  FSEG          |<-->| THE FUNCTION DEFINED ON SEGMENTS
C |  NH_COM_SEG    | -->| ADDRESSES OF INTERFACE SEGMENTS
C |  DIM1NHCOM     | -->| FIRST DIMENSION OF NH_COM_SEG
C |  NB_NEIGHB_SEG | -->| NUMBER OF NEIGHBOUR PROCESSOR (FOR SEGMENTS)
C |NB_NEIGHB_PT_SEG| -->| NUMBER OF INTERFACE SEGMENTS FOR EVERY 
C |                |    | NEIGHBOUR PROCESSOR
C |  XMUL          | -->| THE CONSTANT
C |  NSEG          | -->| NUMBER OF SEGMENTS
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
      INTEGER, INTENT(IN)    :: DIM1NHCOM,NB_NEIGHB_SEG,NSEG
      INTEGER, INTENT(INOUT) :: NH_COM_SEG(DIM1NHCOM,NB_NEIGHB_SEG)
      INTEGER, INTENT(IN)    :: NB_NEIGHB_PT_SEG(NB_NEIGHB_SEG)
      DOUBLE PRECISION, INTENT(INOUT) :: FSEG(NSEG),XMUL
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER ISEG,IPROC,IKA,IADSEG
C
C-----------------------------------------------------------------------
C
C     DONE ONLY IF THERE IS AT LEAST ONE OTHER SUB-DOMAIN SHARING
C     A SEGMENT WITH THIS ONE
C
      IF(NB_NEIGHB_SEG.GT.0) THEN
C
C     LOOP ON ALL NEIGHBOURING SUB-DOMAINS
C
      DO IPROC=1,NB_NEIGHB_SEG
        IKA = NB_NEIGHB_PT_SEG(IPROC)
C
C       LOOP ON ALL SEGMENTS SHARED WITH THIS SUB-DOMAIN
C       WHICH CANNOT BE SHARED WITH ANOTHER SUB-DOMAIN (UNLIKE POINTS)
C
        DO ISEG=1,IKA
C         ADDRESS IN SEGMENT NUMBERING
          IADSEG=NH_COM_SEG(ISEG,IPROC)
          FSEG(IADSEG)=FSEG(IADSEG)*XMUL
        ENDDO 
      ENDDO
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
