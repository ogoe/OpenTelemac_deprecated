C                       *****************
                        SUBROUTINE PARACO
C                       *****************
C
     *(V1,V2,V3,NPOIN,ICOM,IAN,NPLAN,NB_NEIGHB,NB_NEIGHB_PT,LIST_SEND,
     * NH_COM,DIMNHCOM,BUF_SEND,BUF_RECV,DIMBUF)
C
C***********************************************************************
C BIEF VERSION 5.9         18/07/08                      P. VEZOLLE(IBM)
C***********************************************************************
C
C   FONCTION : ASSEMBLING DATA SHARED BY SEVERAL PROCESSORS
C
C
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C | V1,V2,V3       |<-->| VECTORS TO BE COMPLETED
C | NPOIN          | -->| FIRST DIMENSION OF V1,V2,V3
C | ICOM           | -->| OPTION OF COMMUNICATION :
C |                |    | = 1 : VALUE WITH MAXIMUM ABSOLUTE VALUE
C |                |    | = 2 : CONTRIBUTIONS ADDED
C |                |    | = 3 : MAXIMUM CONTRIBUTION RETAINED
C |                |    | = 4 : MINIMUM CONTRIBUTION RETAINED
C | IAN            | -->| NUMBER OF VECTORS TO BE CONDIDERED (1, 2 OR 3)
C | NPLAN          | -->| SECOND DIMENSION OF V1,V2,V3
C | NB_NEIGHB      | -->| NUMBER OF NEIGHBOURING SUB-DOMAINS
C | NB_NEIGHB_PT   | -->| NUMBER OF POINTS SHARED WITH A SUB-DOMAIN
C | LIST_SEND      | -->| LIST OF PROCESSORS NUMBERS
C | NH_COM         | -->| NH_COM(I,IL) : GLOBAL NUMBER IN THIS
C |                |    | SUB-DOMAIN OF THE POINT NUMBER I IN THE LIST
C |                |    | OF POINTS SHARED WITH PROCESSOR NUMBER IL
C |                |    | WHOSE REAL NUMBER IS LIST_SEND(IL)
C | DIM_NH_COM     |    | FIRST DIMENSION OF NH_COM
C | BUF_SEND       |<-->| BUFFER FOR SENDING DATA
C | BUF_RECV       |<-->| BUFFER FOR RECEIVING DATA
C | DIMBUF         | -->| FIRST DIMENSION OF BUFFERS
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :  PLANTE
C
C***********************************************************************
C
      USE BIEF_DEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: NPOIN,ICOM,IAN,NPLAN,NB_NEIGHB
      INTEGER, INTENT(IN) :: DIMNHCOM,DIMBUF
      INTEGER, INTENT(IN) :: NB_NEIGHB_PT(NB_NEIGHB)
      INTEGER, INTENT(IN) :: LIST_SEND(NB_NEIGHB),NH_COM(DIMNHCOM,*)
C
      DOUBLE PRECISION, INTENT(INOUT) :: BUF_SEND(DIMBUF,*)
      DOUBLE PRECISION, INTENT(INOUT) :: BUF_RECV(DIMBUF,*)
      DOUBLE PRECISION, INTENT(INOUT) :: V1(NPOIN,NPLAN)
      DOUBLE PRECISION, INTENT(INOUT) :: V2(NPOIN,NPLAN)
      DOUBLE PRECISION, INTENT(INOUT) :: V3(NPOIN,NPLAN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IKA,IL,II,I,J,K,IPA
C      
      INTRINSIC ABS
C
      INTEGER SEND_REQ(100),RECV_REQ(100)
      INTEGER PARACO_MSG_TAG
      DATA PARACO_MSG_TAG/5000/
C
      SAVE
C
C----------------------------------------------------------------------
C
      IF(IAN.NE.1.AND.IAN.NE.2.AND.IAN.NE.3) THEN
        WRITE(LU,*) 'FALSCHE FREIWERTZAHL BEI KOMMUNIKATION',IAN,
     *              ' AUF PROZESSOR',IPID
        CALL PLANTE(1)
        STOP
      ENDIF
C
C     MESSAGE TAG UPDATE
C
      IF(PARACO_MSG_TAG.LT.1000000) THEN
        PARACO_MSG_TAG = PARACO_MSG_TAG + 1
      ELSE
        PARACO_MSG_TAG = 5001
      ENDIF
!
!== RECEIVE STEP
!
      DO IL=1,NB_NEIGHB
        IKA = NB_NEIGHB_PT(IL)
        IPA = LIST_SEND(IL)
        CALL P_IREAD(BUF_RECV(1,IL),IAN*IKA*NPLAN*8,
     *               IPA,PARACO_MSG_TAG,RECV_REQ(IL))
      ENDDO
! 
!== SEND STEP
!
      DO IL=1,NB_NEIGHB
        IKA = NB_NEIGHB_PT(IL)
        IPA = LIST_SEND(IL)
!
!** INITIALIZATION DES TABLEAUX DE COMMUNICATION
!
       K = 1
       IF(IAN.EQ.3) THEN
            DO J=1,NPLAN
              DO I=1,IKA
                II=NH_COM(I,IL)
                BUF_SEND(K,IL)  =V1(II,J)
                BUF_SEND(K+1,IL)=V2(II,J)
                BUF_SEND(K+2,IL)=V3(II,J)
                K=K+3
              ENDDO
            ENDDO
       ELSEIF(IAN.EQ.2) THEN
            DO J=1,NPLAN
              DO I=1,IKA
                II=NH_COM(I,IL)
                BUF_SEND(K,IL)  =V1(II,J)
                BUF_SEND(K+1,IL)=V2(II,J)
                K=K+2
              ENDDO
            ENDDO
       ELSEIF(IAN.EQ.1) THEN
            DO J=1,NPLAN
              DO I=1,IKA
                II=NH_COM(I,IL)
                BUF_SEND(K,IL)  =V1(II,J)
                K=K+1
              ENDDO
            ENDDO
       ENDIF
!
       CALL P_IWRIT(BUF_SEND(1,IL),IAN*IKA*NPLAN*8,
     *              IPA,PARACO_MSG_TAG,SEND_REQ(IL))
!
      ENDDO
!
!== WAIT RECEIVED MESSAGES (POSSIBILITE DE RECOUVRIR)
!
      DO IL=1,NB_NEIGHB
       IKA = NB_NEIGHB_PT(IL)
       IPA = LIST_SEND(IL)
       CALL P_WAIT_PARACO(RECV_REQ(IL),1)
!
       K=1
!
       IF(ICOM.EQ.1) THEN
            IF(IAN.EQ.3) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(ABS(BUF_RECV(K,IL)).GT.ABS(V1(II,J)))
     *                            V1(II,J)=BUF_RECV(K  ,IL)
                  IF(ABS(BUF_RECV(K+1,IL)).GT.ABS(V2(II,J)))
     *                            V2(II,J)=BUF_RECV(K+1,IL)
                  IF(ABS(BUF_RECV(K+2,IL)).GT.ABS(V3(II,J)))
     *                            V3(II,J)=BUF_RECV(K+2,IL)
                  K=K+3
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.2) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(ABS(BUF_RECV(K,IL)).GT.ABS(V1(II,J)))
     *                            V1(II,J)=BUF_RECV(K  ,IL)
                  IF(ABS(BUF_RECV(K+1,IL)).GT.ABS(V2(II,J)))
     *                            V2(II,J)=BUF_RECV(K+1,IL)
                  K=K+2
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.1) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(ABS(BUF_RECV(K,IL)).GT.ABS(V1(II,J)))
     *                            V1(II,J)=BUF_RECV(K  ,IL)
                  K=K+1
                ENDDO
              ENDDO
            ENDIF
       ELSEIF(ICOM.EQ.2) THEN
            IF(IAN.EQ.3) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  V1(II,J)=V1(II,J)+BUF_RECV(K  ,IL)
                  V2(II,J)=V2(II,J)+BUF_RECV(K+1,IL)
                  V3(II,J)=V3(II,J)+BUF_RECV(K+2,IL)
                  K=K+3
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.2) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  V1(II,J)=V1(II,J)+BUF_RECV(K  ,IL)
                  V2(II,J)=V2(II,J)+BUF_RECV(K+1,IL)
                  K=K+2
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.1) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  V1(II,J)=V1(II,J)+BUF_RECV(K  ,IL)
                  K=K+1
                ENDDO
              ENDDO
            ENDIF
       ELSEIF(ICOM.EQ.3) THEN
            IF(IAN.EQ.3) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(BUF_RECV(K  ,IL).GT.V1(II,J)) 
     *              V1(II,J)=BUF_RECV(K  ,IL)
                  IF(BUF_RECV(K+1,IL).GT.V2(II,J))
     *              V2(II,J)=BUF_RECV(K+1,IL)
                  IF(BUF_RECV(K+2,IL).GT.V3(II,J)) 
     *              V3(II,J)=BUF_RECV(K+2,IL)
                  K=K+3
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.2) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(BUF_RECV(K  ,IL).GT.V1(II,J)) 
     *              V1(II,J)=BUF_RECV(K  ,IL)
                  IF(BUF_RECV(K+1,IL).GT.V2(II,J))
     *              V2(II,J)=BUF_RECV(K+1,IL)
                  K=K+2
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.1) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(BUF_RECV(K  ,IL).GT.V1(II,J)) 
     *              V1(II,J)=BUF_RECV(K  ,IL)
                  K=K+1
                ENDDO
              ENDDO
            ENDIF
        ELSEIF(ICOM.EQ.4) THEN
            IF(IAN.EQ.3) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(BUF_RECV(K  ,IL).LT.V1(II,J)) 
     *              V1(II,J)=BUF_RECV(K  ,IL)
                  IF(BUF_RECV(K+1,IL).LT.V2(II,J))
     *              V2(II,J)=BUF_RECV(K+1,IL)
                  IF(BUF_RECV(K+2,IL).LT.V3(II,J)) 
     *              V3(II,J)=BUF_RECV(K+2,IL)
                  K=K+3
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.2) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(BUF_RECV(K  ,IL).LT.V1(II,J)) 
     *              V1(II,J)=BUF_RECV(K  ,IL)
                  IF(BUF_RECV(K+1,IL).LT.V2(II,J))
     *              V2(II,J)=BUF_RECV(K+1,IL)
                  K=K+2
                ENDDO
              ENDDO
            ELSEIF(IAN.EQ.1) THEN
              DO J=1,NPLAN
                DO I=1,IKA
                  II=NH_COM(I,IL)
                  IF(BUF_RECV(K  ,IL).LT.V1(II,J)) 
     *              V1(II,J)=BUF_RECV(K  ,IL)
                  K=K+1
                ENDDO
              ENDDO
            ENDIF
        ENDIF
C
      ENDDO
!
!== WAIT SENT MESSAGES
!
      CALL P_WAIT_PARACO(SEND_REQ,NB_NEIGHB)
C
C-----------------------------------------------------------------------
C
      RETURN
      END
