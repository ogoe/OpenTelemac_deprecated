C                       *******************
                        SUBROUTINE FRICTION
C                       *******************
C
     *(NS,G,DT,UA,H,QU,QV,CF)
C
C***********************************************************************
C TELEMAC 2D VERSION 5.4               INRIA
C
C***********************************************************************
C
C     FONCTION  : COMPUTATION OF THE FRICTION TERM
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |                | -- |  
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)             :: NS
      DOUBLE PRECISION, INTENT(IN)    :: G,DT
      DOUBLE PRECISION, INTENT(IN)    :: CF(NS)
      DOUBLE PRECISION, INTENT(IN)    :: H(NS),QU(NS),QV(NS)
      DOUBLE PRECISION, INTENT(INOUT) :: UA(3,NS)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C	
      INTEGER IS
      DOUBLE PRECISION AKAP,AKAP1,STRIC2
C                                                          
C-----------------------------------------------------------------------
C
        DO IS =1,NS
C
          STRIC2=CF(IS)**2
C
C FH-FRDATA
!          IF(H(IS).LE.1.D-12.OR.UA(1,IS).LE.1.D-12)  THEN
          IF((H(IS)   .LE.1.D-12).OR.
     *       (UA(1,IS).LE.1.D-12).OR.
     *       (CF(IS)  .LE.1.D-12)    ) THEN
C FH-FRDATA
            AKAP=0.D0
          ELSE
            AKAP= G*DT*SQRT(QU(IS)**2+QV(IS)**2)/
     *           (STRIC2*H(IS)*UA(1,IS)**(4.D0/3.D0))
          ENDIF
C
          AKAP1=1.D0/(1.D0+AKAP)
          UA(2,IS) = AKAP1*UA(2,IS)
          UA(3,IS) = AKAP1*UA(3,IS)
C
        ENDDO
C                                                          
C-----------------------------------------------------------------------
C
       RETURN
       END
