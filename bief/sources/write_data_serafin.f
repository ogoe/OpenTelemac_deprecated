C                       *****************************
                        SUBROUTINE WRITE_DATA_SERAFIN
C                       *****************************
C
     *(NFIC,NVARS,TIME,TIMESTEP,OUTVAR,BVARSOR,FFORMAT,N)
C
C***********************************************************************
C BIEF VERSION 6.0           01/04/2009               R NEBAUER (LNHE) 
C***********************************************************************
C
C FONCTION : WRITES RECORDS OF RESULTS IN A SERAFIN FORMAT FILE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   NFIC         | -->| LOGICAL UNIT OF FILE
C |   NVARS        | -->| NUMBER OF VARIABLES
C |   TIME         | -->| LOGICAL UNIT OF FILE
C |   TIMESTEP     | -->| 
C |   OUTVAR       | -->| INDICATES FOR EACH VARIABLE IF WE SHOULD 
C |                |    | PRINT IT OUT OR NOT
C |   BVARSOR      | -->| BIEF_OBJ BLOCK WITH DATA VALUES
C |   FFORMAT      | -->| FILE FORMAT
C |   N            | -->| NUMBER OF VALUES (MAY BE DIFFERENT FROM
C |                |    | THE NUMBER OF DEGREES OF FREEDOM, E.G. FOR
C |                |    | QUADRATIC ELEMENTS ONLY THE LINEAR VALUES
C |                |    | ARE EXITED)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER        ,  INTENT(IN)          :: NFIC,NVARS,N,TIMESTEP
      DOUBLE PRECISION, INTENT(IN)          :: TIME
      LOGICAL, DIMENSION(NVARS), INTENT(IN) :: OUTVAR
      TYPE(BIEF_OBJ),            INTENT(IN) :: BVARSOR
      CHARACTER(LEN=8), INTENT(IN)          :: FFORMAT
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=2)               :: RF
      DOUBLE PRECISION, DIMENSION(1) :: TTIME
      INTEGER                        :: K,ISTAT
      INTEGER                        :: IBID(1)
      CHARACTER*2                    :: CBID
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_BEGIN('WRITE_DATA_SERAFIN')
!***********************************************************************
!
      IF(FFORMAT.EQ.'SERAFIND') THEN
        RF = 'R8'
      ELSE
        RF = 'R4'
      ENDIF
!
      TTIME(1) = TIME
!
      CALL ECRI2(TTIME,IBID,CBID,1,RF,NFIC,'STD',ISTAT)
!
      DO K=1,NVARS
        IF(OUTVAR(K)) THEN 
          ! EN ESPERANT QUE CA VA MARCHER COMME CA ... 
          ! PUISQUE N N'EST PAS UN ARGUMENT ..
          ! N = BVARSOR%ADR(K)%P%DIM1
!  CORRECTION JMH 21/04/2009 NO, N IS GIVEN AND MAY BE DIFFERENT
!  FROM BVARSOR%ADR(K)%P%DIM1 (QUASI-BUBBLE AND QUADRATIC ELEMENTS)
          IF(ASSOCIATED(BVARSOR%ADR(K)%P%R)) THEN
            CALL ECRI2(BVARSOR%ADR(K)%P%R,IBID,CBID,N,RF,NFIC,'STD',
     *                 ISTAT)
          ELSE
            IF(LNG.EQ.1) THEN
              WRITE(LU,*) 'WRITE_DATA_SERAFIN : VARIABLE NO : ',K
              WRITE(LU,*) '        PAS OU MAL ALLOUEE'
              WRITE(LU,*) '        OU POINTEUR NON ASSOCIE'
            ENDIF
            IF(LNG.EQ.2) THEN
              WRITE(LU,*) 'WRITE_DATA_SERAFIN: VARIABLE NO: ',K
              WRITE(LU,*) '        NOT OR NOT WELL ALLOCATED'
              WRITE(LU,*) '        OR POINTER NOT ASSOCIATED '
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_END('WRITE_DATA_SERAFIN')
!***********************************************************************
!
      RETURN
      END
