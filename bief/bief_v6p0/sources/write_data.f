C                       *********************
                        SUBROUTINE WRITE_DATA
C                       *********************
C
     *(FFORMAT,FILERES,NVARS,TIME,TIMESTEP,OUTVAR,NOMVAR,BVARSOR,N)
C
C***********************************************************************
C BIEF VERSION 6.0                      25/11/08    R NEBAUER (LNHE) 
C***********************************************************************
C
C FUNCTION : WRITE DATA VALUES ON A MESH INTO THE DATA FILE OF THE 
C            GIVEN FILE FORMAT.
C            DATA VALUES ARE STORED IN A BIEF_OBJ BLOCK (BVARSOR),
C            AND THE LOGICAL OUTVAR INDICATES FOR EACH VARIABLE IF 
C            WE SHOULD PRINT IT OUT OR NOT.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   FFORMAT      | -->| FILE FORMAT
C |   FILERES      | -->| LOGICAL UNIT OF FILE
C |   OUTVAR       | -->| VARIABLES TO BE PUT IN THE FILE
C |   NOMVAR       | -->| NAME OF VARIABLES 
C |   BVARSOR      | -->| BIEF BLOCK CONTAINING THE VARIABLES VARIABLES 
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
      USE M_MED
      USE BIEF, EX_WRITE_DATA => WRITE_DATA
!
      IMPLICIT NONE
!
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=8), INTENT(IN)          :: FFORMAT
      INTEGER,          INTENT(IN)          :: FILERES,N
      INTEGER,          INTENT(IN)          :: NVARS
      DOUBLE PRECISION, INTENT(IN)          :: TIME
      INTEGER,          INTENT(IN)          :: TIMESTEP
      CHARACTER(LEN=32),DIMENSION(NVARS), INTENT(IN) :: NOMVAR
      LOGICAL, DIMENSION(NVARS), INTENT(IN) :: OUTVAR
      TYPE(BIEF_OBJ),            INTENT(IN) :: BVARSOR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
     
!***********************************************************************
!     IF(DEBUG) CALL PROC_BEGIN('WRITE_DATA')
!***********************************************************************
      
      SELECT CASE (FFORMAT)
        CASE ('SERAFIN ','SERAFIND')
          CALL WRITE_DATA_SERAFIN(FILERES,NVARS,TIME,TIMESTEP,
     *                            OUTVAR,BVARSOR,FFORMAT,N)

        CASE ('MED     ')
          CALL WRITE_DATA_MED(FILERES,NVARS,TIME,TIMESTEP,
     *                        NOMVAR,OUTVAR,BVARSOR)

        CASE DEFAULT
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'WRITE_DATA : MAUVAIS FORMAT : ',FFORMAT
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'WRITE_DATA: BAD FILE FORMAT : ',FFORMAT
          ENDIF          
          CALL PLANTE(1)
          STOP       
      END SELECT
      
!***********************************************************************
!     IF(DEBUG) CALL PROC_END('WRITE_DATA')
!***********************************************************************
     
      RETURN
      END 
