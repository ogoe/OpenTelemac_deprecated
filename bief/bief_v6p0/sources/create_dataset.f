C                       *************************
                        SUBROUTINE CREATE_DATASET
C                       *************************
C
     *(FFORMAT,NRES,TITLE,NVAR,NOMVAR,OUTVAR)
C
C***********************************************************************
C BIEF VERSION 6.0           25/11/08    R NEBAUER (LNHE) 
C***********************************************************************
C
C FONCTION : CREATE A DATA SET FOR A GIVEN FILE FORMAT IN THE FILE WITH 
C THE LOGICAL UNIT NFILE. THE TITLE OF THE DATASET IS GIVEN AS A 72 
C CHARACTER STRING.
C THE TABLE NOMVAR CONTAINS ALL POSSIBLE VARIABLES TO OUTPUT (IE THE
C NAME OF ALL VARIABLES IN THE OUTPUT BLOCK). THE LOGICAL OUTVAR
C INDICATES FOR EACH VARIABLES WHETHER IT WILL BE WRITTEN OR NOT TO THE
C DATA FILE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   FFORMAT      | -->| FILE FORMAT
C |   NRES         | -->| LOGICAL UNIT OF FILE
C |   TITLE        | -->| TITLE OF FILE
C |   NVAR         | -->| TOTAL NUMBER OF VARIABLES 
C |   NOMVAR       | -->| NAME OF VARIABLES 
C |   OUTVAR       | -->| VARIABLES TO BE PUT IN THE FILE
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE M_MED
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=8)                 , INTENT(IN) :: FFORMAT
      INTEGER                          , INTENT(IN) :: NRES
      CHARACTER(LEN=72)                , INTENT(IN) :: TITLE
      INTEGER                          , INTENT(IN) :: NVAR
      CHARACTER(LEN=32),DIMENSION(NVAR), INTENT(IN) :: NOMVAR
      LOGICAL          ,DIMENSION(NVAR), INTENT(IN) :: OUTVAR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_BEGIN('CREATE_DATASET')
!***********************************************************************
!     
      SELECT CASE (FFORMAT)
        CASE ('SERAFIN ','SERAFIND') !SERAFIN)
            CALL CREATE_DATASET_SERAFIN(
     *                          NRES,
     *                          TITLE,
     *                          NVAR,
     *                          NOMVAR,
     *                          OUTVAR)
!
        CASE ('MED     ') !MED)
            CALL CREATE_DATASET_MED(
     *                          NRES,
     *                          TITLE,
     *                          NVAR,
     *                          NOMVAR,
     *                          OUTVAR)
!
        CASE DEFAULT
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'CREATE_DATASET : MAUVAIS FORMAT : ',FFORMAT
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'CREATE_DATASET: BAD FILE FORMAT : ',FFORMAT
          ENDIF          
          CALL PLANTE(1)
          STOP
      END SELECT
!
!***********************************************************************
!     IF(DEBUG) CALL PROC_END('CREATE_DATASET')
!***********************************************************************
!
      RETURN
      END
