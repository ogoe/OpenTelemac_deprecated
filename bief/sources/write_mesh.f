C                       *********************
                        SUBROUTINE WRITE_MESH
C                       *********************
C
     *(FFORMAT,NFILE,MESH,NPLAN,DATE,TIME,I_ORIG,J_ORIG)
C
C***********************************************************************
C BIEF VERSION 6.0           25/11/08    R NEBAUER (LNHE) 
C***********************************************************************
C
C FUNCTION : WRITE THE MESH, DESCRIBED BY THE BIEF_MESH STRUCTURE INTO 
C            THE FILE. BIEF_MESH STRUCTURE CONTAINS INFORMATIONS ABOUT
C            CONNECTIVITY, COORDINATES, BOUNDARY NODES. OTHER 
C            INFORMATIONS NEEDED : THE DATE AND TIME INFORMATION, AND
C            THE ORIGIN OF THE COORDINATE SYSTEM (X_ORIG,Y_ORIG).
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________
C |   FFORMAT      | -->| FILE FORMAT
C |   NFILE        | -->| LOGICAL UNIT OF FILE
C |   NCSIZE       | -->| NUMBER OF PROCESSORS
C |   NPTIR        | -->| NUMBER OF INTERFACE POINTS
C |   MESH         | -->| MESH STRUCTURE
C |   NPLAN        | -->| NUMBER OF PLANES (3D)
C |________________|____|______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES : RIEN EN STANDARD
C
C***********************************************************************
C
      USE BIEF
      USE M_MED
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CHARACTER(LEN=8) ,     INTENT(IN) :: FFORMAT
      INTEGER          ,     INTENT(IN) :: NFILE,NPLAN
      TYPE(BIEF_MESH),       INTENT(IN) :: MESH
      INTEGER, DIMENSION(3), INTENT(IN) :: DATE
      INTEGER, DIMENSION(3), INTENT(IN) :: TIME
      INTEGER,               INTENT(IN) :: I_ORIG,J_ORIG
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C

!***********************************************************************
!     IF(DEBUG) CALL PROC_BEGIN('WRITE_MESH')
!***********************************************************************

      SELECT CASE (FFORMAT)
        CASE ('SERAFIN','SERAFIND')
           CALL WRITE_MESH_SERAFIN(NFILE,
     *                             MESH,
     *                             NPLAN,
     *                             DATE,
     *                             TIME,
     *                             I_ORIG,J_ORIG,
     *                             FFORMAT)
        CASE ('MED     ') 
           CALL WRITE_MESH_MED(NFILE,MESH,DBLE(I_ORIG),DBLE(J_ORIG))
        CASE DEFAULT
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'WRITE_MESH : MAUVAIS FORMAT : ',FFORMAT
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'WRITE_MESH: BAD FILE FORMAT : ',FFORMAT
          ENDIF          
          CALL PLANTE(1)
          STOP
      END SELECT

!***********************************************************************
!     IF(DEBUG) CALL PROC_END('WRITE_MESH')
!***********************************************************************

      RETURN
      END
