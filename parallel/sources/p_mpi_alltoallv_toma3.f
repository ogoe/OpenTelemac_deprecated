C                       *********************************
                        SUBROUTINE  P_MPI_ALLTOALLV_TOMA3
C                       *********************************
C
     *(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09                 C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:   CALLING MPI_ALLTOALLV
C      ==========
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE FONCTION_TYPE
        SEQUENCE 
        INTEGER :: MYPID ! PARTITION OF THE TRACEBACK ORIGIN (HEAD)
        INTEGER :: NEPID ! THE NEIGHBOUR PARTITION THE TRACEBACK ENTERS TO 
        INTEGER :: INE   ! THE LOCAL 2D ELEMENT NR THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION   
        INTEGER :: KNE   ! THE LOCAL LEVEL THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION   
        INTEGER :: IOR   ! THE POSITION OF THE TRAJECTORY -HEAD- IN MYPID [THE 2D/3D NODE OF ORIGIN]
        INTEGER :: ISP,NSP ! NUMBERS OF RUNGE-KUTTA PASSED AS COLLECTED AND TO FOLLOW AT ALL
        DOUBLE PRECISION :: XP,YP,ZP                ! THE (X,Y,Z)-POSITION NOW 
        DOUBLE PRECISION :: SHP1,SHP2,SHP3,SHZ
        DOUBLE PRECISION :: BP
        DOUBLE PRECISION :: F(6) ! FUNCTION VALUES AT THE 6 POINT OF THE PRISME  
      END TYPE FONCTION_TYPE
C
      TYPE(FONCTION_TYPE) ::  I1(*), I5(*)
      INTEGER, INTENT(IN) ::  I2(*),I3(*),I4,I6(*),I7(*)
      INTEGER, INTENT(IN) :: I8,I9,I10
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_ALLTOALLV(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10)
C
      IF(I10.NE.0) THEN
        WRITE(LU,*) 'P_MPI_ALLTOALLV:'
        WRITE(LU,*) 'MPI ERROR ',I10
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
      


 
     
 
