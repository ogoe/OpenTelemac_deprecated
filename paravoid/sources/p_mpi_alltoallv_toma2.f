C                       *********************************
                        SUBROUTINE  P_MPI_ALLTOALLV_TOMA2
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
      TYPE CHARAC_TYPE_4D
        SEQUENCE 
        INTEGER :: MYPID ! PARTITION OF THE TRACEBACK ORIGIN (HEAD)
        INTEGER :: NEPID ! THE NEIGHBOUR PARTITION THE TRACEBACK ENTERS TO 
        INTEGER :: INE   ! THE LOCAL 2D ELEMENT NR THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION   
        INTEGER :: KNE   ! THE LOCAL LEVEL THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION   
        INTEGER :: FNE   ! THE LOCAL FREQUENCE LEVEL THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION   
        INTEGER :: IOR   ! THE POSITION OF THE TRAJECTORY -HEAD- IN MYPID [THE 2D/3D NODE OF ORIGIN]
        INTEGER :: ISP,NSP ! NUMBERS OF RUNGE-KUTTA PASSED AS COLLECTED AND TO FOLLOW AT ALL
        DOUBLE PRECISION :: XP,YP,ZP,FP                ! THE (X,Y,Z)-POSITION NOW 
        DOUBLE PRECISION :: DX,DY,DZ,DF                ! THE (X,Y,Z)-POSITION NOW 
        DOUBLE PRECISION :: BASKET(10) ! VARIABLES INTERPOLATED AT THE FOOT  
      END TYPE CHARAC_TYPE_4D
C      
      TYPE(CHARAC_TYPE_4D) ::  I1(*), I5(*)
      INTEGER, INTENT(IN) ::  I2(*),I3(*),I4,I6(*),I7(*)
      INTEGER, INTENT(IN) :: I8,I9,I10
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C    
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE  P_MPI_ALLTOALLV  VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_MPI_ALLTOALLV VOID VERSION'
C
C-----------------------------------------------------------------------
C
      RETURN
      END
      


 
     
 
