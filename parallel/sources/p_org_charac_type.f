!                    ****************************
                     SUBROUTINE P_ORG_CHARAC_TYPE 
!                    **************************** 
! 
     &(NOMB,TRACE,CHARACTERISTIC)                      
!
!***********************************************************************
! PARALLEL   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    MPI TYPE FOR TYPE CHARAC_TYPE - CHARACTERISTICS /
!+        USED BY TOMAWAC ONLY
!
!history  C. DENIS
!+        01/07/2011
!+        V6P1
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| NOMB           |<---| NUMBER OF VARIABLES 
!| TRACE          |<---| IF .TRUE. TRACE EXECUTION
!| CHARACTERISTIC |--->| DATATYPE FOR CHARACTERISTIC 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         
      IMPLICIT NONE 
      INCLUDE 'mpif.h' 
! 
      INTEGER, PARAMETER :: MAX_BASKET_SIZE=10 
      INTEGER, INTENT(INOUT) :: NOMB 
      INTEGER, INTENT(INOUT)  :: CHARACTERISTIC 
      LOGICAL, INTENT(IN) ::TRACE 
      TYPE CHARAC_TYPE
        INTEGER :: MYPID ! PARTITION OF THE TRACEBACK ORIGIN (HEAD)
        INTEGER :: NEPID ! THE NEIGHBOUR PARTITION THE TRACEBACK ENTERS TO
        INTEGER :: INE   ! THE LOCAL 2D ELEMENT NR THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION
        INTEGER :: KNE   ! THE LOCAL LEVEL THE TRACEBACK ENTERS IN THE NEIGBOUR PARTITION
        INTEGER :: IOR   ! THE POSITION OF THE TRAJECTORY -HEAD- IN MYPID [THE 2D/3D NODE OF ORIGIN]
        INTEGER :: ISP,NSP ! NUMBERS OF RUNGE-KUTTA PASSED AS COLLECTED AND TO FOLLOW AT ALL
        DOUBLE PRECISION :: XP,YP,ZP                ! THE (X,Y,Z)-POSITION NOW
        DOUBLE PRECISION :: DX,DY,DZ                ! THE (X,Y,Z)-POSITION NOW
        DOUBLE PRECISION :: BASKET(MAX_BASKET_SIZE) ! VARIABLES INTERPOLATED AT THE FOOT
      END TYPE CHARAC_TYPE
      INTEGER, DIMENSION(15) :: CH_BLENGTH=
     &                               (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)     
!     
!     ARRAY OF DISPLACEMENTS BETWEEN BASIC COMPONENTS, HERE INITIALISED ONLY
! 
      INTEGER (KIND=MPI_ADDRESS_KIND), DIMENSION(15) :: CH_DELTA= 
     &                               (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/) 
!
!     ARRAY OF BLOCKLENGTHS OF TYPE COMPONENTS, BASKET INITIALISED TO 1 
!      
!     ARRAY OF COMPONENT TYPES IN TERMS OF THE MPI COMMUNICATION
!  
      INTEGER, DIMENSION(15) :: CH_TYPES 
      INTEGER INTEX, IBASE, IER       
      INTEGER (KIND=MPI_ADDRESS_KIND)  :: EXTENT,ILB,IUB 
      TYPE(CHARAC_TYPE) :: CH 
!          
      INTEGER LNG,LU 
      COMMON/INFO/LNG,LU 
      INTEGER I 
!
!     NOTE JMH : P_MPI_ADDRESS 2 AND 3 ARE IN PARALLEL LIBRARY 
!                THEY ALL CALL MPI_ADDRESS BUT WITH DIFFERENT 
!                DATA TYPES (THIS IS TO ENABLE COMPILING BY NAG)            
!      
      CALL P_MPI_ADDRESS (CH%MYPID,  CH_DELTA(1),  IER)
      CALL P_MPI_ADDRESS (CH%NEPID,  CH_DELTA(2),  IER)
      CALL P_MPI_ADDRESS (CH%INE,    CH_DELTA(3),  IER)
      CALL P_MPI_ADDRESS (CH%KNE,    CH_DELTA(4),  IER)
      CALL P_MPI_ADDRESS (CH%IOR,    CH_DELTA(5),  IER)
      CALL P_MPI_ADDRESS (CH%ISP,    CH_DELTA(6),  IER)
      CALL P_MPI_ADDRESS (CH%NSP,    CH_DELTA(7),  IER)
      CALL P_MPI_ADDRESS2(CH%XP,     CH_DELTA(8),  IER)
      CALL P_MPI_ADDRESS2(CH%YP,     CH_DELTA(9),  IER)
      CALL P_MPI_ADDRESS2(CH%ZP,     CH_DELTA(10), IER)
      CALL P_MPI_ADDRESS2(CH%DX,     CH_DELTA(11),  IER)
      CALL P_MPI_ADDRESS2(CH%DY,     CH_DELTA(12),  IER)
      CALL P_MPI_ADDRESS2(CH%DZ,     CH_DELTA(13), IER)
      CALL P_MPI_ADDRESS3(CH%BASKET, CH_DELTA(14), IER) ! BASKET STATIC
!
      CALL P_MPI_TYPE_GET_EXTENT(MPI_REAL8,INTEX,ILB,IER)
          ! MARKING THE END OF THE TYPE
      CH_DELTA(15) = CH_DELTA(14) + MAX_BASKET_SIZE*INTEX ! MPI_UB POSITION
      IBASE = CH_DELTA(1)
      CH_DELTA = CH_DELTA - IBASE ! RELATIVE ADDRESSES
      IF (NOMB>0.AND.NOMB<=MAX_BASKET_SIZE) THEN
          CH_BLENGTH(14) = NOMB ! CH%BASKET RANGE APPLIED FOR COMMUNICATION
      ELSE
        WRITE(LU,*) ' @STREAMLINE::ORG_CHARAC_TYPE::',
     &        ' NOMB NOT IN RANGE [1..MAX_BASKET_SIZE]'
        WRITE(LU,*) ' MAX_BASKET_SIZE, NOMB: ',MAX_BASKET_SIZE,NOMB
        CALL PLANTE(1)
        STOP
      ENDIF
      CH_TYPES(1)=MPI_INTEGER
      CH_TYPES(2)=MPI_INTEGER
      CH_TYPES(3)=MPI_INTEGER
      CH_TYPES(4)=MPI_INTEGER
      CH_TYPES(5)=MPI_INTEGER
      CH_TYPES(6)=MPI_INTEGER
      CH_TYPES(7)=MPI_INTEGER
      CH_TYPES(8)=MPI_REAL8
      CH_TYPES(9)=MPI_REAL8
      CH_TYPES(10)=MPI_REAL8
      CH_TYPES(11)=MPI_REAL8
      CH_TYPES(12)=MPI_REAL8
      CH_TYPES(13)=MPI_REAL8
      CH_TYPES(14)=MPI_REAL8
      CH_TYPES(15)=MPI_UB       ! THE TYPE UPPER BOUND MARKER
      CALL P_MPI_TYPE_CREATE_STRUCT(15,CH_BLENGTH,CH_DELTA,
     &                         CH_TYPES,CHARACTERISTIC,IER)
      CALL P_MPI_TYPE_COMMIT(CHARACTERISTIC,IER)
!     CALL P_MPI_TYPE_LB (CHARACTERISTIC, ILB, IER)
!     CALL P_MPI_TYPE_UB (CHARACTERISTIC, IUB, IER)
      CALL P_MPI_TYPE_GET_EXTENT(CHARACTERISTIC,INTEX,ILB,IER)
      IUB=INTEX+ILB
!
      IF(TRACE) THEN
        WRITE(LU,*) ' @STREAMLINE::ORG_CHARAC_TYPE:'
        WRITE(LU,*) ' MAX_BASKET_SIZE: ', MAX_BASKET_SIZE
        WRITE(LU,*) ' SIZE(CH%BASKET): ',SIZE(CH%BASKET)
        WRITE(LU,*) ' CH_DELTA: ',CH_DELTA
        WRITE(LU,*) ' CH_BLENGTH: ',CH_BLENGTH
        WRITE(LU,*) ' CH_TYPES: ',CH_TYPES
        WRITE(LU,*) ' COMMITING MPI_TYPE_CREATE_STRUCT: ',
     &                CHARACTERISTIC
        WRITE(LU,*) ' MPI_TYPE_LB, MPI_TYPE_UB: ',ILB, IUB
      ENDIF
      IF (TRACE) WRITE(LU,*) ' -> LEAVING P_ORG_CHARAC_TYPE'
!     
!----------------------------------------------------------------------
!     
      RETURN  
      END 
 
 
