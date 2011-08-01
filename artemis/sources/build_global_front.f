C     ****************
      SUBROUTINE BUILD_GLOBAL_FRONT(mesh)
C     ****************
C
C***********************************************************************
C
C  ARTEMIS VERSION 6.0   16/01/2010    C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTION: RECONSTRUIT LA FRONTIERE GLOBAL DU MAILLAGE 
C      ATTENTION A UTILISER UNIQUEMENT SI ON EST EN TRAITEMENT PARALLELE
C
C
C-----------------------------------------------------------------------
C
C APPELE PAR : ARTEMI
C
C***********************************************************************
C

      
      USE BIEF
      
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'
      INTEGER :: IER,I,POSITION
      INTEGER, ALLOCATABLE :: TEMP1(:)
      INTEGER, ALLOCATABLE :: TEMP2(:)
      INTEGER :: NPOIN_MAX
      INTEGER :: NPTFR_TOT
      INTEGER P_ISUM
      COMMON/INFO/LNG,LU 
      INTEGER LNG,LU
      TYPE(BIEF_MESH), INTENT(INOUT) :: MESH
      IF(NCSIZE.LE.1) THEN
         IF(LNG.EQ.1) WRITE(LU,2018)
         IF(LNG.EQ.2) WRITE(LU,2019)
 2018    FORMAT(1X,'BUILD_GLOBAL_FRONT,',/,1X,
     *        'A UTILISER SEULEMENT EN MODE PARALLELE',///)
 2019    FORMAT(1X,'BULID_GLOBAL_FRONT,',/,1X,
     *        'USE ONLY IN THE PARALLEL VERSION',///)
         CALL PLANTE(1)
         STOP
      ENDIF             
      NPOIN_MAX=P_ISUM(MESH%NPOIN)
      write(35,*) 'NPOIN', MESH%NPOIN, NPOIN_MAX
      stop
      ALLOCATE(TEMP1(NPOIN_MAX),STAT=ier)
      IF (IER .NE. 0) stop 'erreur'
      ALLOCATE(TEMP2(NPOIN_MAX))
      IF (IER .NE. 0) stop 'erreur'
      TEMP1(:)=0
      TEMP2(:)=0
      DO I=1,MESH%NPTFR
         POSITION=MESH%KNOLG%I(MESH%NBOR%I(I))
         TEMP1(POSITION)=1
      END DO 
      CALL MPI_ALLREDUCE(TEMP1,TEMP2,NPOIN_MAX,MPI_INTEGER,
     c     MPI_MAX,
     c     MPI_COMM_WORLD,IER)

      NPTFR_TOT=0
      DO I=1,NPOIN_MAX
         IF (TEMP2(I) .NE. 0)  THEN
            NPTFR_TOT=NPTFR_TOT+1
         END IF
      END DO
c$$$      MESH%NPTFR=NPTFR_TOT
c$$$      DEALLOCATE(MESH%NBOR%I)
c$$$      ALLOCATE(MESH%NBOR%I(NPTFR)
c$$$      NBOR%I(:)=0
c$$$      DO I=1,MESH%NPTFR
c$$$         IF (TEMP2(I) .NE. 0) THEN
c$$$              MESH%NPTFR=NPTFR_TOT
c$$$
c$$$
c$$$
c$$$      NPTFR_TOT=0
      
      WRITE(*,*) 'NPTFR_TOT',NPTFR_TOT
      DEALLOCATE(TEMP1)
      DEALLOCATE(TEMP2)
      END
